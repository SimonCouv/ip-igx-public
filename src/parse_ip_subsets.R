parse_ip_subsets <- function(feat_data, CD_prop_cutoff=0.8, ips_trans_fam_adj, threshold_t_p=0.05){
  # CD_prop_cutoff: proportion of features of a lineage which must have a CD state for that CD state to represent that lineage within the module
  # returns: 
  #   core_CD_list: list(source=list(source_perc, core_CD for that source))
  #   CD_table_list: raw CD data for all features from all sources in the module
  
  lineage_stat <- count(feat_data, composite_lin_source) %>% 
    mutate(prop=n/sum(n))
  
  feat_data_split <- split(feat_data, f=feat_data$composite_lin_source)
  
  subset_regex <- "([^(/+\\-)]+)([+\\-]{1})"
  
  core_CD_list <- core_CD_list_equiv <-list()
  CD_table_list <- CD_table_list_equiv <- list()
  for (lineage in names(feat_data_split)){
    print(lineage)
    
    feat_data_part <- feat_data_split[[lineage]]
    
    # browser()
    
    if (lineage=="Lin"){
      message("skipping module CD subset representation for lineage (i.e. composite_lin_source) 'Lin' because not applicable")
      core_CD_list[[lineage]] <- core_CD_list_equiv[[lineage]] <- NA
      CD_table_list[[lineage]] <- CD_table_list_equiv[[lineage]] <- dplyr::select(feat_data_part, feat)
      next
    }
    if (lineage=="MFI"){
      message("skipping module CD subset representation for lineage (i.e. composite_lin_source) 'MFI' because not applicable")
      core_CD_list[[lineage]] <- core_CD_list_equiv[[lineage]] <- NA
      CD_table_list[[lineage]] <- CD_table_list_equiv[[lineage]] <- dplyr::select(feat_data_part, feat)
      next
    }
    
    # debug
    # lineage <- "B cells/MD-"
    
    
    if (!any(str_detect(feat_data_part$subset_name, pattern="([^(/+\\-)]+)([+\\-]{1})"))){
      message("skipping module CD subset representation for lineage (i.e. composite_lin_source) '", lineage, "' because none of the subset_names match the regex '", subset_regex, "'")
      core_CD_list[[lineage]] <- core_CD_list_equiv[[lineage]] <- NA
      CD_table_list[[lineage]] <- CD_table_list_equiv[[lineage]] <- dplyr::select(feat_data_part, feat)
      next
    }
    
    CD_table <- feat_data_part$subset_name %>% 
      str_match_all(pattern=subset_regex) %>%
      set_names(feat_data_part$feat) %>%
      map(data.table) %>%
      bind_rows(.id="feat") %>%
      set_names(c("feat", "match", "CD", "status")) %>%
      dplyr::select(-match) %>%
      spread(CD, status)
    
    sortvar_order <- map_dbl(dplyr::select(CD_table,-feat), ~mean(is.na(.x))) %>%
      sort() %>%
      names() %>%
      c(., "feat") %>%
      map(rlang::sym)
    
    CD_table <- arrange(CD_table, !!!sortvar_order) %>%
      dplyr::select(!!!sortvar_order)
    CD_table_list[[lineage]] <- CD_table
    
    # browser()
    max_props <- CD_table %>%
      dplyr::select(-feat) %>%
      summarise_all(
        function(x){
          prop.table(table(x, useNA = "always")) %>%
            .[!is.na(names(.))] %>%
            max(.)
        }
      )
    pos <- which(max_props > CD_prop_cutoff)
    
    max_level <- CD_table %>%  #most abundant, non-NA level per core CD
      summarise_all(
        function(x){
          prop.table(table(x, useNA = "always")) %>%
            .[!is.na(names(.))] %>%
            which.max(.) %>%
            names()
        }
      )
    
    core_CD_list[[lineage]] <- max_level[,pos, drop=FALSE]  %>%
      map2_chr(names(.), ., paste0) %>%
      paste0(., collapse = "")
    
    # browser()
    
    # calculate extended core
    # 1. replace NA by equivalent sign where appropriate
    CD_table_equiv <- CD_table
    nccds <- names(CD_table[which(max_props<CD_prop_cutoff)])
    if (length(nccds)>0){     # if there are any non-core CDs
      for (nccd in nccds){      # for each non-core CD
        print(nccd)
        
        # debug
        # nccd <- names(CD_table[which(max_props<CD_prop_cutoff)])[3]
        
        nccd_levels <-unique(CD_table[,nccd])
        
        if (sum(is.na(nccd_levels))==1 & sum(grepl(na.omit(nccd_levels),pattern = "\\+|\\-"))==1){   #if the non-core CD has two levels: NA and either "+" or "-"
          nccd_sign <- str_extract(na.omit(nccd_levels),pattern = "\\+|\\-")
          for (i in which(is.na(CD_table[,nccd]))){   # for each NA value of the non-core CD
            
            feat <- CD_table$feat[i]
            CD_line <- CD_table[i,] %>% dplyr::select(-feat)
            CD_line[,nccd] <- nccd_sign
            equiv_feats <- CD_table %>%
              dplyr::select(-feat) %>% 
              split(x=., f=seq_len(nrow(.))) %>% 
              map_lgl(~all(replaceMissing(.x,"na")==replaceMissing(CD_line, "na"))) %>% 
              CD_table$feat[.]
            
            # equiv_feats can have length>1 if duplicate features 
            # (same CD combination, different color panel) are present in module 
            # see EDA.Rmd, section 'duplicated subset_names'
            if (length(equiv_feats)>0){              #test whether the feature has an equivalent
              equivalence <- t_p <- c()
              for (equiv_feat in equiv_feats){
                t_p[equiv_feat] <- t.test(
                  x=ips_trans_fam_adj[,feat,drop=TRUE],
                  y=ips_trans_fam_adj[,equiv_feat,drop=TRUE],
                  paired = TRUE,
                  na.action="na.omit"
                )$p.value
                
                equivalence[equiv_feat] <- all(replaceMissing(ips_trans_fam_adj[,feat,drop=TRUE],"na") == replaceMissing(ips_trans_fam_adj[,equiv_feat,drop=TRUE],"na") || 
                                                 t_p[equiv_feat]  > threshold_t_p )   #because t.test pv is NA if all pairwise differences are zero
              }

              # test whether feat and the equivalent feature are pairwise highly similar, 
              # which indicates that the feat and equiv_feat measure the same subset of immune cells,
              # i.e. there was no additional discriminative value in splitting the
              # subset of immune cells by nccd
              if (all(equivalence)){     # equiv_feats can have length>1 if duplicate features, see notes 25/06/2019
                CD_table_equiv[i, nccd] <- nccd_sign
                cat(i, feat, "considered equivalent with", paste0(equiv_feats, collapse="/"), "- paired t-test p=", paste0(t_p, collapse = "/"), "\n")
              } else cat(i, feat, "NOT considered equivalent with", paste0(equiv_feats, collapse="/"), "- paired t-test p=", paste0(t_p, collapse = "/"), "\n")
            } else cat(i, feat, "No equivalent feat in module for non-core CD", nccd,"\n")
          }
        }
      }
      CD_table_list_equiv[[lineage]] <- CD_table_equiv
      
      max_props_equiv <- CD_table_equiv %>%
        dplyr::select(-feat) %>%
        summarise_all(
          function(x){
            prop.table(table(x, useNA = "always")) %>%
              .[!is.na(names(.))] %>%
              max(.)
          }
        )
      pos_equiv <- which(max_props_equiv > CD_prop_cutoff)
      max_level_equiv <- CD_table_equiv %>%  
        summarise_all(
          function(x){
            prop.table(table(x, useNA = "always")) %>%
              .[!is.na(names(.))] %>%
              which.max(.) %>%
              names()
          }
        )
      core_CD_list_equiv[[lineage]] <- max_level_equiv[pos_equiv]  %>%
        map2_chr(names(.), ., paste0) %>%
        paste0(., collapse = "")
      
    } else {  # if there are no non-core CDs for this lineage
      core_CD_list_equiv[[lineage]] <- core_CD_list[[lineage]]
      CD_table_list_equiv[[lineage]] <- CD_table_list[[lineage]]
    }
  }
  
  
  # combine results for all lineages in the module
  core_CD_df <-tibble(
    composite_lin_source=names(core_CD_list),
    core_CD=unlist(core_CD_list),
    core_CD_equiv=unlist(core_CD_list_equiv),
    CD_table = CD_table_list,
    CD_table_equiv = CD_table_list_equiv
  ) %>% 
    left_join(lineage_stat, by="composite_lin_source")
  
  return(core_CD_df = core_CD_df)
  
}



