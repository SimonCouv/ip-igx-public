wgcna_parameter_search <- function(data=NULL, correlation_l=NULL, similarity_l=NULL,create_plot=FALSE, network_pars, cut_pars, cross_iteration_pars, feat_modality=NA){
  # valid inputs: 
  #   data
  #   correlation
  #   correlation + similarity 
  
  # assertions
  if (sum(is.null(data),is.null(correlation_l))!=1)
    stop("exactly one of 'data' and 'correlation_l' should be provided.")
  if (!is.null(data) & !is.null(similarity_l))
    stop("only one of 'data' and 'similarity_l' should be provided.")
  if (!is.null(correlation_l)){
    stopifnot(!is.null(names(correlation_l)))
    stopifnot(all(names(correlation_l) %in% c("bicor", "spearman", "pearson")))
    if ("corFnc" %in% names(pars) | "corfnc" %in% names(pars)){
      warning("corfnc determined based on names of correlation_l.
              Column corFnc of pars ignored.")
    }
  }
  if ((!is.null(data)) & ncol(data)>1000)
    warning("number of columns is ", ncol(data), " strongly consider using an optimised procedure to precalculate correlations instead of supplying raw data.")
  
  # initialize variables
  merge_pars_used_l <- list()
  color_overview_l <- list()
  mod_modularities_l <- list()
  plot_l <- list()
  
  # add list arguments to function environment
  list2env(cross_iteration_pars, envir=environment())
  
  # determine required corFnc & networktype combos
  corfnc_nw_combos <- distinct(network_pars, corfnc, networktype)

  
  # calculate similarity_l from raw data if applicable, only for required combos
  if(is.null(correlation_l)){
    message("calculating correlation(-like) matrices")
    correlation_l <- list(
      pearson=WGCNA::cor(x=data, method = "pearson", use = "pairwise.complete.obs"),
      spearman=WGCNA::cor(x=data, method = "spearman", use = "pairwise.complete.obs"),
      bicor=WGCNA::bicor(x=data, use = "pairwise.complete.obs")
    )
  }
  if(is.null(similarity_l)){
    message("calculating similarity matrices")
    similarity_l <- correlation_l %>% 
      keep(., names(.) %in% corfnc_nw_combos$corfnc) %>%   #only use correlation-like functions selected in previous stage (scale-free)
      imap(., 
           function(x,y, corfnc_nw_combos){similarities_from_cor(cormat = x, networktypes = corfnc_nw_combos$networktype[corfnc_nw_combos$corfnc==y])},
           corfnc_nw_combos=corfnc_nw_combos)
  }
  
  corcol <- map(correlation_l,
                ~define_colors(as.vector(.x), breaks = c(-1,0,1), quant_colors = c("blue", "white", "red"))
  )
  

  # one plot per set of network pars
  for (i in 1:nrow(network_pars)){
    print(i)
    
    # debug
    # i <- 8
    
    merge_pars_used_l[[i]] <- list()
    
    
    iter_network_pars <- as.list(network_pars[i,])
    list2env(iter_network_pars, envir=environment())
    
    adj <- similarity_l[[corfnc]][[networktype]]^softpower

    TOM <- TOMsimilarity(adj)
    dissTOM <-1 - TOM
    
    feat_hc_TOM <- fastcluster::hclust(as.dist(dissTOM), method = method)
    # feat_hc_adj <- hclust(as.dist(diss_adj), method = method)
  
    #loop over cut parameters
    network_module_data <- list()
    for (j in 1:nrow(cut_pars)){
      iter_cut_pars <- as.list(cut_pars[j,])
      list2env(iter_cut_pars, envir=environment())
      
    # browser()
      modules <- do.call(define_modules, c(list(feat_hc_TOM=feat_hc_TOM, dissTOM=dissTOM), iter_cut_pars))
      modules$modularities <- get_mod_modularities(adj=adj,TOM=TOM,data_cor=correlation_l[[corfnc]],
                                                   modules=modules$dendrocolors$module)
      network_module_data[[j]] <- modules
      
      
      if ("mergeCriteria" %in% names(modules$dynamic_mods_TOM)){
        merge_pars_used_l[[i]][[j]] <- do.call(rbind, modules$dynamic_mods_TOM$mergeCriteria)
      } else {   #no mergeCriteria returned by cutreeHybrid if all merges are below the cut
        x <- matrix(rep(NA,4), nrow=4)
        rownames(x) <- c("maxCoreScatter", "minGap", "maxAbsCoreScatter", "minAbsGap")
        merge_pars_used_l[[i]][[j]] <- x
      }
    }
    
    modules_df <- map_dfc(network_module_data, "dendrocolors") %>% set_names(paste0("cut_iter_", 1:nrow(cut_pars)))
    color_overview_l[[i]] <- modules_df
    mod_modularities_l[[i]] <- map(network_module_data, "modularities") %>% 
      set_names(paste0("cut_iter_", 1:nrow(cut_pars))) %>%
      map(~bind_cols(.x) %>% structure(., row.names = names(.x[[1]])))
    
    if (create_plot){
      
      
      dend_TOM <- as.dendrogram(feat_hc_TOM)

      
      # TODO: deal with grey module separately
      module_hm <- HeatmapAnnotation(df=modules_df,
                                     which = "row",
                                     show_legend = FALSE, 
                                     width = unit(3,"cm"),
                                     col =  map(as_tibble(modules_df), ~define_colors(.x)),
                                     gp = gpar(col = "white", lwd = 0.5))
      
      modality_df <- data.frame(feat_modality=feat_modality)
      modality_hm <- HeatmapAnnotation(df=modality_df,
                                       which = "row",
                                       show_legend = TRUE,
                                       col=map(as_tibble(modality_df), ~define_colors(.x)),
                                       gp = gpar(col = "white", lwd = 0.5))
      
      core_hm_modality <- ComplexHeatmap::Heatmap(matrix=correlation_l[[corfnc]],
                                                  name="spearman correlation",
                                                  col= corcol[[corfnc]],
                                                  cluster_rows = dend_TOM,
                                                  cluster_columns = rev(dend_TOM),
                                                  row_names_gp = gpar(fontsize = 5),
                                                  column_dend_side = "top",
                                                  show_column_names = FALSE,
                                                  right_annotation = modality_hm)
      # browser()
      hm_list <- module_hm + core_hm_modality
      message("generating module heatmap")
      # plot_l[[i]] <- draw(hm_list)
      plot_l[[i]] <- hm_list
    }
    
  }
  # browser()
  color_overview <- do.call(cbind, color_overview_l) %>% set_names(glue("nw_{nw_range}_cut_{cut_range}", 
                                                                        nw_range=rep(1:nrow(network_pars), each=nrow(cut_pars)),
                                                                        cut_range=rep(1:nrow(cut_pars),nrow(network_pars))))
  
  
  
  # make modularity overview
  modularity_overview_nogrey <- mod_modularities_l %>%
    map(~map(.x, ~rownames_to_column(.x) %>%
               dplyr::filter(rowname != "grey") %>%
               dplyr::select(-rowname) %>%
               dplyr::summarise_all(sum) %>%
               set_names(paste0(names(.), "_no_grey"))) %>% 
          dplyr::bind_rows(.id="cut_iter")) %>%
    dplyr::bind_rows(.id="network_iter")
  modularity_overview_inclgrey <- mod_modularities_l %>%
    map(~map(.x, ~dplyr::summarise_all(.x, sum) %>%
               set_names(paste0(names(.), "_incl_grey"))) %>% 
          dplyr::bind_rows(.id="cut_iter")) %>%
    dplyr::bind_rows(.id="network_iter")
  
  modularity_overview <- inner_join(modularity_overview_nogrey, 
                                    modularity_overview_inclgrey,
                                    by=c("cut_iter", "network_iter")) %>% 
    dplyr::mutate(cut_iter = str_match(cut_iter, 
                                       pattern = "cut_iter_([:digit:]+)")[,2]) %>%
    dplyr::arrange(desc(adj_mod_modularities_no_grey)) %>%
    dplyr::left_join(rownames_to_column(network_pars),
                     by=c("network_iter" = "rowname")) %>% 
    dplyr::left_join(rownames_to_column(cut_pars), 
                     by=c("cut_iter" = "rowname"))
  
  
  return(list(merge_pars_used=merge_pars_used_l, color_overview=color_overview, 
              mod_modularities = mod_modularities_l, modularity_overview=modularity_overview, plots=plot_l))    
}






