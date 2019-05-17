wgcna_parameter_search <- function(data, network_pars, cut_pars, cross_iteration_pars, create_plot=FALSE, feat_modality=NA){
  merge_pars_used_l <- list()
  color_overview_l <- list()
  mod_modularities_l <- list()
  plot_l <- list()
  
  list2env(cross_iteration_pars, envir=environment())
  message("calculating correlation matrix")
  cormat <- cor(data, method = "spearman", use = "pairwise.complete.obs")
  corcol <- define_colors(as.vector(cormat), breaks = c(-1,0,1), quant_colors = c("blue", "white", "red"))
  
  # browser()
  
  # one plot per set of network pars
  for (i in 1:nrow(network_pars)){
    print(i)
    
    # debug
    # i <- 8
    
    merge_pars_used_l[[i]] <- list()
    
    
    iter_network_pars <- as.list(network_pars[i,])
    list2env(iter_network_pars, envir=environment())

        adj <- adjacency(datExpr = data,
                     type=networktype,
                     power = softpower,
                     corFnc = corfnc)
    diss_adj <- 1 - adj
    
    TOM <- TOMsimilarity(adj)
    dissTOM <-1 - TOM
    
    feat_hc_TOM <- hclust(as.dist(dissTOM), method = method)
    # feat_hc_adj <- hclust(as.dist(diss_adj), method = method)
  
    
    
    #loop over cut parameters
    network_module_data <- list()
    for (j in 1:nrow(cut_pars)){
      iter_cut_pars <- as.list(cut_pars[j,])
      list2env(iter_cut_pars, envir=environment())
      
    # browser()
      modules <- do.call(define_modules, c(list(feat_hc_TOM=feat_hc_TOM, dissTOM=dissTOM), iter_cut_pars))
      modules$modularities <- get_mod_modularities(adj,TOM,modules$dendrocolors$module)
      network_module_data[[j]] <- modules
      
      
      if ("mergeCriteria" %in% names(modules$dynamic_mods_TOM)){
        merge_pars_used_l[[i]][[j]] <- do.call(rbind, modules$dynamic_mods_TOM$mergeCriteria)
      } else {   #no mergeCriteria returned by cutreeHybrid if all merges are below the cut
        x <- matrix(rep(NA,4), nrow=4)
        rownames(x) <- c("maxCoreScatter", "minGap", "maxAbsCoreScatter", "minAbsGap")
        merge_pars_used_l[[i]][[j]] <- x
      }
    }
    
    # browser()
    if (create_plot){
      

      dend_TOM <- as.dendrogram(feat_hc_TOM)
      modules_df <- map_dfc(network_module_data, "dendrocolors") %>% set_names(paste0("cut_iter_", 1:nrow(cut_pars)))
      
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
      
      core_hm_modality <- ComplexHeatmap::Heatmap(matrix=cormat,
                                         name="spearman correlation",
                                         col= corcol,
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
    
    color_overview_l[[i]] <- modules_df
    mod_modularities_l[[i]] <- map(network_module_data, "modularities") %>% 
      set_names(paste0("cut_iter_", 1:nrow(cut_pars))) %>%
      map(~bind_cols(.x) %>% structure(., row.names = names(.x[[1]])))
    
  }
  # browser()
  color_overview <- do.call(cbind, color_overview_l) %>% set_names(glue("nw_{nw_range}_cut_{cut_range}", 
                                                                        nw_range=rep(1:nrow(network_pars), each=nrow(cut_pars)),
                                                                        cut_range=rep(1:nrow(cut_pars),nrow(network_pars))))
  
  # make modularity overview
  modularity_overview <- mod_modularities_l %>%
    map(~map(.x, ~rownames_to_column(.x) %>%
               dplyr::filter(rowname != "grey") %>%
               dplyr::select(-rowname) %>%
               dplyr::summarise_all(sum)) %>%
          dplyr::bind_rows(.id="cut_iter")) %>%
    dplyr::bind_rows(.id="network_iter") %>%
    dplyr::mutate(cut_iter = str_match(cut_iter, pattern = "cut_iter_([:digit:]+)")[,2]) %>%
    dplyr::arrange(desc(adj_mod_modularities)) %>%
    dplyr::left_join(rownames_to_column(network_pars), by=c("network_iter" = "rowname")) %>% 
    dplyr::left_join(rownames_to_column(cut_pars), by=c("cut_iter" = "rowname"))
  
  
  return(list(merge_pars_used=merge_pars_used_l, color_overview=color_overview, 
              mod_modularities = mod_modularities_l, plots=plot_l, modularity_overview=modularity_overview))    
}