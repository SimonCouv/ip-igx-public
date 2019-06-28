wgcna_parameter_search <- function(data=NULL, corr_like=NULL, similarity=NULL,create_plot=FALSE, corfnc=NULL, networktype=NULL, method = "average", powers, feat_modality=NA, cut_pars){
  # valid inputs: 
  #   A. data
  #   B. corr_like
  #   C. corr_like + similarity 
  cat_time_pid <- function(message){
    cat(as.character(Sys.time()), "- PID", Sys.getpid(), " -  ", message, "\n")
  }
  
  cat_time_pid("checking assertions")
  # assertions
  if (sum(is.null(data),is.null(corr_like))!=1)
    stop("exactly one of 'data' and 'corr_like' should be provided. Reason: corr_like is always required, data is only required to calc corr_like")
  if (!is.null(data) & !is.null(similarity))
    stop("only one of 'data' and 'similarity' should be provided. Reason: if similarity is available, so should corr_like. Save time by providing pre-calculated corr_like.")
  if (!is.null(corfnc) & !is.null(corr_like))
    stop("corfnc should only be specified if only data is provided")
  if (!is.null(networktype) & !is.null(similarity))
    stop("networktype should only be specified if similarity is not provided")
  if ((!is.null(data)) && ncol(data)>1000)
    warning("number of columns is ", ncol(data), " strongly consider using an optimised procedure to precalculate corr_likes instead of supplying raw data.")
  
  # initialize variables
  merge_pars_used_l <- list()
  color_overview_l <- list()
  mod_modularities_l <- list()
  plot_l <- list()
  
  # calculate similarity from raw data if applicable
  if(is.null(corr_like)){
    cat_time_pid("calculating corr_like(-like) matrices")
    corr_like <- switch (
      corfnc,
      pearson=WGCNA::cor(x=data, method = "pearson", use = "pairwise.complete.obs"),
      spear=WGCNA::cor(x=data, method = "spearman", use = "pairwise.complete.obs"),
      bicor=WGCNA::bicor(x=data, use = "pairwise.complete.obs")
    )
  }
  
  # calculate similarity if required
  if(is.null(similarity)){
    cat_time_pid("calculating similarity matrices")
    similarity <- similarity_from_cor(cordat = corr_like, networktype = networktype)
  }
  

  # one plot per power
  for (softpower in powers){  
    cat_time_pid(paste("softpower:", softpower))
    gc()
    cat_time_pid(paste("softpower:", softpower, "gc completed"))
    mem_str <- format(structure(pryr::mem_used(),class="object_size"), units="auto")
    cat_time_pid(paste("softpower:", softpower, "mem_used: ",mem_str))
    
    
    # each softpower defines a 'network', with adjacency given by:
    adj <- as.matrix(similarity)^softpower
    cat_time_pid(paste("softpower:", softpower, "adj calculated"))
    network_module_data <- list()

    TOM <- TOMsimilarity(adj)
    cat_time_pid(paste("softpower:", softpower, "TOM calculated"))
    dissTOM <-1 - TOM
    feat_hc_TOM <- fastcluster::hclust(as.dist(dissTOM), method = method)
    cat_time_pid(paste("softpower:", softpower, "hc calculated"))
    
    gc()
    cat_time_pid(paste("softpower:", softpower, "gc completed"))
    mem_str <- format(structure(pryr::mem_used(),class="object_size"), units="auto")
    cat_time_pid(paste("softpower:", softpower, "mem_used: ",mem_str))
  
    #loop over cut parameters
    for (j in 1:nrow(cut_pars)){
      mem_str <- format(structure(pryr::mem_used(),class="object_size"), units="auto")
      cat_time_pid(paste("softpower:", softpower, "cut: ",j, "mem_used: ",mem_str))
      iter_cut_pars <- as.list(cut_pars[j,])
      # list2env(iter_cut_pars, envir=environment())
      
    # browser()
      modules <- do.call(define_modules, c(list(feat_hc_TOM=feat_hc_TOM, dissTOM=dissTOM), iter_cut_pars))
      modules$modularities <- get_mod_modularities(adj=adj,TOM=TOM,data_cor=corr_like,
                                                   modules=modules$dendrocolors$module)
      network_module_data[[j]] <- modules
      
      
      if ("mergeCriteria" %in% names(modules$dynamic_mods_TOM)){
        merge_pars_used_l[[as.character(softpower)]][[j]] <- do.call(rbind, modules$dynamic_mods_TOM$mergeCriteria)
      } else {   #no mergeCriteria returned by cutreeHybrid if all merges are below the cut
        x <- matrix(rep(NA,4), nrow=4)
        rownames(x) <- c("maxCoreScatter", "minGap", "maxAbsCoreScatter", "minAbsGap")
        merge_pars_used_l[[as.character(softpower)]][[j]] <- x
      }
    }
    mem_str <- format(structure(pryr::mem_used(),class="object_size"), units="auto")
    cat_time_pid(paste("finalised module cuts; mem_used:",mem_str))
    
    # remove large unneeded objects
    rm(list=c("adj", "TOM"))
    mem_str <- format(structure(pryr::mem_used(),class="object_size"), units="auto")
    cat_time_pid(paste("finalised module cuts; mem_used:",mem_str))
    
    # gc
    gc()
    cat_time_pid(paste("softpower:", softpower, "gc completed"))
    mem_str <- format(structure(pryr::mem_used(),class="object_size"), units="auto")
    cat_time_pid(paste("softpower:", softpower, "mem_used: ",mem_str))
    
    #module df, but in single assignment step
    color_overview_l[[as.character(softpower)]] <- map_dfc(network_module_data, "dendrocolors") %>% set_names(paste0(1:nrow(cut_pars)))
    mem_str <- format(structure(pryr::mem_used(),class="object_size"), units="auto")
    cat_time_pid(paste("finalised building and assigning module_df; mem_used:",mem_str))

    # mod_modularities
    mod_modularities_l[[as.character(softpower)]] <- map(network_module_data, "modularities") %>% 
      set_names(paste0(1:nrow(cut_pars))) %>%
      map(~bind_cols(.x) %>% mutate(colors=names(.x[[1]])))
    mem_str <- format(structure(pryr::mem_used(),class="object_size"), units="auto")
    cat_time_pid(paste("finalised assigning mod_modularities; mem_used:",mem_str))
    
    # plots
    if (create_plot){
      cat_time_pid("generating module heatmap")
      plot_l[[as.character(softpower)]] <- make_heatmap_list(
        corr_like = corr_like,
        feat_hc_TOM=feat_hc_TOM,
        color_overview = color_overview_l[[as.character(softpower)]],
        modality_df = data.frame(feat_modality=feat_modality)
      )
    }
    
  }
  # browser()
  # color_overview <- do.call(cbind, color_overview_l) %>% set_names(glue("power_{power_range}_cut_{cut_range}", 
  #                                                                       power_range=rep(powers, each=nrow(cut_pars)),
  #                                                                       cut_range=rep(1:nrow(cut_pars),length(powers))))
  
  if (create_plot){
    return(list(merge_pars_used=merge_pars_used_l, color_overview=color_overview_l, 
                mod_modularities = mod_modularities_l, plots=plot_l))  
    
  } else {
  return(list(merge_pars_used=merge_pars_used_l, color_overview=color_overview_l, 
              mod_modularities = mod_modularities_l))   
  }
}






