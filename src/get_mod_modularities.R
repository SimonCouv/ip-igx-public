get_mod_modularities <- function(adj, TOM=NULL,data_cor=NULL, modules){
  
  # calculate modularity, based on node degrees (= 'connectivities' in WGCNA-terminology)
  adj_mod_modularities <- calc_modularity(colors=modules,
                                          adj=adj)
  
  # degree of TOM
  #cf Zhang2005, section 5: TOM-based conn outperformed standard conn with respond to significance of outcome association
  if (!is.null(TOM)){
    TOM_mod_modularities <- calc_modularity(colors=modules,
                                            adj=TOM)
  } else
    TOM_mod_modularities <- NA
  
  if (!is.null(data_cor)){
    cor_mod_modularities <- calc_corr_modularity(colors=modules,
                                                 data_cor=as.matrix(data_cor))
  } else
    cor_mod_modularities <- NA
  
  return(list(adj_mod_modularities=adj_mod_modularities, 
              TOM_mod_modularities=TOM_mod_modularities,
              cor_mod_modularities=cor_mod_modularities))
  
}