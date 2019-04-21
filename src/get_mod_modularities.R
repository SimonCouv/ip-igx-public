get_mod_modularities <- function(adj, TOM=NULL, clusters){
  
  # calculate modularity, based on node degrees (= 'connectivities' in WGCNA-terminology)
  adj_mod_modularities <- calc_modularity(colors=clusters,
                                          adj=adj)
  
  # degree of TOM
  #cf Zhang2005, section 5: TOM-based conn outperformed standard conn with respond to significance of outcome association
  if (!is.null(TOM)){
    TOM_mod_modularities <- calc_modularity(colors=clusters,
                                            adj=TOM)
  } else
    TOM_mod_modularities <- NA
  
  
  # browser()
  # connectivity, after softpower-transform
  
  # igraph functions
  # power_graph <- graph_from_adjacency_matrix(adjmatrix = adj,
  #                                            add.colnames = "node_name",
  #                                            mode = "undirected",
  #                                            weighted = TRUE)
  # 
  # modularity(power_graph,
  #            membership = as.numeric(as.factor(dendrocolors$module)))
  
  return(list(adj_mod_modularities=adj_mod_modularities, TOM_mod_modularities=TOM_mod_modularities))
  
}