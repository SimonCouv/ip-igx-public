calc_modularity <- function(colors, adj) {
  # assumes adj is coded with all 1's on the diagonal
  
  # browser()
  # degree of adjacency matrix
  adj_degree <- colSums(adj) - 1
  W <- sum(adj_degree)/2
  # hist(adj_degree)  # degree distribution
  
  # calculate per-module contribution to the modularity, cf Fortunato2010, p90 eq.15 and p107 eq. 36 (for weighted graphs)
  module_modularities <- c()
  # for (c in setdiff(colors, "grey")){
  for (c in unique(colors)){
    members <- which(colors==c)
    adj_degree_intern <- colSums(adj[members,members, drop=FALSE])-1
   
    Wc <- sum(adj_degree_intern)/2   # internal edges
    Sc <- sum(adj_degree[members])  # all edges of modules nodes
    module_modularities[c] <- sum(Wc/W -(Sc/(2*W))^2)
  }
  return(module_modularities)
}