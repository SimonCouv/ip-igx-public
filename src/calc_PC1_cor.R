calc_PC1_cor <- function(IPs, colors, ip_select_scores, ips_trans_fam_adj, ...){
  res <- c()
  
  # WGCNA::cor(ips_trans_fam_adj[,IPs], ip_select_scores[,colors])
  
  for (i in seq_along(IPs)){
    if (i %% 1000 ==0) print(i)
    res[i] <- cor(ip_select_scores[,colors[i]], ips_trans_fam_adj[,IPs[i]], ...)
  }
  res
}