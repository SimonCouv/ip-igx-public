calc_corr_modularity <- function(colors, data_cor){
  
  cor_pos <- ifelse(data_cor>0, data_cor, 0)
  cor_neg <- ifelse(data_cor<0, -data_cor, 0)
  diag(cor_neg) <- rep(1, nrow(cor_neg))
  
  cor_pos_modularity <- calc_modularity(colors, cor_pos)
  cor_neg_modularity <- calc_modularity(colors, cor_neg)
  
  cor_pos_degree <- colSums(cor_pos) - 1
  W_pos <- sum(cor_pos_degree)/2
  cor_neg_degree <- colSums(cor_neg) - 1
  W_neg <- sum(cor_neg_degree)/2
  
  # gomez 2009, but expressed per cluster (i.e. module-wise contribution to the weighted modularity)
  # note that in the true definition, modularity is a global propertiy -> requires modules to form a partition 
  # -> calculate module-wise modularity also for "grey" module, so true total modularity can be calculated if needed
  gomez2009 <- W_pos/(W_pos+W_neg)*cor_pos_modularity - W_neg/(W_neg+W_pos)*cor_neg_modularity
  
  return(gomez2009)
}