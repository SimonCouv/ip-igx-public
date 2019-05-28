lmer_adjust <- function(raw, covar_pos, factor_cols, form){
  
  # perform adjustments
  fam_adj_l <- lmer_adjust_list(raw, covar_pos, factor_cols, form)
  
  # collect results in tibble
  fam_adj <- fam_adj_collect(fam_adj_l, raw)
  
  return(fam_adj)
}
