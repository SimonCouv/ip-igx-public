lmer_adjust_list <- function(raw, covar_pos, factor_cols, form){
  # return early with list in order to avoid loss of expensive computation time in the downstream step
  
  setDT(raw)
  
  fit_model <- function(feat, raw, covar_pos, form){
    # print(covar_pos)
    # nrow(raw)
    message("PID: ", Sys.getpid())
    setkey(raw)
    model_data <- cbind(raw[, ..covar_pos],
                        structure(raw[, ..feat], names="feat"))
    
    tmp <- tryCatch.W.E(
      lmer(formula = form, model_data, 
           control = lmerControl(optimizer = "nlminbwrap"))
    )
    list(
      res=residuals(tmp$value),
      beta=tmp$value@beta,  # fixed effects
      theta=tmp$value@theta,  # random effects
      warns=tmp$warning,
      mess=tmp$message)
  }
  
  # fits <- warns <- msgs <- list()
  
  # convert to factors, efficiently:
  for (col in factor_cols)
    raw[, (col) := as.factor(raw[[col]])]
  
  # fit model per feature
  fam_adj_l <- future_map(
    .x=structure(names(raw)[-covar_pos], names=names(raw)[-covar_pos]),
    .f= fit_model,
    raw=data.table(raw), 
    covar_pos=covar_pos, 
    form=form
  )

  return(fam_adj_l)
}
