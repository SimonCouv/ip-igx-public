lmer_adjust <- function(raw, covar_pos, factor_cols, form){
  
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
    residuals(tmp$value)
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
  
  # for (feat in  names(raw)[-covar_pos]){
  #   
  #   model_data <- cbind(raw[, ..covar_pos],
  #                       structure(raw[, ..feat], names="feat"))
  #   
  #   tmp <- tryCatch.W.E(
  #     lmer(formula = form, model_data, 
  #          control = lmerControl(optimizer = "nlminbwrap"))
  #   )
  #   fits[[feat]] <- tmp$value
  #   warns[[feat]] <- tmp$warning
  #   msgs[[feat]] <- tmp$message
  #   
  #   fam_adj_l[[feat]] <- residuals(fits[[feat]])
  # }
  # 
  
  # # names of residuals correspond to row index of raw (those which are not NA)
  # head(which(!is.na(raw$ENI1H5N4F0S1)))
  # all.equal(as.character(which(!is.na(raw$ENI1H5N4F0S1))),
  #           names(fam_adj_l$ENI1H5N4F0S1))
  
  # collect results in tibble
  fam_adj <-fam_adj_l %>%
    future_map(~tibble(res=.x, index=names(.x)) %>%
          dplyr::right_join(tibble(index=as.character(1:nrow(raw))),by="index") %>%  #make residuals complete (i.e. fill in missing lines)
          dplyr::select(res)) %>%
    dplyr::bind_cols() %>%
    rlang::set_names(names(fam_adj_l))
  
  return(fam_adj)
}
