glycans_lmer_adjust <- function(glycans_raw, form){
  
  # adjust all glycans
  glycans_fam_adj_l <- list()
  fits <- warns <- msgs <- list()
  for (i in 6:ncol(glycans_raw)){
    
    # debug
    # i <- 7
    
    igx <- names(glycans_raw)[i]
    model_data <- cbind(glycans_raw[,1:5],
                        structure(glycans_raw[, ..i], names="glycan"))%>%
      mutate_at(.vars = vars(FID,IID,Plate_NO, Sex), as.factor)
    
    tmp <- tryCatch.W.E(
      lmer(formula = form, model_data, 
           control = lmerControl(optimizer = "nlminbwrap"))
    )
    fits[[igx]] <- tmp$value
    warns[[igx]] <- tmp$warning
    msgs[[igx]] <- tmp$message
    
    glycans_fam_adj_l[[igx]] <- residuals(fits[[igx]])
  }
  
  
  # # names of residuals correspond to row index of glycans_raw (those which are not NA)
  # head(which(!is.na(glycans_raw$ENI1H5N4F0S1)))
  # all.equal(as.character(which(!is.na(glycans_raw$ENI1H5N4F0S1))),
  #           names(glycans_fam_adj_l$ENI1H5N4F0S1))
  
  glycans_fam_adj <-glycans_fam_adj_l %>%
    map(~tibble(res=.x, index=names(.x)) %>%
          dplyr::right_join(tibble(index=as.character(1:nrow(glycans_raw))),by="index") %>%
          dplyr::select(res)) %>%
    dplyr::bind_cols() %>%
    rlang::set_names(names(glycans_fam_adj_l))
  
  return(glycans_fam_adj)
}
