make_sft <- function(data, pars, common_pars){
  sft_df <- data.frame()
  for (i in 1:nrow(pars)){
    print(i)
    iter_pars <- c(list(data=data), 
                   common_pars, 
                   list(corFnc = pars$corFnc[i],
                        networkType= pars$networktype[i],
                        corOptions = pars$corOptions[[i]]))
    capture.output(sft <- do.call(pickSoftThreshold, iter_pars))
    # browser()
    sft_df <- sft$fitIndices %>%
      dplyr::mutate(networktype= pars$networktype[i],
                    corFnc=pars$corFnc[i],
                    corfnc = ifelse("method" %in% names(pars$corOptions[[i]]),
                                       pars$corOptions[[i]]$method,
                                       pars$corFnc[i]),
                    corr_like = paste0("ips_",dplyr::recode(corfnc, cor="pearson", spearman = "spear"))
      ) %>%
      rbind(sft_df,.)
  }
  sft_df
}
