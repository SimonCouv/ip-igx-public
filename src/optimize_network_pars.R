optimize_network_pars <- function(common_pars, pars, powers, hline=0.80, plt_title){
  
  sft_df <- data.frame()
  for (i in 1:nrow(pars)){
    print(i)
    iter_pars <- c(common_pars, list(corFnc = pars$corFnc[i], networkType= pars$networktype[i]))
    sft <- do.call(pickSoftThreshold, iter_pars)
    # browser()
    sft_df <- sft$fitIndices %>%
      dplyr::mutate(networktype= pars$networktype[i],
                    corFnc=pars$corFnc[i]) %>%
      rbind(sft_df,.)
  }
  
  p_rsq <- sft_df %>%
    dplyr::mutate(signed_rsq = -sign(slope)*SFT.R.sq) %>%
    ggplot(aes(x=Power, y= signed_rsq))+
    geom_text(aes(label=Power))+
    facet_grid(networktype~corFnc)+
    scale_y_continuous(breaks = seq(-1,1,by=0.1))+
    geom_hline(yintercept=hline, color="red", alpha=0.5)+
    labs(title = glue("{plt_title} - rsq"))+
    theme_bw()
  
  p_k <- sft_df %>%
    ggplot(aes(x=Power, y= mean.k.))+
    geom_text(aes(label=Power))+
    facet_grid(networktype~corFnc)+
    labs(title = glue("{plt_title} - k"))+
    # scale_y_continuous(breaks = seq(-1,1,by=0.1))+
    theme_bw()
  
  return(list(p_rsq=p_rsq, p_k=p_k, sft_df=sft_df))
}