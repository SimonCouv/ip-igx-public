sf_plot <- function(sft_df, plt_title="scale-free characteristics",hline=0.80){
  #scale free plot
  #sft_df: output of WGCNA::pickSoftThreshold, over corr_like and networktype
  p_rsq <- sft_df %>%
    dplyr::mutate(signed_rsq = -sign(slope)*SFT.R.sq) %>%
    ggplot(aes(x=Power, y= signed_rsq))+
    geom_text(aes(label=Power))+
    facet_grid(networktype~corr_like)+
    scale_y_continuous(breaks = seq(-1,1,by=0.1))+
    geom_hline(yintercept=hline, color="red", alpha=0.5)+
    labs(title = glue("{plt_title} - rsq"))+
    theme_bw()
  p_k <- sft_df %>%
    ggplot(aes(x=Power, y= mean.k.))+
    geom_text(aes(label=Power))+
    facet_grid(networktype~corr_like)+
    labs(title = glue("{plt_title} - k"))+
    # scale_y_continuous(breaks = seq(-1,1,by=0.1))+
    theme_bw()
  return(list(p_rsq=p_rsq, p_k=p_k))
}