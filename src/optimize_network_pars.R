optimize_network_pars <- function(data=NULL, similarity_l=NULL, common_pars, pars, hline=0.80, plt_title){
  # similarity__l: named list of correlation-like matrices names are in c(bicor, pearson, spearman)
  
  # assertions
  if (!is.null(data) & !is.null(similarity_l))
    stop("only one of 'data' and 'similarity_l' should be provided.")
  if (!is.null(similarity_l)){
    stopifnot(!is.null(names(similarity_l)))
    stopifnot(all(names(similarity_l) %in% c("bicor", "spearman", "pearson")))
    if ("corFnc" %in% names(pars)){
      warning("corfnc determined based on names of similarity_l.
              Column corFnc of pars ignored.")
    }
  }
  
  # using similarity_l
  if (!is.null(similarity_l)){
    sft_df <- make_sft_from_similarity(similarity_l, common_pars)
  }
  
  # using raw data
  if(!is.null(data))
    sft_df <- make_sft(data, pars, common_pars)
  
  # create plots
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