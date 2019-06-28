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
  plots <- sf_plot(sft_df, plt_title=plt_title, hline=hline)
  
  return(list(p_rsq=plots$p_rsq, p_k=plots$p_k, sft_df=sft_df))
}