make_similarity_list <- function(correlation_l, pars){

if (!("networktype" %in% names(pars)))
  stop("the selection of networktypes most be provided as column 'networktype' in pars")
  
  map(correlation_l, 
  ~similarities_from_cor(cormat=as.matrix(.x),
                         networktypes = pars$networktype))
  
}