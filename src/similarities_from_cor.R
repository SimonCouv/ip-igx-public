similarities_from_cor <- function(cormat, networktypes){
  stopifnot(is.matrix(cormat))

  if ("unsigned" %in% networktypes)
    unsigned <- abs(cormat)
  if ("signed_hybrid" %in% networktypes){
    signed_hybrid <- cormat
    signed_hybrid[signed_hybrid<0]<-0
  }
  if ("signed" %in% networktypes)
    signed <- (cormat+1)/2
  as.list.environment(environment())[networktypes]
}