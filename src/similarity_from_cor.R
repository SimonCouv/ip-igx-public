similarity_from_cor <- function(cordat, networktype){
  
  if ("unsigned" == networktype)
    return(abs(cordat))
  if ("signed_hybrid" == networktype){
    cordat[cordat<0]<-0
    return(cordat)
  }
  if ("signed" == networktype)
    return((cordat+1)/2)
}