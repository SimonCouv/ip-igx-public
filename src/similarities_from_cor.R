similarities_from_cor <- function(cordat, networktypes){

  if ("unsigned" %in% networktypes)
    return(abs(cordat))
  if ("signed_hybrid" %in% networktypes){
    cordat[cordat<0]<-0
    return(cordat)
  }
  if ("signed" %in% networktypes)
    return((cordat+1)/2)
  environment()[networktypes]
}