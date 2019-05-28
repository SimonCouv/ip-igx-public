min_max_scale <- function(x, na.rm=TRUE){
  (x-min(x, na.rm = na.rm))/(max(x, na.rm=na.rm)-min(x, na.rm=na.rm))
}