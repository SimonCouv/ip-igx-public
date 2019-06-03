topleft <- function(dat, n=6, m=n){
  if (is.null(n))
    n <- nrow(dat)
  if (is.null(m))
    m <- ncol(dat)
  
  dat[1:n,1:m]
}