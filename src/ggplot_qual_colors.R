ggplot_qual_colors <- function(n_colors){
  d <- 360/n_colors
  h <- cumsum(c(15, rep(d,n_colors - 1)))
  hcl(h = h, c = 100, l = 65)
}