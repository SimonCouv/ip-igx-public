define_colors <- function(x, breaks = range(x, na.rm = TRUE), quant_colors = c("white", "darkgreen")){
  if (is.numeric(x)){
    # x_mod <- replace_na(x, max(x, na.rm = TRUE)+1)
    # colorRamp2(breaks = c(range(x, na.rm = TRUE), max(x, na.rm = TRUE)+1), colors = c("white", "red"))(x)
    colors <- circlize::colorRamp2(breaks = breaks, colors = quant_colors)
  } else {
    colors <- structure(ggplot_qual_colors(n_distinct(x)), names=as.character(unique(x)))
    if (any(x=="grey" | x=="gray")){
      colors[x=="grey"]<-"lightgrey"
    }
    # if (any(is.na(x)))
    #   colors <- c(colors, "NA"= "grey")
  }
  colors
}