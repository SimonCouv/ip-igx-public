make_marginal_plots <- function(data, fig_path) {
  data_dist <- list()
  for (i in seq_along(data)){
    var <- names(data)[i]
    print(var)
    data_dist[[var]] <- histogram_n(data, aes(x=!!ensym(var)))
  }
  
  
  # cowplot::plot_grid(plotlist=top_plots, ncol=min(floor(sqrt(ntop)),8))
  
  # pagination: ggforce / gridExtra
  # combine independent plots: cowplot / gridExtra
  # => paginate independent plots: gridExta
  
  # https://cran.r-project.org/web/packages/gridExtra/vignettes/arrangeGrob.html#multiple-pages-output
  # https://stackoverflow.com/questions/39736655/ggplot2-plots-over-multiple-pages/51772409#51772409
  marginal_plots <- marrangeGrob(data_dist, ncol=4, nrow=7)
  ggplot2::ggsave(marginal_plots, filename = fig_path, width = 10, height = 16)
}