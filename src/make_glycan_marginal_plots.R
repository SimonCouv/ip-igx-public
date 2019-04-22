make_glycan_marginal_plots <- function(glycans_raw, fig_path) {
  glycan_dist <- list()
  for (i in 6:ncol(glycans_raw)){
    glycan <- names(glycans_raw)[i]
    print(glycan)
    glycan_dist[[glycan]] <- histogram_n(glycans_raw, aes(x=!!ensym(glycan)))
  }
  
  
  # cowplot::plot_grid(plotlist=top_plots, ncol=min(floor(sqrt(ntop)),8))
  
  # pagination: ggforce / gridExtra
  # combine independent plots: cowplot / gridExtra
  # => paginate independent plots: gridExta
  
  # https://cran.r-project.org/web/packages/gridExtra/vignettes/arrangeGrob.html#multiple-pages-output
  # https://stackoverflow.com/questions/39736655/ggplot2-plots-over-multiple-pages/51772409#51772409
  glycan_marginal_plots <- marrangeGrob(glycan_dist, ncol=4, nrow=7)
  ggplot2::ggsave(glycan_marginal_plots, filename = fig_path, width = 10, height = 16)
}