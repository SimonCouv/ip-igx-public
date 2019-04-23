make_top_assoc_scatters <- function(ntop, ip_igx_univar_good, glycans, ips, file_path){
  
  top_plots_list <- list()
  ip_igx_univar_good <- arrange(ip_igx_univar_good, pvalue)  #repeat for safety
  
  for (i in 1:ntop){
    igx <- ip_igx_univar_good$response[i]
    ip <- ip_igx_univar_good$predictor[i]
    igx_cols <- c("IID",igx)
    ip_cols <- c("IID",ip)
    igx_data <- glycans[, ..igx_cols]
    ip_data <- ips[, ..ip_cols]
    data <- full_join(igx_data, ip_data, by="IID") %>% dplyr::select(-IID)
    # top_plots[[i]] <- ggplot(data, aes(x=!!ensym(ip), y=!!ensym(igx)))+
    #   geom_point()
    top_plots_list[[i]] <- scatter_smooth(data, aes(x=!!ensym(ip), y=!!ensym(igx)))
  }
  
  
  # cowplot::plot_grid(plotlist=top_plots, ncol=min(floor(sqrt(ntop)),8))
  
  # pagination: ggforce / gridExtra
  # combine independent plots: cowplot / gridExtra
  # => paginate independent plots: gridExta
  
  # https://cran.r-project.org/web/packages/gridExtra/vignettes/arrangeGrob.html#multiple-pages-output
  # https://stackoverflow.com/questions/39736655/ggplot2-plots-over-multiple-pages/51772409#51772409
  top_plots <- marrangeGrob(top_plots_list, ncol=4, nrow=7)
  ggplot2::ggsave(top_plots, filename = file_path, width = 10, height = 16)
}