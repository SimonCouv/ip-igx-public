make_module_heatmap_pdf <- function(plots, file_path, pointsize, width, height, network_pars){
  
  dir_create(path_dir(file_path))  # mkdir -p
  
  pdf(file=file_path,
      pointsize = pointsize,
      width = width,
      height = height)
  
  for (i in seq_along(plots)){
    network_string <- glue("{network_pars$corfnc[i]}{network_pars$softpower[i]}/{network_pars$networktype[i]}/{network_pars$method[i]}")
    draw(plots[[i]], column_title = glue("Network {i} - {network_string}"))
  }
  
  dev.off()
}