bind_rows_plan <- function(..., drakeplan, vars){
  # function to extend drake API
  # partially based on function 'post_process' in https://ropenscilabs.github.io/drake-manual/plans.html#the-types-of-transformations
  
  args <- list(...)
  names(args) <- all.vars(substitute(list(...)))
  
  # Get the metadata of the combined targets in the list.
  metadata <- dplyr::filter(drakeplan, target %in% names(args)) %>% .[,c("target", vars)]
  
  bind_rows(args, .id = "target") %>% 
    left_join(metadata)
}