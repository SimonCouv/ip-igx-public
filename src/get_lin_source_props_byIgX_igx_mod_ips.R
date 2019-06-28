get_lin_source_props_byIgX_igx_mod_ips <- function(target_indices, data, igx_var=IgX){
  # proportions for each target_index
  running_lin_source_IgX_prop <-  list()
  igx_var=enquo(igx_var)
  cnt <- 0
  cat("running over indices\n")
  for (i in target_indices){
    if (cnt %% 100 == 0) print(glue("iter {cnt} - index {i}"))
    running_lin_source_IgX_prop[[as.character(i)]] <- data %>%
      head(i) %>% 
      count(!!igx_var, composite_lin_source) %>%
      dplyr::group_by(!!igx_var) %>% 
      dplyr::mutate(prop = n/sum(n))
    cnt <- cnt+1
  }
  
  # browser()
  lin_sources <- running_lin_source_IgX_prop[[length(target_indices)]] %>%
    dplyr::ungroup() %>%
    dplyr::distinct(composite_lin_source) %>% 
    unlist(.)
  
  # collect proportions for each lin_source and IgX combination in the data
  cat("collecting proportions\n")
  lin_source_IgX_props <- list()
  for (lin_source in lin_sources){
    print(lin_source)
    
    for (IgX in c("IgG", "IgA")){
      lin_source_IgX_props[[lin_source]][[IgX]] <- purrr::map(running_lin_source_IgX_prop,
                                                              ~.x$prop[.x$composite_lin_source==lin_source & .x[[igx_var]] == IgX])
    }
  }
  
  # mangle lists to df
  lin_source_IgX_props_df <- purrr::map(lin_source_IgX_props, 
                                        ~purrr::map(.x, ~as.numeric(.x))) %>%
    purrr::transpose(.) %>%
    purrr::map(~bind_cols(.x) %>%
                 dplyr::mutate(top = as.numeric(names(running_lin_source_IgX_prop)))
    ) %>%
    dplyr::bind_rows(.id="IgX")
  
  lin_source_IgX_props_df
}