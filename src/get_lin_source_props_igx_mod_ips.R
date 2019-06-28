get_lin_source_props_igx_mod_ips <- function(target_indices, data, ip_anno){
  
  setDT(data)
  
  data_ipvars <- data[, .(IP, composite_lin_source)]
  
  ip_anno_count <- count(ip_anno, composite_lin_source)
  
  # proportions for each target_index
  running_lin_source_prop <- list()
  cnt <- 0
  for (i in target_indices){
    if (cnt %% 100 == 0) print(glue("iter {cnt} - index {i}"))
    running_lin_source_prop[[as.character(i)]] <-data_ipvars %>% 
      head(i) %>% 
      count(composite_lin_source) %>%
      left_join(ip_anno_count,
                by = "composite_lin_source", 
                suffix = c("_assoc", "_feat")) %>% 
      mutate(n_weighted = n_assoc/n_feat,
             prop_weighted = n_weighted/sum(n_weighted),
             prop = n_assoc/sum(n_assoc))
    cnt <- cnt+1
  }
  
  # get the full list of lin_sources in the analysis, by extracting at the last target index
  lin_sources <- running_lin_source_prop[[length(running_lin_source_prop)]] %>%
    dplyr::select(composite_lin_source) %>% 
    unlist(.)
  
  # collect proportions per lin_source, over target indices
  lin_source_props <- list()
  for (lin_source in lin_sources){
    print(lin_source)
    lin_source_props[[lin_source]] <- list(
      "unweighted" = map(running_lin_source_prop,
                         ~.x$prop[.x$composite_lin_source==lin_source]),
      "weighted" = map(running_lin_source_prop,
                       ~.x$prop_weighted[.x$composite_lin_source==lin_source])
    )
  }
  lin_source_props_df <- lin_source_props %>% 
    purrr::transpose(.) %>% 
    map(~map(.x, ~as.numeric(.x)) %>%
          bind_cols() %>%
          mutate(top = target_indices)) %>%
    bind_rows(.id="weighting")
  
  lin_source_props_df
}