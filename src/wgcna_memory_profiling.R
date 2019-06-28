pm <- profmem(
  glycan_module_stats <- pmap(
    .l = glycan_network_pars,
    .f = wgcna_parameter_search,
    data=glycans_trans_fam_adj,
    cut_pars=expand.grid(deep =c(4, 1),
                         minModuleSize=c(3,5,10,20)),
    create_plot = TRUE,
    feat_modality=glycan_modality
  ) %>% 
    map_dfc(tibble) %>% 
    t() %>%
    as_tibble() %>% 
    set_names(names(.[[1]])) %>% 
    bind_cols(., glycan_network_pars) %>% 
    unnest()
)

profvis(
  glycan_module_stats <- pmap(
    .l = glycan_network_pars,
    .f = wgcna_parameter_search,
    data=glycans_trans_fam_adj,
    cut_pars=expand.grid(deep =c(4, 1),
                         minModuleSize=c(3,5,10,20)),
    create_plot = TRUE,
    feat_modality=glycan_modality
  ) %>% 
    map_dfc(tibble) %>% 
    t() %>%
    as_tibble() %>% 
    set_names(names(.[[1]])) %>% 
    bind_cols(., glycan_network_pars) %>% 
    unnest()
)