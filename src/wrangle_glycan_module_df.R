wrangle_glycan_module_df <- function(modules, network_pars, dat){

  tmp <- modules %>% 
    bind_rows() %>% 
    bind_cols(unnest(network_pars),.) %>% 
    mutate(
      merge_pars_used = map(merge_pars_used,
                            ~map(.x,~t(.x)) %>% 
                              do.call(rbind, .) %>% 
                              as_tibble() %>% 
                              mutate(cut=as.character(1:nrow(.)))),
      feature_stats = map(color_overview, 
                          ~mutate(.x, feat=colnames(dat)) %>% 
                            tidyr::gather(cut, colors, -feat) %>% 
                            nest(-cut, .key="feature_stats")),
      mod_modularities = map(mod_modularities, 
                             ~bind_rows(.x, .id="cut") %>% 
                               nest(-cut, .key="module_stats")),
      cut_stats = pmap(
        .l=list(merge_pars_used, mod_modularities, feature_stats),
        .f=function(df1,df2,df3) df1 %>% full_join(df2, by="cut") %>% full_join(df3, by="cut")
      )
    ) %>% 
    dplyr::select(corfnc, networktype, method, powers, plots, cut_stats)
    
  
  
  
  cat("adding module pc1 scores to module_stats\n")
  tmp %>% 
    add_module_pc1(modules_df=.,
                   dat=dat,
                   impute=FALSE,
                   nPC=1,
                   excludeGrey = FALSE,
                   subHubs = FALSE,
                   softPower = 1,
                   scale=FALSE,
                   trapErrors=FALSE,
                   verbose=5)
}