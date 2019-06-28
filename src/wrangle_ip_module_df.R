wrangle_ip_module_df <- function(ip_module_stats_df, dat){
  
  tmp <- ip_module_stats_df %>% 
    mutate_at(vars(corfnc, nwtype, method), .f=~str_remove_all(., '\\"')) %>% 
    dplyr::rename("networktype"="nwtype") %>% 
    mutate(powers=as.numeric(softpower)) %>% 
    dplyr::select(-one_of(c("target", "softpower"))) %>% 
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
    dplyr::select(corfnc, networktype, method, powers, cut_stats)
  
  
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

