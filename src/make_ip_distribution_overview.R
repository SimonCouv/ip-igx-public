make_ip_distribution_overview <- function(ips){
  #ips: features only, no metadata/ID cols
  
  variances <- apply(ips,2,var, na.rm=TRUE)
  means <- apply(ips,2,mean, na.rm=TRUE)
  CVs <- sqrt(variances)/abs(means)
  
  feature_sumstat_df <- cbind(means, variances, CVs, feature = names(ips)) %>% 
    as_tibble() %>%
    tidyr::gather(statistic, value, -feature) %>% 
    left_join(ip_anno[c("set_name", "composite_lin_source")], by=c("feature"="set_name")) %>% 
    mutate(type=ifelse(composite_lin_source=="MFI", "SPEL", "CSF"),
           value=as.numeric(value))
  
  feature_sumstat_df %>% 
    mutate(value=as.numeric(value)) %>% 
    nest(-one_of("type","statistic")) %>% 
    drop_na() %>% 
    mutate(quantiles = map(data, ~quantile(.x$value, na.rm = TRUE)),
           quantiles2=map(quantiles, ~bind_rows(.x) %>% tidyr::gather())) %>% 
    unnest(quantiles2)  %>%
    kable() %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
  
  sumstat_p <-feature_sumstat_df %>% 
    drop_na() %>% 
    ggplot(aes(x="", y=value))+
    geom_boxplot()+
    facet_wrap(~type+statistic, scales="free")
  
  list(df=feature_sumstat_df,p=sumstat_p)
}