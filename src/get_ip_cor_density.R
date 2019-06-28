get_ip_cor_density <- function(corr_like, ips_trans_fam_adj, size_cutoff=500, ip_anno){
  
  ip_source_size <- ip_anno %>% 
    dplyr::count(composite_lin_source) %>% 
    mutate(size=ifelse(n<size_cutoff, "small", "large")) %>% 
    dplyr::select(-n) %>% 
    inner_join(dplyr::select(ip_anno, set_name, composite_lin_source),
               .,
               by="composite_lin_source") %>% 
    dplyr::select(-composite_lin_source) %>% 
    dplyr::filter(set_name %in% names(ips_trans_fam_adj)) %>% 
    split(., .$size)
  
  
  # for test only
  testdat <- ips_trans_fam_adj[,sample(1:ncol(ips_trans_fam_adj),1000)]
  testcor_dt <- WGCNA::cor(testdat, use="p") %>%
    data.table(keep.rownames = TRUE)
  setnames(testcor_dt, "rn", "IP1")
  small <- intersect(names(testcor_dt), ip_source_size$small$set_name)
  large <- intersect(names(testcor_dt), ip_source_size$large$set_name)
  cor_small <- testcor_dt[IP1 %in% small, ..small] %>% as.matrix() %>% .[upper.tri(.)]
  cor_large <- testcor_dt[IP1 %in% large, ..large] %>% as.matrix() %>% .[upper.tri(.)]
  
  # for real
  setnames(corr_like, names(ips_trans_fam_adj))
  corr_like[, IP1:=names(ips_trans_fam_adj)]
  small <- ip_source_size$small$set_name
  large <- ip_source_size$large$set_name
  cor_small <- corr_like[IP1 %in% small, ..small] %>% as.matrix() %>% .[upper.tri(.)]
  cor_large <- corr_like[IP1 %in% large, ..large] %>% as.matrix() %>% .[upper.tri(.)]
  
  # make df for plot
  cor_data <- data.table(
    corr = c(cor_small, cor_large),
    size = c(rep("small",length(cor_small)), rep("large",length(cor_large)))
  )
  
  # plot
  p_corr_density <- cor_data %>%
    ggplot(aes(x=corr, fill=size))+
    geom_density(alpha=0.5)+
    xlab("correlation")+
    scale_fill_discrete(name="subgroup size")+
    theme_bw()
  
  return(p_corr_density)
}
