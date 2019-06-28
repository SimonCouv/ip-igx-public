get_fdr_nfeat_associated_univar <- function(fdr_thresholds, ip_igx_univar_good, p_adj="pv_adj_global"){
  n_associated <- list()
  for (fdr_threshold in fdr_thresholds){
    n_associated[[as.character(fdr_threshold)]] <- ip_igx_univar_good %>%
      # mutate(total_IgX = n_distinct(response),
      #        total_IgG = unique(response) %>% str_detect(., "IgG") %>% sum(.),
      #        total_IgA = total_IgX - total_IgG,
      #        total_IP = n_distinct(predictor)) %>%
      dplyr::filter(!!p_adj < fdr_threshold) %>%
      summarise(n_IgX = n_distinct(response),
                n_IP = n_distinct(predictor),
                # prop_igx = n_igx/unique(total_igx),
                # prop_ip = n_ip/unique(total_ip),
                n_IgG = unique(response) %>% str_detect(., "IgG") %>% sum(.),
                n_IgA = n_IgX - n_IgG)
  }
  n_associated %>% 
    dplyr::bind_rows(n_associated, .id = "fdr_threshold")%>%
    tidyr::gather(key, n, -fdr_threshold) %>%
    tidyr::separate(key, into = c("type", "omic")) %>%
    dplyr::group_by(omic) %>%
    dplyr::mutate(prop=n/max(n),
           fdr_threshold = as.numeric(fdr_threshold))
}