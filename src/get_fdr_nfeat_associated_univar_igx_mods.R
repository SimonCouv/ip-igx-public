get_fdr_nfeat_associated_univar_igx_mods <- function(fdr_thresholds, data, p_adj="pv_adj_global"){
  n_associated <- list()
  for (fdr_threshold in fdr_thresholds){
    n_associated[[as.character(fdr_threshold)]] <- data %>%
      .[.[[p_adj]] < fdr_threshold, ] %>% 
      summarise(n_modules = n_distinct(module_tag),
                n_IP = n_distinct(subset_name),
                # prop_igx = n_igx/unique(total_igx),
                # prop_ip = n_ip/unique(total_ip),
                n_IgG = unique(module_tag) %>% str_detect(., "G") %>% sum(.),
                n_IgA = n_modules - n_IgG)
  }
  n_associated %>% 
    dplyr::bind_rows(n_associated, .id = "fdr_threshold")%>%
    tidyr::gather(key, n, -fdr_threshold) %>%
    tidyr::separate(key, into = c("type", "omic")) %>%
    dplyr::group_by(omic) %>%
    dplyr::mutate(prop=n/max(n),
                  fdr_threshold = as.numeric(fdr_threshold))
}