compare_input_imputed <- function(raw_df_l, imputed_df_l){
  #inputs: only data, no metadata columns, need same dimensions
  # assumes same ordering of rows in both DFs!
  
  if(!all.equal(map_dfc(raw_df_l, dim),map_dfc(imputed_df_l, dim))) "dimensions are not pairwise equal"
  
  f <- function(raw_df, imputed_df){
    
    mask <- ifelse(is.na(raw_df), "imputed","raw") %>% 
      as_tibble() %>% 
      mutate(pseudo_ID=1:nrow(.)) %>% 
      tidyr::gather(feature, origin, -pseudo_ID)
    
    comp <- imputed_df %>% 
      as_tibble() %>% 
      mutate(pseudo_ID=1:nrow(.)) %>% 
      tidyr::gather(feature, value, -pseudo_ID) %>% 
      left_join(mask, by=c("pseudo_ID", "feature"))
    
    comp
  }
  
  p <- map2(.x=raw_df_l,
       .y=imputed_df_l,
       .f = f) %>% 
    bind_rows(.id="dataset") %>% 
    ggplot(aes(x=feature, y=value, fill=origin, color=origin))+
    geom_boxplot(outlier.size=1, varwidth = TRUE)+
    coord_flip()+
    theme_bw()+
    theme(axis.text=element_text(size=6))+
    facet_wrap(~dataset, nrow=1, scales = "free")
  
  p
}