fam_adj_collect <- function(fam_adj_l, raw){
  fam_adj_l %>%
    future_map(~tibble(res=.x$res, index=names(.x$res)) %>%
                 dplyr::right_join(tibble(index=as.character(1:nrow(raw))),by="index") %>%  #make residuals complete (i.e. fill in missing lines)
                 dplyr::select(res)) %>%
    dplyr::bind_cols() %>%
    rlang::set_names(names(fam_adj_l))
}