make_modularity_overview <- function(mod_modularities_l, glycans_ME_noderiv){
  
  mod_modularities_l %>% 
    map(~dplyr::bind_rows(.x, .id="cut")) %>%
    dplyr::bind_rows(.id="nw") %>% 
      dplyr::mutate(cut = str_match(cut, pattern = "cut_iter_([:digit:]+)")[,2]) %>%
    group_by(nw, cut) %>% 
    nest(.key="modularities") %>% 
    inner_join(glycans_ME_noderiv)
  
  # # make modularity overview
  # modularity_overview_nogrey <- mod_modularities_l %>%
  #   map(~map(.x, ~rownames_to_column(.x) %>%
  #              dplyr::filter(rowname != "grey") %>%
  #              dplyr::select(-rowname) %>%
  #              dplyr::summarise_all(sum) %>%
  #              set_names(paste0(names(.), "_no_grey"))) %>% 
  #         dplyr::bind_rows(.id="cut")) %>%
  #   dplyr::bind_rows(.id="nw")
  # modularity_overview_inclgrey <- mod_modularities_l %>%
  #   map(~map(.x, ~dplyr::summarise_all(.x, sum) %>%
  #              set_names(paste0(names(.), "_incl_grey"))) %>% 
  #         dplyr::bind_rows(.id="cut")) %>%
  #   dplyr::bind_rows(.id="nw")
  # 
  # modularity_overview <- inner_join(modularity_overview_nogrey, 
  #                                   modularity_overview_inclgrey,
  #                                   by=c("cut", "nw")) %>% 
  #   dplyr::mutate(cut = str_match(cut, 
  #                                 pattern = "cut_([:digit:]+)")[,2]) %>%
  #   dplyr::arrange(desc(adj_mod_modularities_no_grey)) %>%
  #   dplyr::left_join(rownames_to_column(network_pars),
  #                    by=c("nw" = "rowname")) %>% 
  #   dplyr::left_join(rownames_to_column(cut_pars), 
  #                    by=c("cut" = "rowname"))
}