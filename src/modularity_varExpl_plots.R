modularity_varExpl_plots <- function(modules_df, varExpl_trh, ...){
  # ... passed to cowplot
  # 
  
  create_cowplot <- function(dat, reordervar2,...){
    rv2 <- enquo(reordervar2)
    # browser()
    # list(
    #   sort_var_expl =
        # dat %>% 
        # arrange(median_var_expl)%>% 
        # ungroup() %>% 
        # mutate(id = factor(id, unique(id))) %>% 
        # ggplot(aes(x=id, y=value, color=criterion))+
        # geom_boxplot(varwidth = TRUE)+
        # coord_flip()+
        # facet_wrap(~criterion, ncol=6, scales = "free_x")+
        # scale_color_discrete(guide=FALSE)
      
      # sort_adj_mod = 
        dat %>%
        arrange(!!rv2)%>%
        ungroup() %>%
        mutate(id = factor(id, unique(id))) %>%
        ggplot(aes(x=id, y=value, color=criterion))+
        geom_boxplot(varwidth = TRUE)+
        coord_flip()+
        facet_wrap(~criterion, ncol=6, scales = "free_x")+
        scale_color_discrete(guide=FALSE)
    # ) %>% 
    #   cowplot::plot_grid(plotlist = ., labels = names(.), ...)
  }
  
  plotdat <- modules_df %>% 
    unnest(cut_stats) %>% 
    unnest(module_stats) %>% 
    unite(corfnc, networktype,powers,cut, col="id") %>% 
    group_by(id) %>% 
    dplyr::rename("module_size"="n") %>% 
    mutate(n_modules = n(),
           median_var_expl = median(varExplained),
           median_adj_mod = median(adj_mod_modularities),
           cut_sum_adj_mod = sum(adj_mod_modularities),
           cut_sum_adj_mod2 = sum(adj_mod_modularities),
           cut_sum_cor_mod = sum(cor_mod_modularities),
           cut_sum_TOM_mod = sum(TOM_mod_modularities),
           cut_sum_adj_mod_constr = sum(adj_mod_modularities*(varExplained>varExpl_trh)),
           cut_sum_adj_mod_constr2 = sum(adj_mod_modularities*(varExplained>varExpl_trh)),
           cut_sum_cor_mod_constr = sum(cor_mod_modularities*(varExplained>varExpl_trh)),
           cut_sum_TOM_mod_constr = sum(TOM_mod_modularities*(varExplained>varExpl_trh)),
           `weighted modularity` = sum(adj_mod_modularities*varExplained),
           cut_sum_adj_mod_weighted2 = sum(adj_mod_modularities*varExplained),
           cut_sum_cor_mod_weighted = sum(cor_mod_modularities*varExplained),
           cut_sum_TOM_mod_weighted = sum(TOM_mod_modularities*varExplained)
    )
  
  # p_mod_distr <- create_cowplot(plotdat %>% 
  #                    tidyr::gather(criterion, value, 
  #                                  c("varExplained", "adj_mod_modularities", "TOM_mod_modularities",
  #                                    "cor_mod_modularities", "n_modules", "module_size")),
  #                  reordervar2 = median_adj_mod, ...)
  
  # p_mod_sum <- plotdat %>%
  #   tidyr::gather(criterion, value,
  #                 c("varExplained", "cut_sum_adj_mod", "cut_sum_cor_mod",
  #                   "cut_sum_TOM_mod", "n_modules", "module_size")) %>%
  #   create_cowplot(reordervar2 = cut_sum_adj_mod2, ...)
  # 
  # p_mod_sum_varexpl_constraint <- plotdat %>%
  #   tidyr::gather(criterion, value,
  #                 c("varExplained", "cut_sum_adj_mod_constr", "cut_sum_cor_mod_constr",
  #                   "cut_sum_TOM_mod_constr", "n_modules", "module_size")) %>%
  #   create_cowplot(reordervar2 = cut_sum_adj_mod_constr2, ...)

  p_mod_sum_varexpl_weighted <- plotdat %>%
    tidyr::gather(criterion, value,
                  c("varExplained",
                    `weighted modularity`, "n_modules", "module_size")) %>%
    mutate(criterion=recode(criterion, 
                            module_size = "module size",
                            n_modules = "#modules",
                            varExplained = "variance explained"
                            )) %>% 
    create_cowplot(reordervar2 = cut_sum_adj_mod_weighted2, ...)

  # list(p_mod_distr = p_mod_distr,
  #      p_mod_sum = p_mod_sum,
  #      p_mod_sum_varexpl_constraint = p_mod_sum_varexpl_constraint,
  #      p_mod_sum_varexpl_weighted = p_mod_sum_varexpl_weighted)
  
  p_mod_sum_varexpl_weighted
  
}

