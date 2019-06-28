add_module_pc1 <- function(modules_df, dat, ...){
  
  # prepare 
  tmp <- modules_df %>% 
    dplyr::select(corfnc, networktype, powers, cut_stats) %>% 
    unnest() %>% 
    ungroup()
  
  # calculate MEs
  cat("calculating module pc1s in parallel")
  MEs = future_map(
    .x = tmp$feature_stats,
    .f = function(x, ...){
      res <- moduleEigengenes(colors=x$colors,
                              expr = dat,
                              ...) %>%  # pass the ellipsis arguments through to moduleEigengenes()
        magrittr::extract(c("varExplained", "eigengenes"))
      res$varExplained <- t(res$varExplained) %>% 
        data.frame(varExplained = .,
                   colors = sort(unique(x$colors)),
                   stringsAsFactors = FALSE)
      res$eigengenes <- as_tibble(res$eigengenes) %>% 
        tidyr::gather(colors, score) %>% 
        nest(score, .key = "score") %>% 
        mutate(colors = str_match(colors, pattern = "^ME(.*)")[,2])
      return(res)
    },
    ...  # pass the ellipsis arguments through to .f
  )

  # add MEs into df
  res <- tmp %>%
    mutate(
      MEs = MEs,
      module_size = map(feature_stats, ~count(.x,colors)),
      module_stats = pmap(
        .l = list(module_stats, MEs, module_size),
        .f = function(x, y, z){
          reduce(list(x, y$varExplained, y$eigengenes, z),
                 .f=full_join,
                 by="colors")
        }
      )
    ) %>% 
    # format for output
    group_by(corfnc, networktype, powers) %>% 
    dplyr::select(-one_of("MEs", "module_size")) %>% 
    nest(.key = "cut_stats") %>% 
    inner_join(dplyr::select(modules_df, -cut_stats),
               by=c("corfnc", "networktype", "powers")) %>% 
    mutate(
      cut_stats = map(cut_stats, 
                      ~dplyr::select(.x, 
                                     c("cut",
                                       "module_stats", "feature_stats",
                                       "maxCoreScatter","minGap",
                                       "maxAbsCoreScatter","minAbsGap")
                      )
      )
    )
    
  # keep plots if applicable  
  if ("plots" %in% names(modules_df)){
    res <- dplyr::select(res, corfnc, networktype, powers, method, plots, cut_stats)
  } else {
    res <- dplyr::select(res, corfnc, networktype, powers, method, cut_stats)
  }
  
  return(res)
}



