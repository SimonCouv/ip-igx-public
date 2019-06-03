make_sft_from_similarity <- function(similarity_l, common_pars){

  capture.output(
    res <- map_depth(.x = similarity_l,
                     .depth = 2,
                     .f = function(x){
                       do.call(
                         what=pickSoftThreshold.fromSimilarity,
                         args = c(list(similarity=similarity), common_pars)
                       )
                     }) %>% 
      bind_rows(., .id="networktype") %>% 
      bind_rows(., .id="corFnc")
  )
  res
}