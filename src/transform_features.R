transform_features <- function(raw){
  # raw: class data.table with features in columns, samples in rows
  # requires all features are in [0,1]
  
  # assertions
  stopifnot(is.data.table(raw))
  stopifnot(all(unlist(raw[,lapply(.SD, is.numeric)])))
  stopifnot(all(unlist(raw)>=0, na.rm = TRUE) & all(unlist(raw)<=1, na.rm = TRUE))
  
  # perform transformations
  for (j in names(raw)){
    set(raw, j=j, value=asin(sqrt(raw[[j]])))          # arcsin(sqrt)
  }
  raw %>% 
    scale(.) %>%                                       # scale
    {t(limma::normalizeQuantiles(t(.)))}               # quantile normalise
}