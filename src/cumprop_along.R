cumprop_along <- function(data, v, along, do.order=TRUE, ...){
  # ... passed down to order(), e.g. set decreasing=TRUE
  
  setDT(data)
  
  if (do.order){
      data <- data[order(data[[along]], ...),]
  }
  set(data, i=NULL, j=v, as.factor(data$along))
  setkey(data, along)
  
  cumprop <- list(length(unique(along)))
  for (step in unique(along)){
    i <- max(which(data$along == step))
    cumprop[[as.character(along)]] <- head(data, i) %>% table()
    
  }
}




debug(cumprop_along)
cumprop_along(iris, "Species", "Petal.Length")
