histogram_n <- function(data, mapping,anno_size=4){
  x_name <- quo_text(mapping$x)
  n <- nrow(data %>% dplyr::select(c(x_name)) %>% drop_na())
  ggplot(data=data,mapping=mapping)+
    geom_histogram(bins = floor(n/25))+
    annotate("text", Inf, Inf, label=sprintf("n=%s",n),
             vjust = "inward", hjust = "inward", size = anno_size)
}