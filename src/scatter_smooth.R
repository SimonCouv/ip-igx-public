scatter_smooth <- function(data, mapping,size=4, corfnc = "spearman"){
  x_name <- quo_text(mapping$x)
  y_name <- quo_text(mapping$y)
  x_pos = 0.5*(max(data[,x_name], na.rm = TRUE)-min(data[,x_name], na.rm = TRUE))
  y_pos =  0.9*(max(data[,y_name], na.rm = TRUE)-min(data[,y_name], na.rm = TRUE))
  sp_cor = round(cor(data[,x_name], data[,y_name], 
                     method = corfnc, 
                     use = "pairwise.complete.obs"),
                 digits=3
  )
  cor_prefix <- substr(corfnc,1,2)
  
  n <- nrow(data %>% dplyr::select(c(x_name, y_name)) %>% drop_na())
  ggplot(data=data,mapping=mapping)+
    geom_point(size=0.5, alpha=0.3)+
    geom_smooth(alpha=0.3, color="black", size=0.5)+
    # geom_smooth(alpha=0.3, color="black", size=0.5, method = "gam", formula = y ~ s(x, bs = "cs"))+
    annotate("text", Inf, Inf, label=sprintf("n=%s\n%s cor=%s",n, cor_prefix,sp_cor),vjust = "inward", hjust = "inward",
             size=size)
  # stat_n_text(aes(x=10,y=10, group="combo"))
}