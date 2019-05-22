filter_ips <- function(ips, ip_anno_dt, row_miss_threshold, col_miss_threshold){
  setDT(ip_anno_dt)
  
  #setup
  ips_excluding_old_p6 <- ip_anno_dt[!(source %in% c("P6", "P6 earlyB")),set_name]
  ips_final_max_qc <- ip_anno_dt[final_qc_max=="OK",set_name]
  
  #filter
  cbind(
    ips[,1:2],
    ips[,-(1:2)] %>% 
      .[, colMeans(is.na(.)) < col_miss_threshold, with=FALSE] %>%
      .[, lapply(., var, na.rm=TRUE)>0, with=FALSE] %>% 
      .[,intersect(names(.),ips_excluding_old_p6), with=FALSE] %>% 
      .[,intersect(names(.),ips_final_max_qc), with=FALSE]
  ) %>% 
    .[rowMeans(is.na(.[,-(1:2)])) < row_miss_threshold, ]
}