cap_ips <- function(ips, ip_anno_dt){
  # for motivation: see IP_CSF_over100_breakdown.Rmd
  # does not require ip_anno_dt to be subsetted to the IPs in ips
  
  stopifnot(is.data.table(ips))
  stopifnot(is.data.table(ip_anno_dt))
  
  # set P3 CSFs for IID 53612 to NA
  message("setting  P3 CSFs for IID 53612 to NA")
  P3_CSFs <- intersect(
    ip_anno_dt[source=="P3", set_name, drop=TRUE],
    names(ips)
  )
  target_row <- which(ips$IID==53612)
  for (j in P3_CSFs){
    set(ips, i=target_row, j=j, value=NA)
  }
  
  # cap all other CSFs to 100% (and convert to proportions in [0,1] instead of
  # percentages in [0,100])
  message("capping all other CSFs to 100%")
  CSF_features <- intersect(
    ip_anno_dt[source!="MFI", set_name, drop=TRUE],
    names(ips)
  )
  for (j in CSF_features){
    set(ips, j=j, value=pmin(ips[[j]]/100,1))
  }
  
  ips
}
