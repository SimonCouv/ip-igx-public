parse_ip_subsets_helper <- function(x,ips_trans_fam_adj){
  map(x$feature_data, ~parse_ip_subsets(feat_data=.x,ips_trans_fam_adj=ips_trans_fam_adj))
}