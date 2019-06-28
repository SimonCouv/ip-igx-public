# NEEDS REVISION BEFORE USE!
# might become redundant if drake plans gets fixed


source("src/imports_slurm.R")
source("src/source_local_functions.R")

loadd(ip_anno)
loadd(ips_trans_fam_adj)

print("pearson")
Sys.time()
loadd(ips_pearson)
p_ips_cor_density_pearson = get_ip_cor_density(corr_like = ips_pearson,
                                               ips_trans_fam_adj = ips_trans_fam_adj,
                                               ip_anno = ip_anno)
saveRDS(p_ips_cor_density_pearson, file = "data/p_ips_cor_density_pearson.RDS")
rm(ips_pearson, p_ips_cor_density_pearson)
gc()

print("pearson")
Sys.time()
loadd(ips_bicor)
p_ips_cor_density_bicor = get_ip_cor_density(corr_like = ips_bicor,
                                             ips_trans_fam_adj = ips_trans_fam_adj,
                                             ip_anno = ip_anno)
saveRDS(p_ips_cor_density_bicor, file = "data/p_ips_cor_density_bicor.RDS")
rm(ips_bicor, p_ips_cor_density_bicor)
gc()

print("pearson")
Sys.time()
loadd(ips_spear)
p_ips_cor_density_spear = get_ip_cor_density(corr_like = ips_spear,
                                             ips_trans_fam_adj = ips_trans_fam_adj,
                                             ip_anno = ip_anno)
saveRDS(p_ips_cor_density_spear, file = "data/p_ips_cor_density_spear.RDS")
rm(ips_spear, p_ips_cor_density_spear)
gc()



