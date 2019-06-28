library(future.batchtools)

future::plan(batchtools_slurm, template = "src/slurm_batchtools.tmpl")

source("src/imports_slurm.R")
source("src/source_local_functions.R")
ip_network_pars <- read_tsv(file_in("WGCNA_parameters/ip_pars_parsimonious.tsv")) %>% 
  dplyr::rename("nwtype"="networktype") %>%
  mutate(
    corr_like = paste0("ips_",dplyr::recode(corfnc, cor="pearson", spearman = "spear")),
    similarity= paste("similarity", nwtype, corr_like, sep="_"),
  ) %>%
  # group_by(corr_like, nwtype, method, similarity) %>%    # don't use looping over powers, to save total memory requirement per SLURM job
  # summarise(powers=list(softpower)) %>% 
  # ungroup() %>% 
  mutate_at(.vars=vars(corr_like, similarity), rlang::syms)

source("src/drake_plan_slurm.R")
config_slurm <- drake_config(plan_slurm)
# vis_drake_graph(config_slurm)
outdated(config_slurm)