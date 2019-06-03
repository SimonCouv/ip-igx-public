plan_slurm = drake_plan(
  
  # see template file specified in make_slurm.R 
  # for valid list names for the resources argument
  
  # correlation matrices -------------------------------------------
  
  ips_pearson = target(fread("data/ips_slurm_blockwise_cor.tsv"),resources = list(ncpus=1, memory=NULL)),
  ips_bicor = target(fread("data/ips_slurm_blockwise_bicor.tsv"),resources = list(ncpus=1, memory=NULL)),
  ips_spear = target(fread("data/ips_slurm_blockwise_spear.tsv"),resources = list(ncpus=1, memory=NULL)),
  
 
  
)

plan_slurm$resources <- lapply(plan_slurm$resources, eval)