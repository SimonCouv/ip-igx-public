library(future.batchtools)

future::plan(batchtools_slurm, template = "src/slurm_batchtools.tmpl")

source("src/imports_slurm.R")
source("src/source_slurm_functions.R")
source("src/drake_plan_slurm.R")

config <- drake_config(plan_slurm)
# vis_drake_graph(config)

make(plan_slurm, 
     parallelism = "future", 
     jobs = 10,
     memory_strategy = "memory")  # https://ropenscilabs.github.io/drake-manual/hpc.html#memory