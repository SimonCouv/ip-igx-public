source("src/imports_slurm.R")
source("src/source_local_functions.R")

# https://github.com/HenrikBengtsson/future/issues/244
cl_login <- makeClusterPSOCK(
  workers = "login1.rosalind.kcl.ac.uk", 
  outfile = NULL, 
  verbose = TRUE, 
  homogeneous=TRUE,
  # rscript = c("module", "load", "general/R/3.5.0;", "Rscript"), 
  user = "k1893262",
  # rshcmd = "/usr/bin/ssh",
  # rshopts = c("-i /home/simon/.ssh/rosalind/rosalind"),
  revtunnel=TRUE,
  manual=FALSE
)
future::plan(
  list(
    future::tweak(future::cluster, workers=cl_login),
    future::tweak(future.batchtools::batchtools_sge, template = "/users/k1893262/brc_scratch/twinsuk/ip-igx/src/sge_gradual.tmpl")
  )
)



ip_network_pars <- read_tsv(file_in("WGCNA_parameters/ip_pars.tsv")) %>% 
  # ip_network_pars <- read_tsv(file_in("WGCNA_parameters/ip_pars_parsimonious.tsv")) %>% 
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
# config_slurm <- drake_config(plan_slurm)
# vis_drake_graph(config_slurm)

drake::make(plan_slurm,
            parallelism = ssh_backend_future,
            cache=drake_cache("/mnt/lustre/users/k1893262/twinsuk/ip-igx/.drake"),
            jobs = 20,
            jobs_preprocess = 2,
            memory_strategy = "speed",
            verbose = 1,
            retries=3,
            caching = "worker"  # ESSENTIAL, otherwise import (reading) of objects is sequential (because performed by single master process), and becomes the bottleneck with large objects
)

# make(plan_slurm, 
#      parallelism = "future", 
#      verbose = 2
# )


# notes on logging/capturing stdout and stderr --------------------
# 
# https://github.com/HenrikBengtsson/future/issues/270  
#   future now returns stdout per future (instead of all together at completion)
# 
# https://www.jottr.org/2018/07/23/output-from-the-future/
#   future captures stdout since v1.9.0, for uniform behaviour across backends
#   BUT: "Disclaimer: A future’s output is relayed only after it is resolved and
#    when its value is retrieved by the master R process. In other words, the 
#    output is not streamed back in a “live” fashion as it is produced. 
#    Also, it is only the standard output that is relayed. See below, for why 
#    the standard error cannot be relayed."
#    
# https://cran.r-project.org/web/packages/future/vignettes/future-2-output.html
#    stdout (cat, print, str, ...) and messages are all relayed and printed
#    only upon evaluation of the future
#    
# the stdout argument of future::Future (and BatchtoolsFuture and 
# batchtools_slurm and future::plan)
# controls whether standard output is captured
# see ?future::Future
# conditions = character(0) can similarly disable capturing of conditions
# (i.e. messages and warnings)
# https://github.com/HenrikBengtsson/future.batchtools/issues/35
# 
# Excellent explainer of workings of future.batchtools in combination with drake
# https://github.com/HenrikBengtsson/future.batchtools/issues/36
#   - object passing via file system
#   - futures are evaluated by main process -> stdout messages go to
#     respective logs in main process
#   - time and memory profiling in future plans is a to-do
#   - future registries are saved in .future/
#   
# "Standard output (e.g. the default for cat() output) has been captured since 
# future 1.9.0 (2018-07-23) and message & warning conditions have been captured
# since future 1.11.0 (2019-01-22). Captured output and conditions are then relayed
#  in the main/parent R process."
# https://stackoverflow.com/questions/54443239/log-r-messages-and-errors-when-using-future-batchtools
#  
# from the batchtools side, it is possible to peek into the registry 
# (and potentially the logs in the registry, is stdout=NA in future?)
# https://github.com/mllg/batchtools/issues/157
# e.g. loadRegistry(".future/20190610_135505-xim77l/ip_module_stats_signed_hybrid_ips_bicor_4_1182180771", writeable=FALSE)
# 




