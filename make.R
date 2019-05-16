# based on https://github.com/wlandau/drake-examples/blob/master/gsp/make.R, 21/04/2019

# set options
options(clustermq.scheduler = "multicore")

source(here::here("src/imports.R"))
source(here::here("src/source_local_functions.R"))
source(here::here("src/dir_defs.R"))
source(here::here("src/my_drake_plan.R"))
# source(here::here())

# see https://ropenscilabs.github.io/drake-manual/hpc.html#parallel-computing-within-targets
# furrr recommends (https://davisvaughan.github.io/furrr/index.html):
# future::plan(multiprocess)
# BUT, this does not work with drake. ERROR: " Having problems with parallel::mclapply(), future::future(), or furrr::future_map() in drake? Try one of the workarounds at https://ropenscilabs.github.io/drake-manual/hpc.html#parallel-computing-within-targets or https://github.com/ropensci/drake/issues/675. "
# instead, use future.callr packages, which provides the futures API (to work with furrr) on top of callr backend (to play nice with drake).

future::plan(future.callr::callr, workers = 3L)

# graph of workflow
config <- drake_config(plan)
vis_drake_graph(config)

# make(plan)
make(plan,  
     parallelism = "clustermq", 
     jobs =4,
     cache_log_file = TRUE)
