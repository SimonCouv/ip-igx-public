# based on https://github.com/wlandau/drake-examples/blob/master/gsp/make.R, 21/04/2019

# set options
options(clustermq.scheduler = "multicore")


source(here::here("src/imports.R"))
source(here::here("src/source_local_functions.R"))
source(here::here("src/dir_defs.R"))
source(here::here("src/my_drake_plan.R"))
# source(here::here())

# future::plan(
#   list(
#     # future::multisession,
#     # future::sequential,
#     future::tweak(future.callr::callr),
#     future::tweak(future.callr::callr, workers = 3L)
#     # future::tweak(future.callr::callr)
#   )
# )

# see https://ropenscilabs.github.io/drake-manual/hpc.html#parallel-computing-within-targets
# furrr recommends (https://davisvaughan.github.io/furrr/index.html):
# future::plan(multiprocess)
# BUT, this does not work with drake. ERROR: " Having problems with parallel::mclapply(), future::future(), or furrr::future_map() in drake? Try one of the workarounds at https://ropenscilabs.github.io/drake-manual/hpc.html#parallel-computing-within-targets or https://github.com/ropensci/drake/issues/675. "
# instead, use future.callr packages, which provides the futures API (to work with furrr) on top of callr backend (to play nice with drake).

# future::plan(future.callr::callr, workers = 3L)

# graph of workflow
config <- drake_config(plan)
# vis_drake_graph(config)

# make(plan)
make(plan,  
     parallelism = "clustermq",
     # parallelism = "future",
     jobs =10,
     jobs_preprocess=4,
     verbose = 2,
     # prework = list(quote(future::plan(future.callr::callr, workers=2L)),
     #                quote(options(future.globals.maxSize =  4*1024^3))),
     # cache_log_file = TRUE,
     keep_going = FALSE)

# reason for passing future.globals.maxSize:  # 4GB  https://github.com/HenrikBengtsson/future/issues/185
# to avoid error:
# simpleError in getGlobalsAndPackages(expr, envir = envir, globals = globals): 
# The total size of the 10 globals that need to be exported for the future 
# expression (‘{; ...future.f.env <- environment(...future.f); if
#  (!is.null(...future.f.env$`~`)) {; if (is_bad_rlang_tilde(...future.f.env$`~`)) 
#  {; ...future.f.env$`~` <- base::`~`; }; ...; .out; }); }’) is 665.12 MiB.
#   This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize'). 
#   The three largest globals are ‘...future.x_ii’ (387.74 MiB of class ‘list’), 
#   ‘raw’ (277.24 MiB of class ‘list’) and ‘%>%’ (101.66 KiB of class ‘function’).
