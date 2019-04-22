# based on https://github.com/wlandau/drake-examples/blob/master/gsp/make.R, 21/04/2019

source(here::here("src/imports.R"))
source(here::here("src/source_local_functions.R"))
source(here::here("src/dir_defs.R"))
source(here::here("src/my_drake_plan.R"))
# source(here::here())


# graph of workflow
config <- drake_config(plan) # nolint
vis_drake_graph(config)

make(plan)
# make(plan, jobs = 2)
