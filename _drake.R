# based on https://ropenscilabs.github.io/drake-manual/projects.html#safer-interactivity, 21/04/2019

source(here::here("src/imports.R"))
source(here::here("src/source_local_functions.R"))
source(here::here("src/my_drake_plan.R"))
source(here::here("src/dir_defs.R"))
# source(here::here())


# graph of workflow
# config <- drake_config(plan) # nolint
# vis_drake_graph(config) 

drake_config(plan, verbose = 2)
