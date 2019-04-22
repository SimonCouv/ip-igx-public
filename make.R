# based on https://github.com/wlandau/drake-examples/blob/master/gsp/make.R, 21/04/2019

# explicit library imports and functions definitions defined
# here in the master script for sake of simplicity


library(qqman)
library(drake)
library(digest)

make_qq <- function(v_pval, fig_path){
  png(fig_path, res=600,
      width = 6, height = 6, units = "in")
  qqman::qq(v_pval)
  dev.off()
  
}

source(here::here("my_drake_plan.R"))

# graph of workflow
config <- drake_config(plan) # nolint
vis_drake_graph(config)

make(plan)
