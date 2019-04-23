# based on https://ropenscilabs.github.io/drake-manual/projects.html#safer-interactivity, 21/04/2019

# source(here::here())


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

drake_config(plan, verbose = 2)
