library(future.batchtools)
library(furrr)
library(magrittr)
library(dplyr)
future::plan(
  list(
    tweak(batchtools_slurm, 
          template="src/slurm-simple.tmpl",
          resources = list(ncpus = 10, 
                           memory = 1000,
                           ntasks = 3,
                           walltime=1000)),
    multicore
  ))

cat("available cores: ", availableCores(methods = "Slurm"), "\n")
cat("available workers: ", availableWorkers(methods = "Slurm"), "\n")

test %<-% 
  future_map(1:500, ~{Sys.sleep(1);data.frame(
    pid=Sys.getpid(), 
    nodename=Sys.info()['nodename'])})
# test %<-% 
#   future_map(1:5, ~data.frame(
#     pid=Sys.getpid(), 
#     nodename=Sys.info()['nodename']))

bind_rows(test) %>% 
  count(pid,nodename)
# x %<-% Sys.getpid()
# y %<-% Sys.getpid()
# list(x,y)