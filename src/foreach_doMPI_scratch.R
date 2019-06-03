library(doMPI)
library(foreach)
cl <- startMPIcluster()
registerDoMPI(cl)


df <- foreach(i=1:500,
              .combine="rbind",
              .options.mpi=list(chunkSize=10)) %dopar% {
  
  # Sys.sleep(1)
  data.frame(
    pid=Sys.getpid(), 
    nodename=Sys.info()['nodename'],
    n_cores_slurm = as.numeric(future::availableCores(methods = "Slurm")),
    n_cores = as.numeric(future::availableCores()),
    mpi_size = mpi.comm.size(0)
  )
}


# df
table(df[,1:2])
summary(df$n_cores_slurm)
summary(df$n_cores)

closeCluster(cl)