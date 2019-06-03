library(doMPI)
# library(future)
cl <- startMPIcluster()
registerDoMPI(cl)
library(data.table)
library(WGCNA)
library(purrr)                
library(readr)

# parse args ----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
data_fp <- args[1]
output_fp <- args[2]
corfnc <- args[3]

stopifnot(corfnc %in% c("bicor", "pearson", "spearman"))

# setup future  ----------------------------------------------------------------
# makeClusterFunctionsSlurm(template = "slurm", array.jobs = TRUE,
#                           nodename = "localhost", scheduler.latency = 1, fs.latency = 65)

# cat("available cores: ", availableCores(methods = "Slurm"), "\n")
# cat("available workers: ", availableWorkers(methods = "Slurm"), "\n")

# load data  ----------------------------------------------------------------
message(data_fp)
message(output_fp)
data <- fread(data_fp, header = TRUE)
message("data loaded.")
# data <- fread("/home/simon/OneDrive/KCL/Falchi/ip-igx/data/ips_trans_fam_adj.tsv", header = TRUE)

# functions ----------------------------------------------------------------
# 
assign_cols <- function(n_proc, n_feat){
  n_opt <- n_feat*(n_feat-1)/(2*n_proc)
  col <- 1
  block <- 1
  block_col_list <- list()
  while (col <= n_feat){
    n_pairs <- 0
    block_cols <- c()
    while (n_pairs < n_opt & col <= n_feat){
      block_cols <- c(block_cols, col)
      n_pairs <- n_pairs + col   # number of rows in col == col number
      col <- col + 1
    }
    block_col_list[[block]] <- block_cols
    block <- block + 1
  }
  block_col_list
}

block_correlatation <- function(mat, block_cols, ...){
  # correlations are calculated between columns of mat in block_cols, filling up 
  # the corresponding columns of the upper triangular correlation matrix
  
  stopifnot(is.matrix(mat))
  
  # correlation between block cols
  block_cor <- WGCNA::cor(mat[, block_cols], ...)
  # .[upper.tri(., diag = TRUE)]
  
  # correlation between block cols and cols with lower index
  if (min(block_cols) > 1){
    prev_col_cor <- WGCNA::cor(x = mat[, 1:(min(block_cols)-1)], 
                               y = mat[, block_cols],
                               ...)
  } else
    prev_col_cor <- NA
  
  return(list(block_cor=block_cor, prev_col_cor = prev_col_cor))
}

block_bicorrelatation <- function(mat, block_cols, ...){
  # correlations are calculated between columns of mat in block_cols, filling up 
  # the corresponding columns of the upper triangular correlation matrix
  
  stopifnot(is.matrix(mat))
  
  # correlation between block cols
  block_cor <- WGCNA::bicor(mat[, block_cols], ...)
  # .[upper.tri(., diag = TRUE)]
  
  # correlation between block cols and cols with lower index
  if (min(block_cols) > 1){
    prev_col_cor <- WGCNA::bicor(x = mat[, 1:(min(block_cols)-1)], 
                               y = mat[, block_cols],
                               ...)
  } else
    prev_col_cor <- NA
  
  return(list(block_cor=block_cor, prev_col_cor = prev_col_cor))
}

reconstruct_cor_mat <- function(block_cor_list){
  
  # calculate based on inputs
  n_block_cols <- map_dbl(block_cor_list, ~nrow(.x$block_cor))
  n_feat <- sum(n_block_cols)
  
  # initialise
  cor_mat <- matrix(0, n_feat, n_feat)
  
  # fill in first block, no off-diagonal blocks
  cor_mat[1:n_block_cols[1], 1:n_block_cols[1]] <- block_cor_list[[1]]$block_cor
  
  # fill in later blocks
  for (i in 2:length(block_cor_list)){
    block <- block_cor_list[[i]]
    nc <- n_block_cols[[i]]
    last_filled_col <- sum(n_block_cols[1:i-1])
    fill_start <- last_filled_col+1
    fill_end <- last_filled_col + nc
    
    cor_mat[fill_start:fill_end, fill_start:fill_end] <- block$block_cor
    cor_mat[1:(last_filled_col), fill_start:fill_end] <- block$prev_col_cor
    cor_mat[fill_start:fill_end, 1:(last_filled_col)] <- t(block$prev_col_cor)
  }
  cor_mat
}

# execute  ----------------------------------------------------------------

# for testing:
# data <- data[,1:100]

n_proc = mpi.comm.size(0)-1
message("splitting calculations over ", n_proc, " processes")
col_blocks <- assign_cols(n_proc = n_proc,
                          n_feat = ncol(data))
# number of pairs per block
# map(col_blocks, ~sum(.x))

if (corfnc %in% c("pearson", "spearman")){
  cor_l <- foreach(block_cols=col_blocks) %dopar%{
    block_correlatation(
      mat= as.matrix(data),
      block_cols = block_cols,
      use='pairwise.complete.obs',
      method=corfnc
    )
  }
}
if (corfnc == "bicor"){
  cor_l <- foreach(block_cols=col_blocks) %dopar%{
    block_bicorrelatation(
      mat= as.matrix(data),
      block_cols = block_cols,
      use='pairwise.complete.obs',
    )
  }
}


cormat <- reconstruct_cor_mat(cor_l)
write_tsv(as.data.frame(cormat), path=output_fp)

closeCluster(cl)


