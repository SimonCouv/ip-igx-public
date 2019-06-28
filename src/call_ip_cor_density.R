#!/bin/bash

# NEEDS REVISION BEFORE USE!
# might become redundant if drake plans gets fixed

#SBATCH -n 1
##SBATCH --cpus-per-task 1
#SBATCH -N 1
#SBATCH --mem-per-cpu 100G


#SBATCH --time=01-00:00:00    ## "days-hours:minutes:seconds"

#SBATCH -e /trinity/home/couvreurs/ip-igx/slurm_logs/slurm-%j.err-%N
#SBATCH -o /trinity/home/couvreurs/ip-igx/slurm_logs/slurm-%j.out-%N

module purge > /dev/null 2>&1
module load R

cd /trinity/home/couvreurs/ip-igx/
Rscript src/ip_cor_density_batch.R
