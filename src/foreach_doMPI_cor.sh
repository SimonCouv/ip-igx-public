#!/usr/bin/bash

# Number of MPI tasks
#SBATCH -n 60
##SBATCH -N 4
#
# Number of tasks per node
##SBATCH --tasks-per-node=12
#
# Runtime of this job
#SBATCH --time=07:00:00
# Stderr and stdout
#SBATCH -e /trinity/home/couvreurs/ip-igx/slurm_logs/slurm-%j.err-%N
#SBATCH -o /trinity/home/couvreurs/ip-igx/slurm_logs/slurm-%j.out-%N

# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1

# Load the module environment suitable for the job
module load R
module load openmpi

cd /trinity/home/couvreurs/ip-igx/src
# And finally run the job​
mpirun -np 60 Rscript mpi_blockwise_cor.R ../data/ips_trans_fam_adj.tsv ../data/ips_slurm_blockwise_cor.tsv pearson



# run as:
# sbatch foreach_doMPI_scratch.sh
