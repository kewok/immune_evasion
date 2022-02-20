#!/bin/bash 
# Submit with sbatch --array=1-330400 simulateModel.sh
#SBATCH --time=96:00:00
#SBATCH --partition=amdsmall
#SBATCH --cpus-per-task=1
#SBATCH --mem=100gb
module load R/4.1.0
cd DIRECTORY_WHERE_spawn_folders_300.R_IS_STORED
Rscript spawn_folders_300.R $SLURM_ARRAY_TASK_ID
