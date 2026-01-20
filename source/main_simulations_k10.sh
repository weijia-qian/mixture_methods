#!/bin/bash
#SBATCH --array=1-6%6
#SBATCH --job-name=k10_2
#SBATCH --partition=week-long-cpu
#SBATCH --output=main_simulations_k10_2.out
#SBATCH --error=main_simulations_k10_2.err

module purge
module load R/4.4.0

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript main_simulations_k10.R $JOBID


