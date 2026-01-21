#!/bin/bash
#SBATCH --array=1-6%6
#SBATCH --job-name=n2000_1
#SBATCH --partition=week-long-cpu
#SBATCH --output=main_simulations_n2000_1.out
#SBATCH --error=main_simulations_n2000_1.err

module purge
module load R/4.4.0

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript main_simulations_n2000.R $JOBID


