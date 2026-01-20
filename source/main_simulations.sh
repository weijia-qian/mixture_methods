#!/bin/bash
#SBATCH --array=1-54%54
#SBATCH --job-name=simulations_4
#SBATCH --partition=week-long-cpu
#SBATCH --output=main_simulations_4.out
#SBATCH --error=main_simulations_4.err

module purge
module load R/4.4.0

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript main_simulations.R $JOBID


