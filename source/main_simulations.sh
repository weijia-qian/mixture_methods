#!/bin/bash
#SBATCH --array=1-54%54
#SBATCH --job-name=wqs2_qg_3
#SBATCH --partition=week-long-cpu
#SBATCH --output=wqs2_qg_3.out
#SBATCH --error=wqs2_qg_3.err

module purge
module load R/4.4.0

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript main_simulations.R $JOBID


