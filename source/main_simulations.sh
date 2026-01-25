#!/bin/bash
#SBATCH --array=1-54%54
#SBATCH --job-name=wqs2_qg_1
#SBATCH --partition=wrobel
#SBATCH --output=wqs2_qg_1.out
#SBATCH --error=wqs2_qg_1.err

module purge
module load R/4.4.0

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript main_simulations.R $JOBID


