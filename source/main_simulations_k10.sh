#!/bin/bash
#SBATCH --array=1-12%12
#SBATCH --job-name=k10_20_0
#SBATCH --partition=wrobel
#SBATCH --output=k10_20_0.out
#SBATCH --error=k10_20_0.err

module purge
module load R/4.4.0

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript main_simulations_k10.R $JOBID


