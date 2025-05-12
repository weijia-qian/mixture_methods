#!/bin/bash
#SBATCH --array=1-21%21
#SBATCH --job-name=simulations
#SBATCH --partition=wrobel
#SBATCH --output=main_simulations.out
#SBATCH --error=main_simulations.err

module purge
module load R/4.4.0

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript main_simulations.R $JOBID


