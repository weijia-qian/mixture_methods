#!/bin/bash
#SBATCH --array=1-9%9
#SBATCH --job-name=inter_6
#SBATCH --partition=chang
#SBATCH --output=inter_6.out
#SBATCH --error=inter_6.err

module purge
module load R/4.4.0

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript main_simulations.R $JOBID


