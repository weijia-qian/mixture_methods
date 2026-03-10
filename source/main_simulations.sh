#!/bin/bash
#SBATCH --array=1-18%18
#SBATCH --job-name=inter_0
#SBATCH --partition=wrobel
#SBATCH --output=inter_0.out
#SBATCH --error=inter_0.err

module purge
module load R/4.4.0

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript main_simulations.R $JOBID


