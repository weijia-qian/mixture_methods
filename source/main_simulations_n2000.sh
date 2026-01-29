#!/bin/bash
#SBATCH --array=1-18%18
#SBATCH --job-name=n2000_4
#SBATCH --partition=chang
#SBATCH --output=n2000_4.out
#SBATCH --error=n2000_4.err

module purge
module load R/4.4.0

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript main_simulations_n2000.R $JOBID


