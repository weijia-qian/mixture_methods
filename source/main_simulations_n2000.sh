#!/bin/bash
#SBATCH --array=1-6%6
#SBATCH --job-name=n2000_2
#SBATCH --partition=wrobel
#SBATCH --output=n2000_2.out
#SBATCH --error=n2000_2.err

module purge
module load R/4.4.0

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript main_simulations_n2000.R $JOBID


