#!/bin/bash

#SBATCH --job-name=application
#SBATCH --time=2:00:00       
#SBATCH --ntasks=1
#SBATCH --mem=1G      
#SBATCH --error=logs/application.error
#SBATCH --output=logs/application.out
#SBATCH --array=2,7,12,17 # random split: 2,7,12,17, stabilize: 4,9,14,19
#SBATCH -p shared

module load R/4.2.2-fasrc01
Rscript application.R $SLURM_ARRAY_TASK_ID