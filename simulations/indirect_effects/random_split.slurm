#!/bin/bash

#SBATCH --job-name=IAPOs_random_split
#SBATCH --time=10:00:00       
#SBATCH --ntasks=1
#SBATCH --mem=1G      
#SBATCH --error=logs/IAPOs_random_split.error
#SBATCH --output=logs/IAPOs_random_split.out
#SBATCH --array=1-100
#SBATCH -p shared

module load R/4.2.2-fasrc01
Rscript random_split.R $SLURM_ARRAY_TASK_ID