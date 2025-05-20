#!/bin/bash

#SBATCH --job-name=IAPOs_all_data
#SBATCH --time=10:00:00       
#SBATCH --ntasks=1
#SBATCH --mem=1G      
#SBATCH --error=logs/IAPOs_all_data.error
#SBATCH --output=logs/IAPOs_all_data.out
#SBATCH --array=1-100
#SBATCH -p shared

module load R/4.2.2-fasrc01
Rscript all_data.R $SLURM_ARRAY_TASK_ID