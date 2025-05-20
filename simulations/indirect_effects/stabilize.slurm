#!/bin/bash

#SBATCH --job-name=IAPOs_stabilize
#SBATCH --time=10:00:00       
#SBATCH --ntasks=1
#SBATCH --mem=1G      
#SBATCH --error=logs/IAPOs_stabilize.error
#SBATCH --output=logs/IAPOs_stabilize.out
#SBATCH --array=1-100
#SBATCH -p shared

module load R/4.2.2-fasrc01
Rscript stabilize.R $SLURM_ARRAY_TASK_ID