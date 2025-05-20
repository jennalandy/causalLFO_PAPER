#!/bin/bash

#SBATCH --job-name=all_data_boot_0.2
#SBATCH --time=50:00:00       
#SBATCH --ntasks=1
#SBATCH --mem=5G      
#SBATCH --error=logs/all_data_boot_0.2.error
#SBATCH --output=logs/all_data_boot_0.2.log
#SBATCH -p shared

# echo "SLURM_JOBID: $SLURM_JOB_ID"
# echo "Running on host: $(hostname)"
# echo "Working directory: $(pwd)"

module load R/4.2.2-fasrc01
Rscript all_data_boot_0.2.R