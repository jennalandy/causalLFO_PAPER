#!/bin/bash

#SBATCH --job-name=stabilize_boot_0.8
#SBATCH --time=50:00:00       
#SBATCH --ntasks=1
#SBATCH --mem=5G      
#SBATCH --error=logs/stabilize_boot_0.8.error
#SBATCH --output=logs/stabilize_boot_0.8.log
#SBATCH -p shared

# echo "SLURM_JOBID: $SLURM_JOB_ID"
# echo "Running on host: $(hostname)"
# echo "Working directory: $(pwd)"

module load R/4.2.2-fasrc01
Rscript stabilize_boot_0.8.R