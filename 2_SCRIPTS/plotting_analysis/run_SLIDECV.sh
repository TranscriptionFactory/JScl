#!/bin/bash
#SBATCH -t 3-00:00
#SBATCH --mail-user=aar126@pitt.edu
#SBATCH --mail-type=FAIL
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=150g
#SBATCH --cpus-per-task=16
#SBATCH --job-name=Allre

module load gcc/10.2.0
module load r/4.2.0

Rscript run_SLIDECV_LoSAI.R
