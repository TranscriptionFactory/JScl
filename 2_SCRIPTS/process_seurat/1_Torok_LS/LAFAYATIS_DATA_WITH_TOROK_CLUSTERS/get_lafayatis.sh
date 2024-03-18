#!/bin/bash
#SBATCH -t 3-00:00
#SBATCH --job-name=mall
# partition (queue) declaration
#SBATCH --mail-user=aar126@pitt.edu
#SBATCH --partition=smp
#SBATCH --cluster=smp
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=150g
#SBATCH --cpus-per-task=16

module load gcc/12.2.0
module load r/4.3.0


# module load gcc/10.2.0
# module load r/4.2.0

Rscript /ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/scripts/getData/get_lafayatis_expr_mat.R

