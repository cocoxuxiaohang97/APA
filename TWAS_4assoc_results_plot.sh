#!/bin/bash
#SBATCH -J assoc_plot
#SBATCH -n 1
#SBATCH --error=./log_sbatch/assoc_normal_%A_%a.err
#SBATCH --output=./log_sbatch/assoc_normal_%A_%a.out

module purge
module load gcc-10.2.0
module load plink-1.9
source ~/.bashrc

cd 6_Pred_TWAS/PrediXcan 
Rscript ./74assoc_results_plot.R  > ./log_sbatch/lung_tissue_normal_assoc_plot.Rout
