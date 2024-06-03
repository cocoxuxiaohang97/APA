#!/bin/bash
#SBATCH -J TrainModel
#SBATCH -n 1
#SBATCH --array=1-195
#SBATCH --error=./log_sbatch/train_normal_%A_%a.err
#SBATCH --output=./log_sbatch/train_normal_%A_%a.out

module purge
module load gcc-10.2.0
module load plink-1.9
source ~/.bashrc

main_dir=6_Pred_TWAS/PrediXcan

plink_file=0_ref/freeze_2_auto_phased_final_pass_lung_tissue_normal_ID

expr_file=6_Pred_TWAS/pdui.peer.residuals_N.dna.txt

gencode_file=6_Pred_TWAS/3UTR_location.txt

Rscript ./7_PrediXcan_r.R \
        --model_training \
        --main_dir ${main_dir} \
        --plink_file_name ${plink_file} \
        --expression_file_name ${expr_file} \
        --subjob_id ${SLURM_ARRAY_TASK_ID} \
        --annotation_file_name ${gencode_file} >./log_sbatch/train_normal${SLURM_ARRAY_TASK_ID}.Rout

