#!/bin/bash
#SBATCH -J Pred_DB
#SBATCH -n 1
#SBATCH --error=./log_sbatch/DB_%J.err
#SBATCH --output=./log_sbatch/DB_%J.out

module purge
module load gcc-10.2.0
module load plink-1.9
source ~/.bashrc

main_dir=6_Pred_TWAS/PrediXcan

plink_file=0_ref/freeze_2_auto_phased_final_pass_lung_tissue_normal_ID

expr_file=6_Pred_TWAS/pdui.peer.residuals_N.dna.txt

gencode_file=6_Pred_TWAS/3UTR_location.txt


Rscript ./7_PrediXcan_r.R \
        --generate_db_and_cov \
        --main_dir ${main_dir} \
        --plink_file_name ${plink_file} \
        --expression_file_name ${expr_file} \
        --annotation_file_name ${gencode_file} \
	--output_file_name lung_tissue_normal > ./log_sbatch/lung_tissue_normal_db.Rout

