#!/bin/bash
#SBATCH -J Pred_Asso
#SBATCH -n 1
#SBATCH --error=./log_sbatch/Asso_%J.err
#SBATCH --output=./log_sbatch/Asso_%J.out

module purge
module load gcc-10.2.0
module load plink-1.9
source ~/.bashrc

main_dir=6_Pred_TWAS/PrediXcan

plink_file=0_ref/freeze_2_auto_phased_final_pass_lung_tissue_normal_ID

expr_file=6_Pred_TWAS/pdui.peer.residuals_N.dna.txt

gencode_file=6_Pred_TWAS/3UTR_location.txt

db_path=6_Pred_TWAS/PrediXcan/output/lung_tissue_normal.db

cov_path=6_Pred_TWAS/PrediXcan/output/lung_tissue_normal.cov

gwas_path=0_ref/EAS_weightN_META_OUR_Discovery_Replication_BBJ_RESULT1_polish.txt

asso_out_path=6_Pred_TWAS/PrediXcan/results/lung_tissue_normal_asso.csv

Rscript ./7_PrediXcan_r.R \
        --asso_test \
        --main_dir ${main_dir} \
        --plink_file_name ${plink_file} \
        --expression_file_name ${expr_file} \
        --annotation_file_name ${gencode_file} \
        --db_path  ${db_path} \
        --cov_path ${cov_path} \
        --gwas_path  ${gwas_path} \
        --asso_out_path  ${asso_out_path} \
        --gwas_variant_col  MarkerName \
        --gwas_beta_col BETA \
        --gwas_se_col  SE \
        --gwas_eff_allele_col Allele1  \
        --gwas_ref_allele_col Allele2   > ./log_sbatch/lung_tissue_normal_asso.Rout

