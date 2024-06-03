rm(list = ls())
# -------------------- prepare covariates matrix T
load("24curate_covs_T.RData")
library(dplyr)
factors.df <- cbind(rownames(factors),factors)
colnames(factors.df) <- c("id",colnames(factors))
factors.df <- as.data.frame(factors.df)
factors.reorder <- factors.df %>% dplyr::select("id",all_of(id_order))
covariate_file <- "../4_aQTL_mapping/Covariate_matrix_T.txt"
write.table(factors.reorder, file=covariate_file, row.names=FALSE,quote=FALSE, sep='\t',col.names=T)
# -------------------- prepare covariates matrixN
load("24curate_covs_N.RData")
library(dplyr)
factors.df <- cbind(rownames(factors),factors)
colnames(factors.df) <- c("id",colnames(factors))
factors.df <- as.data.frame(factors.df)
factors.reorder <- factors.df %>% dplyr::select("id",all_of(id_order))
covariate_file <- "../4_aQTL_mapping/Covariate_matrix_N.txt"
write.table(factors.reorder, file=covariate_file, row.names=FALSE,quote=FALSE, sep='\t',col.names=T)

#!/opt/app/languages/R-3.6.3/bin/Rscript
rm(list=ls())
library(optparse)
library(MatrixEQTL)

option_list <- list(
  make_option(c("-p","--phenotype"),type="character",default="../4_aQTL_mapping/Phenotype_matrix_N.txt",action="store",help="APA expression data for MatrixEQTL"),
  make_option(c("-g","--genotype"),type="character",default="../4_aQTL_mapping/Genotype_matrix_N.txt",action="store",help="Genotype data for MatrixEQTL"),
  make_option(c("-c","--covariate"),type="character",default="../4_aQTL_mapping/Covariate_matrix_N.txt",action="store",help="Covariates for MatrixEQTL"),
  make_option(c("-s","--snp_location"),type="character",default="../4_aQTL_mapping/snp_location.txt",action="store",help="SNP locations"),
  make_option(c("-u","--utr_location"),type="character",default="../4_aQTL_mapping/3utr_location.txt",action="store",help="3UTR locations"),
  make_option(c("-w","--window"),type="numeric",default=1e6,action="store",help="window size"),
  make_option(c("-q","--cis_pvalue"),type="numeric",default=1e-2,action="store",help="p value threshold for cis-3'aQTL"),
  make_option(c("-Q","--trans_pvalue"),type="numeric",default=1e-5,action="store",help="p value trheshold for trans-3'aQTL")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))
PHENO <- opt$phenotype
GENO <- opt$genotype
COVARIATE <- opt$covariate
SNPLOC <- opt$snp_location
UTRLOC <- opt$utr_location
CIS_DISTANCE <- as.numeric(opt$window)
CIS_P_CUTOFF <- as.numeric(opt$cis_pvalue)
TRANS_P_CUTOFF <- as.numeric(opt$trans_pvalue)

cat('Options:\n','Phenotype:',PHENO,'\n','Genotype:',GENO,'\n','Covariates:',COVARIATE,'\n','CIS_DISTANCE:',CIS_DISTANCE,'\n',
    'CIS_P_CUTOFF:',CIS_P_CUTOFF,'\n','TRANS_P_CUTOFF:',TRANS_P_CUTOFF,'\n')
# - Use linear model
useModel = modelLINEAR # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
# - Genotype file name
SNP_file_name = GENO
snps_location_file_name = SNPLOC
# - APA expression file name
expression_file_name = PHENO
gene_location_file_name = UTRLOC
# - Covariates file name
covariates_file_name = COVARIATE
# - output file name
output_file_name_cis = "../4_aQTL_mapping/Cis_3aQTL_all_control_gene_exprs_N.txt"
output_file_name_tra = "../4_aQTL_mapping/Trans_3aQTL_all_control_gene_exprs_N.txt"
output_figure_name_cis = "../4_aQTL_mapping/Cis_3aQTL_genotype_info_control_gene_exprs_N.pdf"
pdf(output_figure_name_cis)
# - threshold
pvOutputThreshold_cis = CIS_P_CUTOFF;
pvOutputThreshold_tra = TRANS_P_CUTOFF;
# - Error covariance matrix
# set to numeric() for identity
errorCovariance = numeric();
# - Distance for local gene-SNP pairs
cisDist = CIS_DISTANCE;

# -- load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";
snps$fileOmitCharacters = "NA";
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

# -- load apa expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

# -- load covariates data
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE);
gz1 <- "25data.RDataw2_N"
save.image(gz1)
#unlink(output_file_name_tra);
#unlink(output_file_name_cis);
# -- Results
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')
cat('Detected local aQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant aQTLs:', '\n');
show(me$trans$eqtls)
write.table(me$cis$min.pv.gene, "../4_aQTL_mapping/cis.min.pv.gene_N.txt")
write.table(me$trans$min.pv.gene, "../4_aQTL_mapping/trans.min.pv.gene_N.txt")
save(gene,snps,file="../4_aQTL_mapping/Gene_SNP_N.RData")
save(snps,genepos,file="../4_aQTL_mapping/permutation_N.RData")
# -- plot the Q-Q plot of local and distant p-values
plot(me)
dev.off()

rm(list=ls())
library(optparse)
library(MatrixEQTL)
option_list <- list(
  make_option(c("-p","--phenotype"),type="character",default="../4_aQTL_mapping/Phenotype_matrix_T.txt",action="store",help="APA expression data for MatrixEQTL"),
  make_option(c("-g","--genotype"),type="character",default="../4_aQTL_mapping/Genotype_matrix_T.txt",action="store",help="Genotype data for MatrixEQTL"),
  make_option(c("-c","--covariate"),type="character",default="../4_aQTL_mapping/Covariate_matrix_T.txt",action="store",help="Covariates for MatrixEQTL"),
  make_option(c("-s","--snp_location"),type="character",default="../4_aQTL_mapping/snp_location.txt",action="store",help="SNP locations"),
  make_option(c("-u","--utr_location"),type="character",default="../4_aQTL_mapping/3utr_location.txt",action="store",help="3UTR locations"),
  make_option(c("-w","--window"),type="numeric",default=1e6,action="store",help="window size"),
  make_option(c("-q","--cis_pvalue"),type="numeric",default=1e-2,action="store",help="p value threshold for cis-3'aQTL"),
  make_option(c("-Q","--trans_pvalue"),type="numeric",default=1e-5,action="store",help="p value trheshold for trans-3'aQTL")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))
PHENO <- opt$phenotype
GENO <- opt$genotype
COVARIATE <- opt$covariate
SNPLOC <- opt$snp_location
UTRLOC <- opt$utr_location
CIS_DISTANCE <- as.numeric(opt$window)
CIS_P_CUTOFF <- as.numeric(opt$cis_pvalue)
TRANS_P_CUTOFF <- as.numeric(opt$trans_pvalue)

cat('Options:\n','Phenotype:',PHENO,'\n','Genotype:',GENO,'\n','Covariates:',COVARIATE,'\n','CIS_DISTANCE:',CIS_DISTANCE,'\n',
    'CIS_P_CUTOFF:',CIS_P_CUTOFF,'\n','TRANS_P_CUTOFF:',TRANS_P_CUTOFF,'\n')
# - Use linear model
useModel = modelLINEAR # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
# - Genotype file name
SNP_file_name = GENO
snps_location_file_name = SNPLOC
# - APA expression file name
expression_file_name = PHENO
gene_location_file_name = UTRLOC
# - Covariates file name
covariates_file_name = COVARIATE
# - output file name
output_file_name_cis = "../4_aQTL_mapping/Cis_3aQTL_all_control_gene_exprs_T.txt"
output_file_name_tra = "../4_aQTL_mapping/Trans_3aQTL_all_control_gene_exprs_T.txt"
output_figure_name_cis = "../4_aQTL_mapping/Cis_3aQTL_genotype_info_control_gene_exprs_T.pdf"
pdf(output_figure_name_cis)
# - threshold
pvOutputThreshold_cis = CIS_P_CUTOFF;
pvOutputThreshold_tra = TRANS_P_CUTOFF;
# - Error covariance matrix
# set to numeric() for identity
errorCovariance = numeric();
# - Distance for local gene-SNP pairs
cisDist = CIS_DISTANCE;

# -- load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";
snps$fileOmitCharacters = "NA";
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

# -- load apa expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

# -- load covariates data
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE);
gz1 <- "25data.RDataw2_T"
save.image(gz1)
#unlink(output_file_name_tra);
#unlink(output_file_name_cis);
# -- Results
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')
cat('Detected local aQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant aQTLs:', '\n');
show(me$trans$eqtls)
write.table(me$cis$min.pv.gene, "../4_aQTL_mapping/cis.min.pv.gene_T.txt")
write.table(me$trans$min.pv.gene, "../4_aQTL_mapping/trans.min.pv.gene_T.txt")
save(gene,snps,file="../4_aQTL_mapping/Gene_SNP_T.RData")
save(snps,genepos,file="../4_aQTL_mapping/permutation_T.RData")
# -- plot the Q-Q plot of local and distant p-values
plot(me)
dev.off()
