##keeping genes with an expression value of at least 1 TPM in at least 20% (30/156) of tumour samples in the TRACERx multi-region RNAseq dataset.
##In total, 16,286 genes were filtered out of the 25,343 unique genes outputted by RSEM.
keep <- rowSums(tpm_raw>=1) >= floor(0.20*ncol(tpm_raw))
table(keep)
##22539 16687 
filter_tpm_raw <- tpm_raw[keep,]
filter_tpm_raw[1:4,1:4]
dim(filter_tpm_raw)
filter_tpm_adj <- tpm_adj[keep,]
filter_tpm_raw[1:4,1:4]
dim(filter_tpm_adj)
filter_counts <- counts[keep,]
filter_counts[1:4,1:4]
dim(filter_counts)
save(filter_counts, file = "R_results/01filterCount_join747.Rdata")
save(filter_tpm_adj, file = "R_results/01filterAdjTPM_join747.Rdata")
save(filter_tpm_raw, file = "R_results/01filterRawTPM_join747.Rdata")

##---------------------
#RawCounts_Combatseq_batcheffects
library(sva)
count_matrix <- as.matrix(counts)
batch <- as.character(group$batch)
group_sva <- group$group
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group_sva)
load("R_results/02adjusted_counts.Rdata")
# counts to tpm
genelenth <- read.table("genelength.txt",
                        header = T,sep = "\t")
adjusted_counts <- as.data.frame(adjusted_counts)
adjusted_counts <- cbind(rownames(adjusted_counts),adjusted_counts)
colnames(adjusted_counts)[1] <- "GeneId"

mycounts <- inner_join(adjusted_counts, genelenth, by = "GeneId")
head(mycounts)
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
adjusted_counts <- adjusted_counts[,-1]
head(mycounts)

kb <- mycounts$GeneLength / 1000
head(kb)

countdata <- mycounts[,1:747]
head(countdata)
rpk <- countdata / kb
head(rpk)
adjusted_tpm <- t(t(rpk)/colSums(rpk) * 1000000)
adjusted_tpm <- as.data.frame(adjusted_tpm)
