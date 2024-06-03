rm(list = ls())
ExM <- read.table("Phenotype_matrix_pdui_747.txt",header=T,sep="\t",stringsAsFactors=F, check.names =FALSE)
rownames(ExM) <- ExM[,1]
ExM <- ExM[,-1]
ID <- read.table("finalsample1108.csv",stringsAsFactors=F, check.names =FALSE,sep=",",header=T)
sT <-  ID$ID[ID$group=="T"]
sN <- ID$ID[ID$group=="N"]
tumor <- signif(apply(ExM[rownames(ExM), sT], 1, mean), 4)
normal <- signif(apply(ExM[rownames(ExM), sN], 1, mean), 4)
PDUI <- tumor - normal
rPDUI = data.frame(tumor, normal, PDUI)
cat("mark change\n") 
rPDUI$change <- rep("NC", nrow(rPDUI))
rPDUI$change[rPDUI$PDUI > 0.1] = "lengthened_in_tumor"
rPDUI$change[rPDUI$PDUI < -0.1] = "shortened_in_tumor"
library(ggplot2)
hist(rPDUI$PDUI)
p1 <- ggplot(data = rPDUI,aes(x=tumor,y=normal,color=change))+
  geom_point(alpha=1, size=1.5) + 
  theme_set(theme_set(theme_bw(base_size=20))) + theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border = element_rect(color = "black", size = 0.1, fill = NA)) + 
  theme(axis.title.x = element_text(size = 15, vjust = 0.5, hjust = 0.5, colour = "black"),axis.text.x = element_text(size = 15, vjust = 0.5, hjust = 0.5, colour = "black"))+
  theme(axis.title.y = element_text(size = 15, vjust = 0.5, hjust = 0.5, colour = "black"),axis.text.y = element_text(size = 15, vjust = 0.5, hjust = 0.5, colour = "black"))+
  theme(legend.text=element_text(size=15))+
  xlab("PDUI_Tumor") + ylab("PDUI_Non_Tumor") + geom_abline(intercept = 0.1, slope = 1, linetype = "dashed", size = 0.3) + geom_abline(intercept = -0.1, slope = 1, linetype = "dashed", size = 0.3) +
  scale_colour_manual(values = c(lengthened_in_tumor='firebrick',NC='grey',shortened_in_tumor='dodgerblue4')) 
ggsave(plot = p1, '../8_diff_TNOUTCOME/80_diff_pdui.pdf',height = 8,width = 10)
rm(list = ls())
ExM <- read.table("Phenotype_matrix_pdui_747.txt",header=T,sep="\t",stringsAsFactors=F, check.names =FALSE)
rownames(ExM) <- ExM[,1]
ExM <- ExM[,-1]
log2ExM <- log2(ExM+1)
log2ExM <- na.omit(log2ExM)
OUT <- read.table("finalsample1108.csv",stringsAsFactors=F, check.names =FALSE,sep=",",header=T)
sT <-  OUT$ID[OUT$group=="T"] 
sN <- OUT$ID[OUT$group=="N"]
save(log2ExM, OUT,file = "../8_diff_TNOUTCOME/80_log2Pheno_ID_pdui.RData")

# cat("wilcox.test\n")
PValue = FDR = logFC = matrix(0, nrow(log2ExM), 1)
for(i in 1:nrow(log2ExM)){
  PValue[i, 1] = p.value = wilcox.test(as.numeric(log2ExM[i, sT]), as.numeric(log2ExM[i, sN]))$p.value
  logFC[i, 1] = mean(as.numeric(log2ExM[i, sT]))-mean(as.numeric(log2ExM[i, sN]))
}
gene = rownames(log2ExM)
rTable = data.frame(gene, logFC, PValue)
# rTable2 = rTable[order(rTable$PValue), ]
# PValue2 <- rTable2$PValue
# rTable2$FDR = p.adjust(as.vector(PValue2), "BH", n = length(PValue2))
FDR = p.adjust(as.vector(rTable$PValue), "fdr")
rTable = data.frame(logFC, PValue, FDR, row.names = rownames(log2ExM))
tumor <- signif(apply(log2ExM[rownames(rTable), sT], 1, mean), 4)
normal <- signif(apply(log2ExM[rownames(rTable), sN], 1, mean), 4)
cat("mark change\n") 
change <- rep("NC", nrow(log2ExM))
logFC_cutoff_exactT <- with(rTable, mean(abs(logFC))+2*sd(abs(logFC))) ##0.07695772
change[((rTable$FDR) < 0.05) & (rTable$logFC > logFC_cutoff_exactT)] = "lengthened_in_tumor"
change[((rTable$FDR) < 0.05) & (rTable$logFC < -logFC_cutoff_exactT)] = "shortened_in_tumor"

rTable = data.frame(gene, tumor, normal, rTable[, c("logFC", "PValue", "FDR")], change)
table(rTable$change)
write.table(rTable,file="../8_diff_TNOUTCOME/81_rTable_diff_utr_pdui.txt",row.names=F,col.names=T,quote=F,sep="\t")
