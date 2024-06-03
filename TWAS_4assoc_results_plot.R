rm(list = ls())
library(readr)
library(data.table)
normal <- as.data.frame(fread("results/lung_tissue_normal_asso.csv"))
normal$`type`='normal'

gene_pos<-as.data.frame(fread("6_Pred_TWAS/3UTR_location.txt"))
#gene_pos<-gene_pos[which(gene_pos$genetype %in% c('lncRNA','protein_coding','miRNA')),]
gene_pos<-gene_pos[,-6]

normal <- merge(normal,gene_pos,by.x = "gene",by.y = "geneid")
pooled_pfdr_within_all <- normal
pooled_pfdr_within_all $`p_fdr` <- p.adjust(pooled_pfdr_within_all $pvalue, method = "fdr")
pooled_pfdr_within_all $`p_bonf` <- p.adjust(pooled_pfdr_within_all $pvalue, method = "bonferroni")
colnames(pooled_pfdr_within_all)[11] <- "pname"
write.table(pooled_pfdr_within_all, file="pooled_pfdr_within_all_normal.txt",quote = F, row.names = F, sep = '\t')
pooled <- pooled_pfdr_within_all
pooled$mergylabel <- paste(pooled$gene, pooled$genename, sep = "|")

library(data.table)
require(dplyr)
library(foreach)   #new added
library(doParallel)   #new added
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
options(warn = -1)
#data_list <- lapply(lf, read.table,header= T)

manhattan_fun = function(group, title = '', i_tag = 1){

  df = pooled[,c('genename','pvalue','p_bonf','chr','left','right','p_fdr','zscore','pname','mergylabel')]
  #df$dir = df[,paste0(group,'.beta')] * df[,'zscore']
  df$label = NA
  label_idx = which(df$p_fdr < 0.1)
  #label_idx = intersect(label_idx, which(df[,paste0(group, '.p')] < 0.05))
  df$label=as.character(df$label)
  df[label_idx,'label'] <- df[label_idx, 'pname']

  df$plog = -log(df$pvalue,10)
  df$chr = as.numeric(sub('chr','',df$chr))

  #pch
  df$pch = 'cir'
  df[which(df$p_fdr<0.1),'pch'] = 'dim'
  df = df[order(df$pch, decreasing = T),]

  #pos
  df$pos=(df$left+df$right)/2
  chr_box<-as.data.frame(matrix(data=NA,nrow=22,ncol=3))
  chr_df<-df[df$chr==1,]
  chr_box[1,1]<-max(chr_df$pos)
  chr_box[1,2]<-0
  chr_box[1,3]<-chr_box[1,1]/2

  for (i in 2:22){
    chr_df<-df[df$chr==i,]
    chr_box[i,1]<-max(chr_df$pos)
    chr_box[i,2]<-chr_box[i-1,1]+chr_box[i-1,2]
    df[which(df$chr==i),'pos']<-df[which(df$chr==i),'pos']+chr_box[i,2]
    chr_box[i,3]<-chr_box[i,2]+chr_box[i,1]/2}

  df$group<-ifelse(df$chr%%2==0,'odd','even')

  set.seed(1)
  manhattan <- ggplot(df, aes(x=pos, y=zscore, label = label)) +
    scale_x_continuous(breaks=chr_box$V3, labels = c(as.character(seq(1,19)),' ','21',' ')) +  #x axis 
    geom_point(data = df, aes(x = pos, y = zscore,
                              color = group,
                              alpha = pch,
                              shape = pch,
                              size = pch), show.legend = FALSE) + #points
    scale_color_manual(values = c("odd" = 'firebrick', "even" = 'dodgerblue4')) +
    scale_shape_manual(values = c("cir" = 20, 'dim' = 5))+
    scale_alpha_manual(values = c("cir" = 0.3, 'dim' = 0.95))+
    scale_size_manual(values = c("cir" = 7, 'dim' = 6))+

    labs(x = "Chromosome", y = "zscore") +
    theme(legend.position="none") +
    theme_bw() +  
    theme(panel.grid =element_blank(),panel.border = element_rect(color = "black", size = 0.1, fill = NA))+ 
    theme(axis.title.x = element_text(size = 20, vjust = 0.5, hjust = 0.5),axis.text.x = element_text(size = 20, vjust = 0.5, hjust = 0.5))+
    theme(axis.title.y = element_text(size = 20, vjust = 0.5, hjust = 0.5),axis.text.y = element_text(size = 20, vjust = 0.5, hjust = 0.5))+
    scale_y_continuous(
                       limits = c(min(df$zscore),max(df$zscore))) +  #y axis log
    # labs(tag = letters[i_tag])+
    ggtitle(title) +theme(plot.title = element_text(size = 20, vjust = 0, hjust = 0))+ geom_hline(aes(yintercept=0),  size=1)+

    geom_label_repel( #non overlapped labels
      size=6,
      #color = df$color,
      #nudge_x=5e4, #shift to the righ                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     t
      segment.alpha = 0.2,  #transparent of segment
      min.segment.length = 1,
      segment.size = 0.2,
      fill=rgb(255, 255, 255, 210, maxColorValue=255),
      fontface='italic',
      seed = 2022
    )
  return(manhattan)
}

pdf('pooled_pfdr_within_all_normal_zscore.pdf', width = 16, height = 12)
manhattan_fun(group = '', title = 'pooled_pfdr_within_all_normal',i_tag = 1)
dev.off()

rm(list = ls())
setwd("6_Pred_TWAS/PrediXcan")
pooled <- read.table("pooled_pfdr_within_all_normal.txt", header=T,sep="\t",stringsAsFactors=F, check.names =FALSE)
pooled_sig_bof <- pooled[pooled$p_bonf<0.1,]
pooled_sig_fdr <- pooled[pooled$p_fdr<0.1,]
pooled_sig_p <- pooled[pooled$pvalue<0.05,]
save(pooled_sig_bof, pooled_sig_fdr, pooled_sig_p, pooled, file = "pooled_sig_within_all_normal.RData")
write.table(pooled_sig_p, file="pooled_pfdr_within_all_normal_p0.05.txt",quote = F, row.names = F, sep = '\t')

load("pooled_sig_within_all_normal.RData")

rm(list = ls())
library(readr)
library(data.table)
normal <- as.data.frame(fread("results/TWAS_EAS_Meta_based_lung_tissue_normal.txt"))
normal$`type`='normal'

gene_pos<-as.data.frame(fread("0_ref/gene_location_gencode_v32_GRCh38.txt"))

normal <- merge(normal,gene_pos,by.x = "gene",by.y = "geneid")
pooled_pfdr_within_all <- normal
pooled_pfdr_within_all $`p_fdr` <- p.adjust(pooled_pfdr_within_all $pvalue, method = "fdr")
pooled_pfdr_within_all $`p_bonf` <- p.adjust(pooled_pfdr_within_all $pvalue, method = "bonferroni")
colnames(pooled_pfdr_within_all)[2] <- "pname"
write.table(pooled_pfdr_within_all, file="pooled_pfdr_within_all_normal_exp.txt",quote = F, row.names = F, sep = '\t')
pooled <- pooled_pfdr_within_all


library(data.table)
require(dplyr)
library(foreach)   #new added
library(doParallel)   #new added
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
options(warn = -1)
#data_list <- lapply(lf, read.table,header= T)

manhattan_fun = function(group, title = '', i_tag = 1){
  
  df = pooled[,c('pname','pvalue','p_bonf','chr','left','right','p_fdr','zscore')]
  #df$dir = df[,paste0(group,'.beta')] * df[,'zscore']
  df$label = NA
  label_idx = which(df$p_fdr < 0.1)
  #label_idx = intersect(label_idx, which(df[,paste0(group, '.p')] < 0.05))
  df$label=as.character(df$label)
  df[label_idx,'label'] <- df[label_idx, 'pname']
  
  df$plog = -log(df$pvalue,10)
  df$chr = as.numeric(sub('chr','',df$chr))
  
  #pch
  df$pch = 'cir'
  df[which(df$p_fdr<0.1),'pch'] = 'dim'
  df = df[order(df$pch, decreasing = T),]
  
  #pos
  df$pos=(df$left+df$right)/2
  chr_box<-as.data.frame(matrix(data=NA,nrow=22,ncol=3))
  chr_df<-df[df$chr==1,]
  chr_box[1,1]<-max(chr_df$pos)
  chr_box[1,2]<-0
  chr_box[1,3]<-chr_box[1,1]/2
  
  for (i in 2:22){
    chr_df<-df[df$chr==i,]
    chr_box[i,1]<-max(chr_df$pos)
    chr_box[i,2]<-chr_box[i-1,1]+chr_box[i-1,2]
    df[which(df$chr==i),'pos']<-df[which(df$chr==i),'pos']+chr_box[i,2]
    chr_box[i,3]<-chr_box[i,2]+chr_box[i,1]/2}
  
  df$group<-ifelse(df$chr%%2==0,'odd','even')
  
  set.seed(1)
  manhattan <- ggplot(df, aes(x=pos, y=zscore, label = label)) +
    scale_x_continuous(breaks=chr_box$V3, labels = c(as.character(seq(1,19)),' ','21',' ')) +  #x axis 
    geom_point(data = df, aes(x = pos, y = zscore,
                              color = group,
                              alpha = pch,
                              shape = pch,
                              size = pch), show.legend = FALSE) + #points
    scale_color_manual(values = c("odd" = 'firebrick', "even" = 'dodgerblue4')) +
    scale_shape_manual(values = c("cir" = 20, 'dim' = 5))+
    scale_alpha_manual(values = c("cir" = 0.3, 'dim' = 0.95))+
    scale_size_manual(values = c("cir" = 7, 'dim' = 6))+
    
    labs(x = "Chromosome", y = "zscore") +
    theme(legend.position="none") +
    theme_bw() +  
    theme(panel.grid =element_blank(),panel.border = element_rect(color = "black", size = 0.1, fill = NA))+ #去除网格线
    theme(axis.title.x = element_text(size = 20, vjust = 0.5, hjust = 0.5),axis.text.x = element_text(size = 20, vjust = 0.5, hjust = 0.5))+
    theme(axis.title.y = element_text(size = 20, vjust = 0.5, hjust = 0.5),axis.text.y = element_text(size = 20, vjust = 0.5, hjust = 0.5))+
    scale_y_continuous(
      limits = c(min(df$zscore),max(df$zscore))) +  #y axis log
    # labs(tag = letters[i_tag])+
    ggtitle(title) +theme(plot.title = element_text(size = 20, vjust = 0, hjust = 0))+ geom_hline(aes(yintercept=0),  size=1)+
    
    geom_label_repel( #non overlapped labels
      size=6,
      #color = df$color,
      #nudge_x=5e4, #shift to the righ                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     t
      segment.alpha = 0.2,  #transparent of segment
      min.segment.length = 1,
      segment.size = 0.2,
      fill=rgb(255, 255, 255, 210, maxColorValue=255),
      fontface='italic',
      seed = 2022
    )
  return(manhattan)
}

  
pdf('pooled_pfdr_within_all_normal_exp_zscore.pdf', width = 16, height = 12)
manhattan_fun(group = '', title = 'pooled_pfdr_within_all_normal',i_tag = 1)
dev.off()
