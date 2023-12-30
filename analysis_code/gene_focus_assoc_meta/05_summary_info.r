library(reshape2)
library(dplyr)

rm(list=ls())

wd="/stanley/huang_lab/home/ychen/proj-CNV/08CNV_meta_analysis/02gene_assoc_meta"
setwd(wd)

rs_del <- read.table("output/META_PGC_EAS.del_1.txt", header=T)
rs_del$CNV_type <- "del"
rs_dup <- read.table("output/META_PGC_EAS.dup_1.txt", header=T)
rs_dup$CNV_type <- "dup"
rs <- rbind(rs_del,rs_dup)
rs$key <- paste0(rs$MarkerName,rs$CNV_type)
rs_info <- read.table("output/gene_assoc_META.txt", header=T)
rs_info$key <- paste0(rs_info$META_GENE,rs_info$META_CNV_type)
rs <- merge(rs,rs_info,all.x=T)
rs$key <- paste0(rs$MarkerName,rs$CNV_type)

EAS_freq <- read.csv("input/EAS_freq_output.csv", header=T)
EAS_freq$File.Name <- gsub(".log","",EAS_freq$File.Name)
names(EAS_freq) <- paste0("EAS_",names(EAS_freq))
EAS_freq$key <- paste0(EAS_freq$EAS_File.Name,EAS_freq$EAS_CopyN) 

EUR_freq <- read.csv("input/EUR_freq_output.csv", header=T)
EUR_freq$File.Name <- gsub(".log","",EUR_freq$File.Name)
EUR_freq$CNV_type <- ifelse(EUR_freq$CopyN=="1","del","dup")
names(EUR_freq) <- paste0("EUR_",names(EUR_freq))
EUR_freq$key <- paste0(EUR_freq$EUR_File.Name,EUR_freq$EUR_CNV_type) 

EAS_ori_del <- read.table("input/EAS.del.assoc",header=T)
EAS_ori_del$CNV_type <- "del"
names(EAS_ori_del) <- paste0("EAS_",names(EAS_ori_del))
EAS_ori_dup <- read.table("input/EAS.dup.assoc",header=T)
EAS_ori_dup$CNV_type <- "dup"
names(EAS_ori_dup) <- paste0("EAS_",names(EAS_ori_dup))
EAS_ori <- rbind(EAS_ori_del,EAS_ori_dup)
EAS_ori$key <- paste0(EAS_ori$EAS_GENE, EAS_ori$EAS_CNV_type)

EUR_ori_del <- read.table("input/PGC.del.assoc",header=T)
EUR_ori_del$CNV_type <- "del"
names(EUR_ori_del) <- paste0("EUR_",names(EUR_ori_del))
EUR_ori_dup <- read.table("input/PGC.dup.assoc",header=T)
EUR_ori_dup$CNV_type <- "dup"
names(EUR_ori_dup) <- paste0("EUR_",names(EUR_ori_dup))
EUR_ori <- rbind(EUR_ori_del,EUR_ori_dup)
EUR_ori$key <- paste0(EUR_ori$EUR_GENE, EUR_ori$EUR_CNV_type)

result_ori <- rs %>%
  left_join(EAS_freq, by="key") %>%
  left_join(EUR_freq, by="key") %>%
  left_join(EAS_ori, by="key") %>%
  left_join(EUR_ori, by="key")

#result <- result_ori[,c("MarkerName", "locus", "qvalue", "adjust.p", "loci_fdr", "loci_FWER", "Effect", "StdErr", "P.value", "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal", "CNV_type", "EAS_Case", "EAS_Control", "EUR_Case", "EUR_Control", "EAS_BETA", "EAS_SE", "EAS_P", "EUR_BETA", "EUR_P", "EUR_SE")]
result <- result_ori[,c("META_GENE","locus","META_P","META_Zscore","META_CNV_type","start","end","seqname","gene_name", "EAS_Case", "EAS_Control", "EUR_Case", "EUR_Control", "EAS_BETA", "EAS_SE", "EAS_P", "EUR_BETA", "EUR_P", "EUR_SE")]
result$EAS_Case_freq <- ifelse(result$EAS_Case=="NA",0,as.numeric(result$EAS_Case)/15475)
result$EAS_Control_freq <- ifelse(result$EAS_Control=="NA",0,as.numeric(result$EAS_Control)/17163)
result$EUR_Case_freq <- ifelse(result$EUR_Case=="NA",0,as.numeric(result$EUR_Case)/17506)
result$EUR_Control_freq <- ifelse(result$EUR_Control=="NA",0,as.numeric(result$EUR_Control)/16751)

write.table(result,"output/gene_assoc_META_info.txt",row.names=F,col.names=T,sep='\t',quote=F)
