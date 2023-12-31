library(ggplot2)
library(ggman)
library(tidyverse)
library(ggsci)
library(cowplot)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(reshape2)
library(gtools)

rm(list=ls())
wd="/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/01main_figure/figure3_gene_associate/"
setwd(wd)
set.seed(123)

#color Platte

cohort_colors <- c(
  EUR = "#BE1E2D",
  #EAS = "#0087FF",
  EAS = "#3182BD",
  'Deletion' = "#ff7f0e",
  'Duplication' = "#2ca02c",
  'Novel Regions' = "#7ec0ee",
  'All Regions' = "#A9A9A9")

#only run once
# gene_annot <- read.csv("/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/misc/gencode.v19.gene.annot.csv")
# cytoband <- read.table("/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/misc/cytoBand.hg19.txt.gz")
# names(cytoband) <- c("chr","start","end","pos","type")
# cytoband$name <- paste0(cytoband$chr,cytoband$pos)
# eas_association_all <- read.table("input/EAS.all.assoc",header=T)
# eas_association_all$CNV_type <-"all"
# eas_association_all_permute <- read.table("input/all_permute_thresh.txt",header=T,na.strings = NA)
# eas_association_all <- merge(eas_association_all,eas_association_all_permute,all.x=T)
# eas_association_all$permute_p <- 10^(-(eas_association_all$THRESHOLD))
# 
# eas_association_del <- read.table("input/EAS.del.assoc",header=T)
# eas_association_del$CNV_type <- "del"
# eas_association_del_permute <- read.delim("input/del_permute_thresh.txt",header=T,na.strings = NA)
# eas_association_del <- merge(eas_association_del,eas_association_del_permute,all.x=T)
# eas_association_del$permute_p <- 10^(-(eas_association_del$THRESHOLD))
# 
# eas_association_dup <- read.table("input/EAS.dup.assoc",header=T)
# eas_association_dup$CNV_type <- "dup"
# eas_association_dup_permute <- read.delim("input/dup_permute_thresh.txt",header=T,na.strings = NA)
# eas_association_dup <- merge(eas_association_dup,eas_association_dup_permute,all.x=T)
# eas_association_dup$permute_p <- 10^(-(eas_association_dup$THRESHOLD))
# 
# eas_association_all <- rbind(eas_association_all,eas_association_del,eas_association_dup)
# names(eas_association_all) <- paste0("EAS_",names(eas_association_all))
# eas_association_all <- merge(eas_association_all,gene_annot,by.x="EAS_GENE",by.y="gene_id")
# eas_association_all$locus <- NA
# 
# # Iterate over each gene and check if it is within any cytoband range
# for (i in 1:nrow(eas_association_all)) {
#   gene_chromosome <- eas_association_all$seqname[i]
#   gene_start <- eas_association_all$start[i]
#   gene_end <- eas_association_all$end[i]
# 
#   for (j in 1:nrow(cytoband)) {
#     cytoband_chromosome <- cytoband$chr[j]
#     cytoband_start <- cytoband$start[j]
#     cytoband_end <- cytoband$end[j]
# 
#     if (gene_chromosome == cytoband_chromosome && gene_start >= cytoband_start && gene_end <= cytoband_end) {
#       #print(paste("Gene", i, "is located within cytoband", cytoband$name[j]))
#       eas_association_all[i,"locus"] <- cytoband$name[j]
#       break  #
#     }
#   }
# }
# 
# res <- eas_association_all[,c("EAS_GENE","locus","EAS_P","EAS_permute_p","EAS_BETA","EAS_SE","EAS_CNV_type","start","end","seqname","gene_name")]
# 
# write.table(res, "output/gene_assoc_EAS.txt", col.names = T, row.names = F, quote = F, sep = '\t' )
# write.table(res, "input/gene_assoc_EAS.txt", col.names = T, row.names = F, quote = F, sep = '\t' )

eas_association_all <- read.table("input/gene_assoc_EAS.txt", header=T)
eas_association_all$BP <- eas_association_all$start + ((eas_association_all$end - eas_association_all$start)/2)
eas_association_all$length <- eas_association_all$end - eas_association_all$start
eas_association_all$CHR <- gsub("chr","",eas_association_all$seqname)
eas_association_all$CHR <- ifelse(eas_association_all$CHR=="X",23,as.numeric(eas_association_all$CHR))
eas_association_all_res <- eas_association_all

eas_association_all$CNV_type <- eas_association_all$EAS_CNV_type
eas_association_all$CNV_type <- gsub("all","All rCNVs",eas_association_all$CNV_type)
eas_association_all$CNV_type <- gsub("del","Deletion",eas_association_all$CNV_type)
eas_association_all$CNV_type <- gsub("dup","Duplication",eas_association_all$CNV_type)

#threshold_FWER_p <- (0.05/(dim(eas_association_all)[1]))
#threshold_FWER_p <- (0.05/(18531*2))
#threshold_permute_p <- quantile(na.omit(eas_association_all$EAS_permute_p), 0.05)
#threshold_permute_p <- 4.08488655632183e-05*0.97
threshold_permute_p <- 7.119180431199928e-05 * 0.97
#threshold_permute_p <- (0.05/(18531*2))
threshold_p <- threshold_permute_p

## PLOT
# PLOT deletion
# data.frame prepare
eas_association <- eas_association_all[eas_association_all$EAS_CNV_type=="del",]
#threshold_p <- (0.05/length(unique(eas_association$EAS_GENE)))
eas_association.sig <- eas_association[eas_association$EAS_P<threshold_p,]
#eas_association.sig <- eas_association[eas_association$loci_permute=="Y",]

#x axis index
eas_association <- eas_association[order(eas_association$BP),]
eas_association <- eas_association[mixedorder(eas_association$CHR),]
eas_association$index <- 1:nrow(eas_association)

# caculate max y value
ymax_val <- max(-log10(min(eas_association$EAS_P)))

#prepare label - p-threshold
thresh_lable <- data.frame(
  x = 0,
  y = (-log10(signif(threshold_permute_p)))+1,
  text = paste0(" Permutation threshold = ",signif(threshold_permute_p,digits = 2))
)

#thresh_lable <- rbind(thresh_permutation_lable, thresh_FWER_lable)

P1 <- ggman(eas_association, snp = "locus", bp = "start", chrom = "CHR", pvalue = "EAS_P",  pointSize = 0.6, cex.axis=0.8, ymax = ymax_val+2, sigLine = NA) + 
  theme_classic() + 
  #geom_hline(yintercept=-log10(signif(threshold_FWER_p, digits = 3)),linetype=2) +
  geom_hline(yintercept=-log10(signif(threshold_permute_p, digits = 3)),linetype=2) +
  xlab("Chromosome") + 
  ylab("-log10(p-value)") +
  theme(axis.text.x = element_text(angle = 45, size=8, hjust=1)) +  
  geom_text(data = thresh_lable, aes(x, y, label = text, col="black"), hjust = 0) +
  ggtitle("Deletion") + 
  scale_colour_manual(values = c("darkgrey","lightgrey","black","lightskyblue")) +
  scale_fill_manual(values = c("lightskyblue2","midnightblue","black","IndianRed1")) 

#prepare label- significant loci 
dfm <- P1[[1]]
dfm_novel <- dfm[dfm$locus %in% c("chr11p15.5", "chr11q13.1", "chr11q13.2", "chr17q25.3", "chr18q23", "chr19p13.3", "chr19q13.42", "chr5q31.3", "chr8p21.3", "chr8q24.3"),]
dfm_novel_sig <- dfm_novel[dfm_novel$EAS_P < threshold_p,]

filtered_data <- eas_association.sig %>% 
  group_by(locus) %>% 
  dplyr::filter(EAS_P== min(EAS_P))
filtered_data$snp <- filtered_data$locus
filtered_data <- filtered_data[!duplicated(filtered_data$locus),]
filtered_data_del <- filtered_data

labelDfm <- as.data.frame(filtered_data)
labelDfm$snp <- labelDfm[,"EAS_GENE"]
labelDfm$label <- labelDfm[,"EAS_GENE"]
dfm.sub <- dfm[dfm$EAS_GENE %in% labelDfm$EAS_GENE,]
dfm.sub$label <- paste0(dfm.sub$locus)
dfm.label.sub <- dfm.sub[!duplicated(dfm.sub$locus),]
dfm.label.sub <- dfm.label.sub[dfm.label.sub$pvalue<threshold_p,]
dfm.label.sub <- na.omit(dfm.label.sub)

association_del <- P1  +  
  geom_label_repel(data = dfm.label.sub, aes(label = label, col="black"), max.overlaps = Inf, min.segment.length = 0, label.size = 0, box.padding = 0.5,
                   fill = NA) +  
  geom_point(data = dfm_novel_sig, aes(col="lightskyblue")) +
  geom_point(data = dfm.sub,aes(col="black"),shape = 5)

#duplication 
eas_association <- eas_association_all[eas_association_all$EAS_CNV_type=="dup",]
eas_association.sig <- eas_association[eas_association$EAS_P<threshold_p,]

#x axis index
eas_association <- eas_association[order(eas_association$BP),]
eas_association <- eas_association[mixedorder(eas_association$CHR),]
eas_association$index <- 1:nrow(eas_association)

# caculate max y value
ymax_val <- max(-log10(min(eas_association$EAS_P)))

#prepare label - p-threshold
thresh_lable <- data.frame(
  x = 0,
  y = (-log10(signif(threshold_permute_p)))+1,
  text = paste0(" Permutation threshold = ",signif(threshold_permute_p,digits = 2))
)

P1 <- ggman(eas_association, snp = "locus", bp = "start", chrom = "CHR", pvalue = "EAS_P",  pointSize = 0.6, cex.axis=0.8, ymax = ymax_val+2, sigLine = NA) + 
  theme_classic() + 
  geom_hline(yintercept=-log10(signif(threshold_permute_p, digits = 3)),linetype=2) +
  xlab("Chromosome") + 
  ylab("-log10(p-value)") +
  theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1)) +  
  geom_text(data = thresh_lable, aes(x, y, label = text, col = "black"), hjust = 0) +  # Left align text
  ggtitle("Duplication") + 
  scale_colour_manual(values = c("darkgrey","lightgrey","black","lightskyblue")) 

#prepare label- significant loci 
dfm <- P1[[1]]
dfm_novel <- dfm[dfm$locus %in% c("chr11p15.5", "chr11q13.1", "chr11q13.2", "chr17q25.3", "chr18q23", "chr19p13.3", "chr19q13.42", "chr5q31.3", "chr8p21.3", "chr8q24.3"),]
dfm_novel_sig <- dfm_novel[dfm_novel$EAS_P < threshold_p,]

filtered_data <- eas_association.sig %>% 
  group_by(locus) %>% 
  dplyr::filter(EAS_P== min(EAS_P))
filtered_data$snp <- filtered_data$locus
filtered_data <- filtered_data[!duplicated(filtered_data$locus),]
filtered_data_dup <- filtered_data

labelDfm <- as.data.frame(filtered_data)
labelDfm$snp <- labelDfm[,"EAS_GENE"]
labelDfm$label <- labelDfm[,"EAS_GENE"]
dfm.sub <- dfm[dfm$EAS_GENE %in% labelDfm$EAS_GENE,]
dfm.sub$label <- paste0(dfm.sub$locus)
dfm.label.sub <- dfm.sub[!duplicated(dfm.sub$locus),]
dfm.label.sub <- dfm.label.sub[dfm.label.sub$pvalue<threshold_p,]

association_dup <- P1  +  
  geom_label_repel(data = dfm.label.sub, aes(label = label, col="black"), max.overlaps = Inf, min.segment.length = 0,  label.size = 0, box.padding = 0.5,
                   fill = NA) +  
  geom_point(data = dfm_novel_sig, aes(col="lightskyblue")) +
  geom_point(data = dfm.sub,aes(col="black"),shape = 5)

plot_legend <- ggplot(eas_association, aes(x = EAS_P, y = EAS_BETA)) + 
  geom_point(aes(col="Novel Regions"))+
  geom_point(aes(col="All Regions"))+
  theme(legend.position="bottom",legend.title = element_blank()) +
  scale_color_manual(values = cohort_colors, name="") +
  scale_linetype_manual(values=c("Odds Ratio = 1"="dashed","95% Confidence Interval"="solid"), name="") +
  scale_shape_manual(values=c("p-value <  0.05" = 16,"p-value >  0.05" = 21), name="") + theme_classic() 

legend <- get_legend(
  # create some space to the left of the legend
  plot_legend+
    theme(legend.position="bottom"))

aligned_plots <- align_plots(association_del, association_dup, align = 'v')
plot_association <- plot_grid(plotlist = aligned_plots,ncol=1)


loci_result <-rbind(filtered_data_del,filtered_data_dup)
loci_result <- loci_result[,c("seqname", "start", "end", "locus", "CNV_type", "EAS_P", "EAS_BETA")]
loci_result$seqname <- gsub("chr","",loci_result$seqname)
loci_result <- write.table(loci_result,"output/loci_result.txt",row.names=F,col.names=T,quote=F,sep='\t')

#freq-OR comparison
eas_association_del_regenie <- read.table("input/del.step2_AFF.regenie",header=T)
eas_association_del_regenie$CNV_type <- "Deletion"
eas_association_dup_regenie <- read.table("input/dup.step2_AFF.regenie",header=T)
eas_association_dup_regenie$CNV_type <- "Duplication"

eas_association_all_regenie <- rbind(eas_association_del_regenie,eas_association_dup_regenie)
eas_association_all_res <- eas_association_all_res[eas_association_all_res$EAS_CNV_type!="all",]
eas_association_all <- merge(eas_association_all_regenie,eas_association_all_res,by.y="EAS_GENE",by.x="ID")

threshold_p <- threshold_p
EAS_sig <- eas_association_all[eas_association_all$LOG10P>(-log10(threshold_p)),]
eas_association_all$CNV_type <- factor(eas_association_all$CNV_type,level=c("all","del","dup"))
EAS_sig$CNV_type <- factor(EAS_sig$CNV_type,level=c("all CNVs","Deletion","Duplication"))
EAS_sig$locus <- factor(EAS_sig$locus,level=unique(EAS_sig$locus))
EAS_sig_novel <- EAS_sig[EAS_sig$locus %in% c("chr19p13.3", "chr8p21.3","chr11q13.1"),]

EAS_sig_filtered_del <- EAS_sig[EAS_sig$CNV_type=="Deletion",] %>% 
  group_by(locus) %>% 
  dplyr::filter(BETA == max(BETA))
EAS_sig_filtered_del <- EAS_sig_filtered_del[!duplicated(EAS_sig_filtered_del$locus),]

EAS_sig_filtered_dup <- EAS_sig[EAS_sig$CNV_type=="Duplication",] %>% 
  group_by(locus) %>% 
  dplyr::filter(BETA == max(BETA))
EAS_sig_filtered_dup <- EAS_sig_filtered_dup[!duplicated(EAS_sig_filtered_dup$locus),]

EAS_sig_filtered <- rbind(EAS_sig_filtered_del,EAS_sig_filtered_dup)

OR_Freq_del <- ggplot(EAS_sig[EAS_sig$CNV_type=="Deletion",],aes(x=A1FREQ,y=2^BETA,col=CNV_type)) + 
  geom_point() +
  geom_hline(yintercept=1,linetype=2) +
  theme_classic() +
  theme(plot.margin = unit(c(0.9, 0.1, 0.1, 0.1), "cm")) +
  scale_color_d3("category20c") +
  theme(legend.position="none") +
  xlab("rCNV Freqency") +
  ylab("Odds ratio") + facet_grid(~CNV_type)

OR_Freq_dup <- ggplot(EAS_sig[EAS_sig$CNV_type=="Duplication",],aes(x=A1FREQ,y=2^BETA,col=CNV_type)) + 
  geom_point() +
  geom_hline(yintercept=1,linetype=2) +
  theme_classic() +
  theme(plot.margin = unit(c(0.8, 0.2, 0.2, 0.2), "cm")) +
  scale_color_d3("category20c") +
  theme(legend.position="none") +
  xlab("CNV Freqency") +
  ylab("Odds ratio") + facet_grid(~CNV_type)

OR_Freq <- ggplot(EAS_sig,aes(x=A1FREQ,y=2^BETA,col=CNV_type)) + 
  geom_point() +
  geom_label_repel(data = EAS_sig_filtered[!is.na(EAS_sig_filtered$locus),], aes(x=A1FREQ, y=2^BETA,label = paste0(locus)), max.overlaps = Inf,inherit.aes = FALSE) + 
  geom_hline(yintercept=1,linetype=2) +
  theme_classic() +
  scale_x_log10(breaks = c(0.0005, 0.001, 0.002, 0.003),
                labels = c("0.0005", "0.001", "0.002", "0.003")) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.4), "cm")) +
  scale_color_manual(values = cohort_colors) +
  theme(legend.position="none",legend.title = element_blank()) +
  xlab("CNV Freqency") +
  ylab("Odds ratio") 

OR_Freq_ld <- ggplot(EAS_sig,aes(x=A1FREQ,y=2^BETA,col=CNV_type)) + 
  geom_label_repel(data = EAS_sig_filtered[!is.na(EAS_sig_filtered$locus),], aes(x=A1FREQ, y=2^BETA,label = paste0(locus,"-",gene_name)), max.overlaps = Inf,inherit.aes = FALSE) + 
  geom_point() +theme_classic() +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.4), "cm")) +
  scale_x_log10() + 
  scale_color_manual(values = cohort_colors) +
  theme(legend.position="bottom",legend.title = element_blank()) 
  
OR_Freq_legend <- get_legend(
  # create some space to the left of the legend
  OR_Freq_ld+
    theme(legend.position="bottom"))

OR_Freq_final <- plot_grid( OR_Freq, OR_Freq_legend, nrow = 2, rel_heights = c(1,0.1))

#EAS-EUR frequency
freq <- read.table("input/EAS_EUR_freq.txt",header=T)
freq_sig <- freq[freq$Gene %in% EAS_sig$ID,]
freq_sig$EAS_freq <- (freq_sig$EAS_Case + freq_sig$EAS_Control)/44161
freq_sig$EUR_freq <- (freq_sig$EUR_Case + freq_sig$EUR_Control)/34257
freq_sig_novel <- freq_sig[freq_sig$locus %in% c("chr11p15.5","chr11q13.1","chr11q13.2","chr17q25.3","chr18q23","chr19p13.3","chr19q13.42","chr5q31.3","chr8p21.3","chr8q24.3"),]
to_plot <- melt(freq_sig_novel[,c("Gene","EAS_freq","EUR_freq")])
to_plot$Population <- gsub("_freq","",to_plot$variable)
freq_diff <- wilcox.test(to_plot[to_plot$Population=="EAS","value"],to_plot[to_plot$Population=="EUR","value"])
Freq_com <- ggplot(to_plot,aes(x=variable,y=value,col=Population, group=Gene)) + 
  geom_jitter(width = 0.2) +
  theme_classic() +
  scale_color_manual(values = cohort_colors) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.4), "cm")) +
  theme(legend.position="none",legend.title = element_blank()) +
  xlab("Datasets")  +
  ylab("Frequency") +
  geom_text(aes(x=1.5, y=max(value) * 1.05, label=paste0("p-value = ",signif(freq_diff$p.value,2)), col=NULL), show.legend = FALSE)

Freq_com_ld <- ggplot(to_plot,aes(x=variable,y=value,col=Population, group=Gene)) + 
  geom_jitter(width = 0.2) +
  theme_classic() +
  scale_color_manual(values = cohort_colors) +
  theme(legend.position="none",legend.title = element_blank()) 
  
Freq_com_legend <- get_legend(
  # create some space to the left of the legend
  Freq_com_ld+
    theme(legend.position="bottom"))

Freq_com_final <- plot_grid( Freq_com, Freq_com_legend, nrow = 2, rel_heights = c(1,0.1))


plot_feature <- plot_grid( OR_Freq_final, Freq_com_final, labels = c("b","c"), nrow = 1)
plot_final <- plot_grid(plot_association,legend, plot_feature, labels = c("a",""), nrow = 3,rel_heights = c(2,0.2,1))


ggsave("output/Figure_3_EAS_CNV_gene_association_same_thresh_final_v3.pdf", plot_final, width=9, height=9, units = "in")
