library(ggplot2)
library(ggman)
library(tidyverse)
library(ggsci)
library(cowplot)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(gtools)

rm(list=ls())
wd="/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/01main_figure/figure4_meta_analysis"
setwd(wd)

cohort_colors <- c(
  EUR = "#BE1E2D",
  EAS = "#3182BD",
  'Deletion' = "#ff7f0e",
  'Duplication' = "#2ca02c",
  'Novel Regions' = "#7ec0ee",
  'All Regions' = "#A9A9A9")

## only run once
# meta_association_file_name <- list.files("input/",pattern="*gene_assoc_meta*.tsv$")
# meta_association_all <- read.table(paste0("input/",meta_association_file_name),header=T)
# gene_annot <- read.csv("/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/misc/gencode.v19.gene.annot.csv")
# cytoband <- read.table("/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/misc/cytoBand.hg19.txt.gz")
# names(cytoband) <- c("chr","start","end","pos","type")
# cytoband$name <- paste0(cytoband$chr,cytoband$pos)
# meta_association_all <- read.table("input/META_PGC_EAS.all_1.txt",header=T)
# meta_association_all$META_CNV_type <-"all"
# meta_association_del <- read.table("input/META_PGC_EAS.del_1.txt",header=T)
# meta_association_del$META_CNV_type <- "del"
# meta_association_dup <- read.table("input/META_PGC_EAS.dup_1.txt",header=T)
# meta_association_dup$META_CNV_type <- "dup"
# 
# meta_association_all <- rbind(meta_association_del,meta_association_dup)
# names(meta_association_all) <- paste0("META_",names(meta_association_all))
# meta_association_all <- merge(meta_association_all,gene_annot,by.x="META_MarkerName",by.y="gene_id")
# meta_association_all$locus <- NA
# 
# # Iterate over each gene and check if it is within any cytoband range
# for (i in 1:nrow(meta_association_all)) {
#   gene_chromosome <- meta_association_all$seqname[i]
#   gene_start <- meta_association_all$start[i]
#   gene_end <- meta_association_all$end[i]
# 
#   for (j in 1:nrow(cytoband)) {
#     cytoband_chromosome <- cytoband$chr[j]
#     cytoband_start <- cytoband$start[j]
#     cytoband_end <- cytoband$end[j]
# 
#     if (gene_chromosome == cytoband_chromosome && gene_start >= cytoband_start && gene_end <= cytoband_end) {
#       #print(paste("Gene", i, "is located within cytoband", cytoband$name[j]))
#       meta_association_all[i,"locus"] <- cytoband$name[j]
#       break  #
#     }
#   }
# }
# res <- meta_association_all[,c("META_MarkerName","locus","META_P.value","META_Zscore","META_META_CNV_type","start","end","seqname","gene_name")]
# 
# names(res) <- c("META_GENE","locus","META_P","META_Zscore","META_CNV_type","start","end","seqname","gene_name")
# write.table(res, "input/gene_assoc_META.txt", col.names = T, row.names = F, quote = F, sep = '\t' )
# write.table(res, "output/gene_assoc_META.txt", col.names = T, row.names = F, quote = F, sep = '\t' )

meta_association_all <- read.table("input/gene_assoc_META.txt", header=T)
meta_association_all[meta_association_all$META_GENE=="ENSG00000005075",]$locus <- NA
meta_association_all$BP <- meta_association_all$start + ((meta_association_all$end - meta_association_all$start)/2)
meta_association_all$length <- meta_association_all$end - meta_association_all$start
meta_association_all$CHR <- gsub("chr","",meta_association_all$seqname)
meta_association_all$CHR <- ifelse(meta_association_all$CHR=="X",23,as.numeric(meta_association_all$CHR))
meta_association_all$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = meta_association_all$META_GENE, column = "SYMBOL", keytype = "ENSEMBL")

#threshold_exome_p <-  (0.05/(18531*2))
threshold_exome_p <-   2.7202000000000005e-05*(1667/1678)
threshold_permute_p <- 2.7202000000000005e-05*(1667/1678)

# PLOT deletion
# data.frame prepare
meta_association <- meta_association_all[meta_association_all$META_CNV_type=="del",]
meta_association$META_CNV_type <- gsub("del", "Deletion", meta_association$META_CNV_type)
#x axis index
meta_association <- meta_association[order(meta_association$BP),]
meta_association <- meta_association[mixedorder(meta_association$CHR),]
meta_association$index <- 1:nrow(meta_association)
# caculate max y value
ymax_val <- max(-log10(min(meta_association$META_P)))

#prepare label - p-threshold
threshold_FWER_p <- threshold_exome_p 
#threshold_BH_p <- threshold_del_permute_p 
threshold_p <- threshold_exome_p

thresh_lable <- data.frame(
  x = 0,
  y = (-log10(signif(threshold_p)))+2,
  text = paste0(" Permutation threshold = ",signif(threshold_BH_p,digits = 2))
)

meta_association.sig <- meta_association[meta_association$META_P<threshold_p,]
filtered_data <- meta_association.sig %>% 
  group_by(locus) %>% 
  dplyr::filter(META_P == min(META_P))
filtered_data$snp <- filtered_data$locus

P1 <- ggman(meta_association, snp = "locus", bp = "start", chrom = "CHR", pvalue = "META_P",  pointSize = 0.6, cex.axis=0.8, ymax = ymax_val+2, sigLine = NA) + 
  theme_classic() + 
  geom_hline(yintercept=-log10(signif(threshold_BH_p, digits = 3)),linetype=2) +
  xlab("Chromosome") + 
  ylab("-log10(p-value)") +
  theme(axis.text.x = element_text(angle = 45, size=8, hjust=1)) +  
  geom_text(data = thresh_lable, aes(x, y, label = text, col="black"), hjust = 0) +
  ggtitle("Deletion") + 
  scale_colour_manual(values = c("darkgrey","lightgrey","black","lightskyblue"))  +
  scale_fill_manual(values = c("lightskyblue2","midnightblue","black","IndianRed1")) 

#prepare label- significant loci 
dfm <- P1[[1]]
dfm_novel <- dfm[dfm$locus %in% c("chr1q21.2","chr11p15.5", "chr11q13.1", "chr11q13.2", "chr17q25.3", "chr18q23", "chr19p13.3", "chr19q13.42", "chr5q31.3", "chr7q22.1", "chr8p21.3", "chr8q24.3"),]
dfm_novel_sig <- dfm_novel[dfm_novel$META_P < threshold_p,]
labelDfm <- as.data.frame(filtered_data)
labelDfm$snp <- labelDfm[,"META_GENE"]
labelDfm$label <- labelDfm[,"META_GENE"]
dfm.sub <- dfm[dfm$META_GENE %in% labelDfm$META_GENE,]
#dfm.sub$label <- paste0(dfm.sub$locus,"_",dfm.sub$Gene_Symbol)
dfm.sub$label <- paste0(dfm.sub$locus)
dfm.sub$label <- gsub("NA","",dfm.sub$label)
dfm.label.sub <- dfm.sub[!duplicated(dfm.sub$locus),]
print(paste0("del sig loci number:", dim(dfm.sub)[1]))


association_del <- P1  +  
  geom_label_repel(data = dfm.label.sub, aes(label = label, col="black"), label.size = 0, max.overlaps = Inf,box.padding = 0.5,
                   fill = NA) +  
  geom_point(data = dfm_novel_sig, aes(col="lightskyblue")) +
  geom_point(data = dfm.sub,aes(col="black"),shape = 5)

#duplication 
meta_association <- meta_association_all[meta_association_all$META_CNV_type=="dup",]
meta_association$META_CNV_type <- gsub("dup", "Duplication", meta_association$META_CNV_type)

#x axis index
meta_association <- meta_association[order(meta_association$BP),]
meta_association <- meta_association[mixedorder(meta_association$CHR),]
meta_association$index <- 1:nrow(meta_association)

# caculate max y value
ymax_val <- max(-log10(min(meta_association$META_P)))

#prepare label - p-threshold
threshold_FWER_p <- threshold_exome_p 
#threshold_BH_p <- threshold_dup_permute_p 
threshold_p <- threshold_exome_p

thresh_lable <- data.frame(
  x = 0,
  y = (-log10(signif(threshold_p)))+1,
  text = paste0(" Permutation threshold = ",signif(threshold_BH_p,digits = 2))
)

meta_association.sig <- meta_association[meta_association$META_P<threshold_p,]
filtered_data <- meta_association.sig %>% 
  group_by(locus) %>% 
  dplyr::filter(META_P == min(META_P))
filtered_data$snp <- filtered_data$locus

P1 <- ggman(meta_association, snp = "locus", bp = "start", chrom = "CHR", pvalue = "META_P",  pointSize = 0.6, cex.axis=0.8, ymax = ymax_val+2, sigLine = NA) + 
  theme_classic() + 
  #geom_hline(yintercept=-log10(signif(threshold_FWER_p, digits = 3)),linetype=2) +
  geom_hline(yintercept=-log10(signif(threshold_BH_p, digits = 3)),linetype=2) +
  xlab("Chromosome") + 
  ylab("-log10(p-value)") +
  theme(axis.text.x = element_text(angle = 45, size=8, hjust=1)) +  
  geom_text(data = thresh_lable, aes(x, y, label = text, col="black"), hjust = 0) +
  #geom_text(data = thresh_BH_lable, aes(x, y, label = text, col="black")) +
  ggtitle("Duplication") + 
  scale_colour_manual(values = c("darkgrey","lightgrey","black","lightskyblue"))  +
  scale_fill_manual(values = c("lightskyblue2","midnightblue","black","IndianRed1")) 

#prepare label- significant loci 
dfm <- P1[[1]]
dfm_novel <- dfm[dfm$locus %in% c("chr11p15.5", "chr11q13.1", "chr11q13.2", "chr17q25.3", "chr18q23", "chr19p13.3", "chr19q13.42", "chr5q31.3", "chr7q22.1", "chr8p21.3", "chr8q24.3"),]
dfm_novel_sig <- dfm_novel[dfm_novel$META_P < threshold_p,]
labelDfm <- as.data.frame(filtered_data)
labelDfm$snp <- labelDfm[,"META_GENE"]
labelDfm$label <- labelDfm[,"META_GENE"]
dfm.sub <- dfm[dfm$META_GENE %in% labelDfm$META_GENE,]
#dfm.sub$label <- paste0(dfm.sub$locus,"_",dfm.sub$Gene_Symbol)
dfm.sub$label <- paste0(dfm.sub$locus)
dfm.sub$label <- gsub("NA","",dfm.sub$label)
dfm.label.sub <- dfm.sub[!duplicated(dfm.sub$locus),]
print(paste0("del sig loci number:", dim(dfm.sub)[1]))


association_dup <- P1  +  
  geom_label_repel(data = dfm.label.sub, aes(label = label, col="black"), label.size = 0, max.overlaps = Inf,box.padding = 0.5,
                   fill = NA) +  
  geom_point(data = dfm_novel_sig, aes(col="lightskyblue")) +
  geom_point(data = dfm.sub,aes(col="black"),shape = 5)

plot_legend <- ggplot(meta_association, aes(x = META_P, y = META_Zscore)) + 
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


#plot_association <- plot_grid(association_all, association_del, association_dup, nrow = 3)
plot_association <- plot_grid( association_del, association_dup, legend, nrow = 3, labels = c("","",""), rel_heights = c(3,3,0.2))

freq <- read.table("input/EAS_EUR_freq.txt",header=T)
freq[freq$Gene=="ENSG00000005075",]$locus <- NA

threshold_p <- threshold_exome_p
meta_association.sig <- meta_association_all[meta_association_all$META_P<threshold_p&meta_association_all$META_CNV_type != "all",]

freq_sig <- freq[paste0(freq$Gene,freq$CNV_type) %in% paste0(meta_association.sig$META_GENE,meta_association.sig$META_CNV_type),]
freq_sig$EAS_freq <- (freq_sig$EAS_Case + freq_sig$EAS_Control)/44161
freq_sig$EUR_freq <- (freq_sig$EUR_Case + freq_sig$EUR_Control)/34257
#freq_sig_novel <- freq_sig[freq_sig$locus %in% c("chr11p15.5","chr11q13.1","chr11q13.2","chr17q25.3","chr18q23","chr19p13.3","chr19q13.42","chr5q31.3","chr7q22.1","chr8p21.3","chr8q24.3"),]

freq_sig$EAS_OR <-ifelse(freq_sig$EAS_BETA=="-",0,2^as.numeric(freq_sig$EAS_BETA))
freq_sig$EUR_OR <-ifelse(freq_sig$EUR_BETA!="-",2^as.numeric(freq_sig$EUR_BETA),0)

freq_sig$EAS_VE <- (2*freq_sig$EAS_freq*(1-freq_sig$EAS_freq)) * log(freq_sig$EAS_OR+1)^2
freq_sig$EUR_VE <- (2*freq_sig$EUR_freq*(1-freq_sig$EUR_freq)) * log(freq_sig$EUR_OR+1)^2

freq_sig_del <- freq_sig[freq_sig$CNV_type=="del",c("locus","Gene","EAS_VE","EUR_VE")]
VE_plot <- melt(freq_sig_del)
# sorted VE based on locus
VE_plot <- VE_plot %>% arrange(locus, Gene)
VE_plot$Gene <- as.factor(VE_plot$Gene)
VE_plot$pop <- gsub("_VE","",VE_plot$variable)
VE_del <- ggplot(na.omit(VE_plot),aes(x = Gene, y = value)) +
  geom_bar( aes(fill = pop),stat = "identity",position="fill", colour = NA) +
  theme_classic() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),strip.text = element_text(size = 10)) +
  labs(x = "", y = "% of Variance Explained", title = "Deletion") + 
  scale_fill_manual(values = cohort_colors, name="") +
  theme(legend.position="none") +
  facet_grid(~locus, scales = "free_x")

freq_sig_dup <- freq_sig[freq_sig$CNV_type=="dup",c("locus","Gene","EAS_VE","EUR_VE")]
VE_plot <- melt(freq_sig_dup)
# sorted VE based on locus
VE_plot <- VE_plot %>% arrange(locus, Gene)
VE_plot$Gene <- as.factor(VE_plot$Gene)
VE_plot$pop <- gsub("_VE","",VE_plot$variable)
VE_dup <- ggplot(na.omit(VE_plot),aes(x = Gene, y = value)) +
  geom_bar( aes(fill = pop),stat = "identity",position="fill", colour = NA) +
  theme_classic() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "", y = "% of Variance Explained", title = "Duplication") + 
  scale_fill_manual(values = cohort_colors, name="") +
  theme(legend.position="none") +
  facet_grid(~locus, scales = "free_x")

VE_legend <- ggplot(na.omit(VE_plot),aes(x = Gene, y = value)) +
  geom_bar( aes(fill = pop),stat = "identity",position="fill", colour = NA) +
  scale_fill_manual(values = cohort_colors, name="",
  labels = c(expression(VE[EAS]/(VE[EAS] + VE[EUR])), expression(VE[EUR]/(VE[EAS] + VE[EUR])))
  ) 

legend <- get_legend(
  # create some space to the left of the legend
  VE_legend+
    theme(legend.position="bottom"))

plot_VE <- plot_grid(VE_del,VE_dup,legend,ncol=1,rel_heights = c(3,3,0.3))
plot_final <- plot_grid(plot_association,plot_VE,ncol = 1,labels = 'auto')
ggsave("output/Figure4_META_CNV_gene_association_new.pdf", plot_final, width=9, height=10, units = "in")

