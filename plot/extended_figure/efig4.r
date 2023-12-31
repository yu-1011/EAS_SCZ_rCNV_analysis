setwd("/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/02extended_figure/extended_fig2_CNV_QC")
rm(list=ls())
library("ggsci")
library("ggplot2")
library("gridExtra")
library("cowplot")
library("reshape2")

#data prep
quality_eur <- read.table("PGC_41K_QC.cnv.indiv",header=T)
quality_eur$pop_lable <- "European"
quality_eur$Population <- "EUR"

quality_eas <- read.table("cfile_final.datMeta",header=T)
quality_eas$pop_lable <- "East Asian"
quality_eas$Population <- "EAS"

quality_eas <- quality_eas[c("FID","IID","PHE","NSEG","KB","KBAVG","pop_lable","Population")]

quality_ori <- rbind(quality_eas,quality_eur)
quality <- quality_ori[quality_ori$NSEG!=0,]


#color setting
cohort_colors <- c(
  EUR = "#BE1E2D",
  EAS = "#0087FF")

p_NSEG <- ggplot(quality, aes(y = NSEG)) +
  geom_boxplot(aes(x=pop_lable,fill=Population)) +
  labs(title="CNV num per carrier \n-  remove non-CNV carriers", y="CNV number", x = "datasets") +
  theme_classic() + 
  ylim(c(0,30)) +
  theme(legend.position="none") +
  scale_fill_manual(values = cohort_colors) + 
  coord_flip() 

p_NSEG_all <- ggplot(quality_ori, aes(y = NSEG)) +
  geom_boxplot(aes(x=pop_lable,fill=Population)) +
  labs(title="CNV num per carrier \n- all samples", y="CNV number", x = "datasets") +
  theme_classic() + 
  theme(legend.position="none") +
  ylim(c(1,30)) +
  scale_fill_manual(values = cohort_colors) + 
  coord_flip() 

p_KBAVG <- ggplot(quality, aes(y = log2(KBAVG))) +
  geom_boxplot(aes(x=pop_lable,fill=Population)) +
  labs(title="Average CNV length \n- remove non-CNV carriers", y="average CNV length", x = "datasets") +
  theme_classic() + 
  scale_fill_manual(values = cohort_colors) + 
  coord_flip() 

p_KBAVG_all <- ggplot(quality_ori, aes(y = log2(KBAVG))) +
  geom_boxplot(aes(x=pop_lable,fill=Population)) +
  labs(title="Average CNV length \n- all samples", y="average CNV length", x = "datasets") +
  theme_classic() + 
  theme(legend.position="none") +
  scale_fill_manual(values = cohort_colors) + 
  coord_flip() 

legend <- get_legend(
  # create some space to the left of the legend
  p_KBAVG)

p_cnv_all <- plot_grid(
  p_NSEG_all, p_KBAVG_all+theme(legend.position="none"), p_NSEG, p_KBAVG+theme(legend.position="none"), 
  labels = "AUTO", nrow = 2
)
p_cnv_final <- plot_grid(p_cnv_all,legend, rel_widths = c(4, 1))
ggsave("output/EF_04_boxplot_PGC_comparison_log2_CNV_QCed_results.pdf" , p_cnv_final, width=12, height=6)
