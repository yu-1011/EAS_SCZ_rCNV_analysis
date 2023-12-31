setwd("/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/02extended_figure/extended_fig2_CNV_QC")
rm(list=ls())
library("ggridges")
library("ggsci")
library("ggplot2")
library("gridExtra")
library("cowplot")
library("reshape2")

quality_ori <- read.table("cfile_final.datMeta",header=T)
array_file <- read.delim("/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/misc/cohort_details.txt")
array_file <- array_file[array_file$dataset_id!="six1",]
array_file <- array_file[array_file$dataset_id!="umc1",]
quality_ori <- merge(quality_ori,array_file,by.x="Cohort",by.y="cohort")
quality_ori$CNV_average_length <- quality_ori$CNV_total_length / quality_ori$CNV_num

array_file <- array_file[order(array_file$array_type),]
quality_ori$dataset_id <- factor(quality_ori$dataset_id, levels=unique(array_file$dataset_id))
quality <- quality_ori[quality_ori$NSEG!=0,]

p_LRR_mean <- ggplot(quality_ori, aes(x = LRR_mean)) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  labs(title="", x ="Average Log R Ratio", y = "datasets") +
  theme_classic() + 
  theme(legend.position="none") +
  scale_fill_d3("category20c") 

p_LRR_SD <- ggplot(quality_ori, aes(x = LRR_SD)) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  labs(title="", x="Standard error of Log R Ratio", y = "datasets") +
  theme_classic() + 
  theme(legend.position="none") +
  scale_fill_d3("category20c") 

p_BAF_mean <- ggplot(quality_ori, aes(x = BAF_mean)) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  labs(title="", x ="Average B Allele Frequency", y = "datasets") +
  theme_classic() + 
  theme(legend.position="none") +
  scale_fill_d3("category20c") 

p_BAF_SD <- ggplot(quality_ori, aes(x = BAF_SD)) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  labs(title="", x ="Standard error of B Allele Frequency", y = "datasets") +
  theme_classic() + 
  theme(legend.position="none") +
  scale_fill_d3("category20c") 

p_WF <- ggplot(quality_ori, aes(x = WF)) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  labs(title="", x ="Waviness Factor", y = "datasets") +
  theme_classic() + 
  theme(legend.position="none") +
  scale_fill_d3("category20c") 

p_GCWF <- ggplot(quality_ori, aes(x = GCWF)) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  labs(title="", x ="WF adjusted by gcmodel", y = "datasets") +
  theme_classic() + 
  scale_fill_d3("category20c") 

p_CNV_num <- ggplot(quality, aes(x = as.numeric(NSEG))) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  labs(title="CNV num per carrier \n- remove non-CNV carriers", x ="CNV number", y = "datasets") +
  theme_classic() + 
  theme(legend.position="none") +
  scale_fill_d3("category20c") 

p_CNV_num_all <- ggplot(quality_ori, aes(x = NSEG)) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  labs(title="CNV num per carrier \n- all samples", x = "CNV number", y = "datasets") +
  theme_classic() + 
  theme(legend.position="none") +
  scale_fill_d3("category20c") 

p_CNV_average_length <- ggplot(quality, aes(x = log2(CNV_average_length))) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  labs(title="Average CNV length \n- remove non-CNV carriers", x ="average CNV length", y = "datasets") +
  theme_classic() + labs(fill = "Array") +
  scale_fill_d3("category20c") 

p_CNV_average_length_all <- ggplot(quality, aes(x = log2(CNV_average_length))) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  labs(title="Average CNV length \n- all samples", x ="average CNV length", y = "datasets") +
  theme_classic() + 
  theme(legend.position="none") +
  scale_fill_d3("category20c") 

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  p_CNV_average_length)

p_quality_all <- plot_grid(
  p_LRR_mean, p_BAF_mean, p_WF, p_LRR_SD, p_BAF_SD , p_GCWF+theme(legend.position="none"), 
  labels = "auto", nrow = 2
)
p_quality_final <- plot_grid(p_quality_all,legend, rel_widths = c(4, 1))
ggsave("output/EF_02_Ridges_EAS_CNV_QCed_quality_results.pdf" , p_quality_final, width=12, height=6)
