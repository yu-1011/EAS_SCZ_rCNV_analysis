setwd("/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/02extended_figure/extended_fig2_CNV_QC")
rm(list=ls())
library("ggridges")
library("ggsci")
library("ggplot2")
library("gridExtra")
library("cowplot")
library("reshape2")
library("scales")

quality_ori <- read.table("cfile_final.datMeta",header=T)
quality_before <- read.table("02_4_rawcnv_quality_summary_all.txt",header=T)
array_file <- read.delim("/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/misc/cohort_details.txt")
quality_ori <- merge(quality_ori,array_file,by.x="Cohort",by.y="cohort")
quality_ori <- quality_ori[!(quality_ori$dataset_id%in%c("six1","umc1")),]
quality_ori$CNV_average_length <- quality_ori$CNV_total_length / quality_ori$CNV_num
quality_before <- merge(quality_before,array_file,by.x="Cohort",by.y="cohort")
quality_before <- quality_before[!(quality_before$dataset_id%in%c("six1","umc1")),]
quality_before$CNV_average_length <- quality_before$CNV_total_length / quality_before$CNV_num


array_file <- array_file[order(array_file$array_type),]
quality_ori$dataset_id <- factor(quality_ori$dataset_id, levels=unique(array_file$dataset_id))
quality_before$dataset_id <- factor(quality_before$dataset_id, levels=unique(array_file$dataset_id))
quality <- quality_ori[quality_ori$NSEG!=0,]

quality_ori$type <- "after"
quality_before$type <- "before"

array_file <- array_file[order(array_file$array_type),]
quality <- quality_ori[quality_ori$NSEG!=0,]

p_CNV_num <- ggplot(quality, aes(x = as.numeric(CNV_num))) +
  scale_x_continuous(trans = log2_trans(),
                       breaks = trans_breaks("log2", function(x) 2^x),
                       labels = trans_format("log2", math_format(2^.x))) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  labs(title="CNV num per carrier \n- remove non-CNV carriers", x ="CNV number", y = "datasets") +
  theme_classic() + 
  theme(legend.position="none") +
  scale_fill_d3("category20c") 

p_CNV_num_all <- ggplot(quality_ori, aes(x = log2(as.numeric(CNV_num)+1))) +
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  labs(title="CNV num per carrier \n- all samples", x = "CNV number", y = "datasets") +
  theme_classic() + 
  theme(legend.position="none") +
  scale_fill_d3("category20c") 

p_CNV_average_length <- ggplot(quality, aes(x = log2(KBAVG))) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(title="Average CNV length \n- remove non-CNV carriers", x ="average CNV length", y = "datasets", fill="Array Type") +
  theme_classic() + 
  scale_fill_d3("category20c") 

p_CNV_average_length_all <- ggplot(quality, aes(x = log2(KBAVG))) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(title="Average CNV length \n- all samples", x ="average CNV length", y = "datasets") +
  theme_classic() + 
  theme(legend.position="none") +
  scale_fill_d3("category20c") 

model <- lm(KBAVG ~ array_type, data = quality)
quality$CNV_KBAVG_adjusted <- predict(model, newdata = quality)
p_KBAVG_adjusted <- ggplot(quality, aes(x = log2(CNV_KBAVG_adjusted))) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(title="After adjustment", x ="average CNV length", y = "datasets", fill = "Array Type") +
  theme_classic() + 
  scale_fill_d3("category20c") 

p_KBAVG <- ggplot(quality, aes(x = log2(KBAVG))) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(title="Before adjustment", x ="average CNV length", y = "datasets") +
  theme_classic() + 
  scale_fill_d3("category20c") 

model <- lm(KB ~ array_type, data = quality)
quality$CNV_KB_adjusted <- predict(model, newdata = quality)
p_KB_adjusted <- ggplot(quality, aes(x = log2(CNV_KB_adjusted))) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(title="After adjustment", x ="total CNV length", y = "datasets") +
  theme_classic() + 
  scale_fill_d3("category20c") 

p_KB <- ggplot(quality, aes(x = log2(KB))) +
  geom_density_ridges(aes(y=dataset_id,fill=array_type)) +
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(title="Before adjustment", x ="total CNV length", y = "datasets") +
  theme_classic() + 
  scale_fill_d3("category20c") 

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  p_CNV_average_length)

p_cnv_all <- plot_grid(
  p_KBAVG+theme(legend.position="none"), p_KBAVG_adjusted+theme(legend.position="none"), p_KB+theme(legend.position="none"), p_KB_adjusted+theme(legend.position="none"), 
  labels = "auto", nrow = 2
)
p_cnv_final <- plot_grid(p_cnv_all,legend, rel_widths = c(4, 1))
ggsave(filename = "output/EF_03_Ridges_EAS_log2_CNV_QCed_results.pdf" , plot=p_cnv_final, width=12, height=6)
