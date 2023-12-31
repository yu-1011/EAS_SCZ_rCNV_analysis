library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggstance)
library(ggforestplot)
library(ggthemes)
library(ggsci)
library(cowplot)
library(metafor)

rm(list=ls())

wd="/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/02extended_figure/extended_fig3_CNV_burden/"
setwd(wd)
file_name <- list.files("input/",pattern="*study*")
burden_study <- read.table(paste0("input/",file_name),header=T)
array_file <- read.delim("/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/misc/cohort_annotation.txt")

# study plot -all CNV
psignif <- 0.05 
burden_study <- burden_study[burden_study$CNV_type=="allCNV",]
burden_study$data_set <- factor(burden_study$data_set,level=unique(burden_study$data_set))
burden_study <- merge(burden_study,array_file,by.x="data_set",by.y="cohort",all.x = T)
burden_study[burden_study$data_set=="combined","array_type"] <- "combined"
burden_study[burden_study$data_set=="combined","dataset_id"] <- "combined"
burden_study$array_type <- factor(burden_study$array_type,level=c("ASA","GSAA","Il1M","Om25","PSYC","combined"))
burden_study <-burden_study[order(burden_study$array_type),]
burden_study$dataset_id <- factor(burden_study$dataset_id ,level=unique(burden_study$dataset_id))
burden_study$CNV_size <- factor(burden_study$CNV_size,level=c(">20kb",">100kb",">200kb",">300kb",">400kb",">500kb",">600kb"))
burden_study <- burden_study[burden_study$CNV_type=="allCNV"&burden_study$CNV_size!=">20kb",]
burden_study <-
  burden_study %>%
  dplyr::mutate(.filled_gene = NGENE_glm_pval < psignif)

plot_study_gene <- ggplot(data = burden_study, aes(x = NGENE_glm_OR, y = dataset_id)) +
  geom_effect(ggplot2::aes(xmin = NGENE_glm_lowerCI,xmax = NGENE_glm_upperCI,color = array_type,filled =.filled_gene),position = ggstance::position_dodgev(height = 0.5)) +
  theme_classic() +
  geom_stripes(odd = "#33333333", even = "#00000000") +
  geom_vline(xintercept = 1, linetype=2) +
  xlab("Odds ratio") +
  #ggtitle("Number of gene affected by CNVs") +
  ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, size=8, hjust=1)) +
  theme(legend.position="right") +
  theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 10, face = "plain")) +
  theme(plot.tag = element_text(face='bold')) +
  labs(color = "Array") +
  scale_color_d3("category20c") +
  facet_wrap(~CNV_size)


ggsave("output/EF_07_CNV_burden_study_allCNV.pdf", plot_study_gene, width=12, height=8 ,units = "in")
