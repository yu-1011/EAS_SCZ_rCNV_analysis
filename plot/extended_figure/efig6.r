library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggstance)
library(ggforestplot)
library(ggthemes)
library(ggsci)
library(cowplot)
library(ggrepel)
library(ggalt)

rm(list=ls())

wd="/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/02extended_figure/extended_fig3_CNV_burden/"
setwd(wd)
file_name <- list.files("input/",pattern="*platform*")
burden_all <- read.table(paste0("input/",file_name),header=T)
array_file <- read.delim("/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/misc/cohort_annotation.txt")

array_file <- array_file[!(array_file$fam_file_cohort %in% c("six1","umc1","xyscz")),]
array_sample <- array_file %>%
  group_by(array_type) %>%
  summarise(array_sample_size = sum(sample_size))
array_sample <- array_sample %>%
  bind_rows(summarise(., across(where(is.numeric), sum),
                      across(where(is.character), ~'combined')))
array_sample$array_sample <- paste0(array_sample$array_type,"\n(n = ",array_sample$array_sample_size," )")
burden_all <- merge(burden_all,array_sample, by.x="data_set", by.y="array_type")

# combined plot - data ini
psignif <- 0.05
burden_all$CNV_type <- factor(burden_all$CNV_type,level=unique(burden_all$CNV_type))
burden_all$CNV_type_plot <- burden_all$CNV_type
burden_all$CNV_type_plot <- gsub("allCNV","All rCNVs",burden_all$CNV_type_plot)
burden_all$CNV_type_plot <- gsub("del","Deletion",burden_all$CNV_type_plot)
burden_all$CNV_type_plot <- gsub("dup","Duplication",burden_all$CNV_type_plot)
burden_all$CNV_type_plot <- factor(burden_all$CNV_type_plot,level=unique(burden_all$CNV_type_plot))

burden_all$Array <- factor(burden_all$array_sample,level=array_sample$array_sample)
burden_all$data_set <- factor(burden_all$data_set,level=c("ASA","GSAA","Il1M","Om25","PSYC","combined"))
burden_all$CNV_size <- factor(burden_all$CNV_size,level=unique(burden_all$CNV_size))
burden_all$CNV_type <- factor(burden_all$CNV_type,level=unique(burden_all$CNV_type))
burden_all <-
  burden_all %>%
  dplyr::mutate(.filled_gene = NGENE_glm_pval< psignif)

burden_all <-
  burden_all %>%
  dplyr::mutate(.filled_count = NSEG_GENIC_glm_pval< psignif)

burden_all <-
  burden_all %>%
  dplyr::mutate(.filled_length = KB_GENIC_glm_pval< psignif)

# platform plot
psignif <- 0.05
#burden_array <- burden_all[burden_all$CNV_size==">20kb",]
burden_array <- burden_all[burden_all$CNV_type=="allCNV"&burden_all$CNV_size!=">20kb",]

plot_array_gene <- ggplot(data = burden_array, aes(x = NGENE_glm_OR, y = Array)) +
  geom_effect(ggplot2::aes(xmin = NGENE_glm_lowerCI,xmax = NGENE_glm_upperCI,color = Array,filled =.filled_gene),position = ggstance::position_dodgev(height = 0.5)) +
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
  scale_color_d3("category20c") +
  facet_wrap(~CNV_size)

plot_array_counts <- ggplot(data = burden_array, aes(x = NSEG_GENIC_glm_OR, y = Array)) +
  geom_effect(ggplot2::aes(xmin = NSEG_GENIC_glm_lowerCI,xmax = NSEG_GENIC_glm_upperCI,color = Array,filled =.filled_count),position = ggstance::position_dodgev(height = 0.5)) +
  theme_classic() +
  geom_stripes(odd = "#33333333", even = "#00000000") +
  geom_vline(xintercept = 1, linetype=2) +
  xlab("Odds ratio\n95% Confidence Interval") +
  ggtitle("Number of CNVs") +
  ylab("Platform") +
  theme(axis.text.x = element_text(angle = 45, size=8, hjust=1)) +
  theme(legend.position="none") +
  theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 10, face = "plain")) +
  theme(plot.tag = element_text(face='bold')) +
  scale_color_d3("category20c") +
  facet_wrap(~CNV_size)

plot_array_length <- ggplot(data = burden_array, aes(x = KB_GENIC_glm_OR, y = Array)) +
  geom_effect(ggplot2::aes(xmin = KB_GENIC_glm_lowerCI,xmax = KB_GENIC_glm_upperCI,color = Array,filled =.filled_length),position = ggstance::position_dodgev(height = 0.5)) +
  theme_classic() +
  geom_stripes(odd = "#33333333", even = "#00000000") +
  geom_vline(xintercept = 1, linetype=2) +
  xlab("Odds ratio\n95% Confidence Interval") +
  ggtitle("Length of CNVs") +
  ylab("Platform") +
  theme(axis.text.x = element_text(angle = 45, size=8, hjust=1)) +
  theme(legend.position="none") +
  theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 10, face = "plain")) +
  theme(plot.tag = element_text(face='bold')) +
  scale_color_d3("category20c") +
  facet_wrap(~CNV_size)

plot_legend<- ggplot(burden_array, aes(x = data_set, y = KB_GENIC_glm_OR ,color = data_set)) + 
  geom_hline(aes(yintercept = 1, linetype="Odds Ratio = 1")) + 
  #geom_hline(aes(yintercept = 1, linetype="Confidence Interval")) + 
  geom_point(aes(shape="p.value <  0.05"))+
  geom_point(aes(shape="p.value >  0.05"))+
  scale_color_d3("category20c", name="") +
  scale_linetype_manual(values=c("Odds Ratio = 1"="dashed","Confidence Interval"="solid"), name="") +
  scale_shape_manual(values=c("p.value <  0.05" = 16,"p.value >  0.05" = 21), name="") + theme_classic() 

legend <- get_legend(
  # create some space to the left of the legend
  plot_legend+
    theme(legend.position="bottom"))

cnv_event <- cowplot::plot_grid(plot_array_length, plot_array_counts, plot_array_gene+theme(legend.position="none"),NULL,legend,ncol=3,rel_heights = c(5,0.8))
plot_array_gene_final  <- cowplot::plot_grid(plot_array_gene+theme(legend.position="none"),legend,ncol=1,rel_heights = c(5,0.8))
ggsave("output/EF_06_CNV_burden_platform_update.pdf", plot_array_gene_final, width=12, height=5 ,units = "in")
