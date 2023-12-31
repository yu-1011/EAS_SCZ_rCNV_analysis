library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggstance)
library(ggforestplot)
library(ggthemes)
library(ggsci)
library(cowplot)
library(reshape2)

rm(list=ls())

cohort_colors <- c(
  EUR = "#BE1E2D",
  EAS = "#3182BD",
  'Novel Regions' = "#7ec0ee",
  'All Regions' = "#A9A9A9")

wd="/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/02extended_figure/extended_fig5_CNV_burden_pop_com"
setwd(wd)
eur_file_name <- list.files("input/",pattern="*EUR_platform.*.burden$")
eas_file_name <- list.files("input/",pattern="*EAS_platform.*.burden$")

burden_eur <- read.table(paste0("input/",eur_file_name),header=T)
burden_eas <- read.table(paste0("input/",eas_file_name),header=T)

burden_eas$data_set <- gsub("combined","EAS",burden_eas$data_set)
idx <- intersect(names(burden_eur),names(burden_eas))
burden_pop <- rbind(burden_eur[burden_eur$data_set=="EUR",idx],burden_eas[burden_eas$data_set=="EAS",idx])

burden_pop$plot_color <- factor(burden_pop$data_set,level=c("EAS","EUR"))
burden_pop$CNV_size <- factor(burden_pop$CNV_size,level=unique(burden_pop$CNV_size))
burden_pop$CNV_type <- gsub("allCNV","All rCNVs",burden_pop$CNV_type)
burden_pop$CNV_type <- gsub("del","Deletion",burden_pop$CNV_type)
burden_pop$CNV_type <- gsub("dup","Duplication",burden_pop$CNV_type)

psignif <- 0.05
burden_pop <-
  burden_pop %>%
  dplyr::mutate(.filled_gene = NGENE_glm_pval < psignif)

burden_pop <-
  burden_pop %>%
  dplyr::mutate(.filled_count = NSEG_GENIC_glm_pval < psignif)

burden_pop <-
  burden_pop %>%
  dplyr::mutate(.filled_KB = KB_GENIC_glm_pval < psignif)

burden_pop$CNV_size <- gsub("kb","",burden_pop$CNV_size)
burden_pop$CNV_size <- factor(burden_pop$CNV_size,levels = unique(burden_pop$CNV_size))
#burden_pop <- burden_pop[burden_pop$CNV_type=="All rCNVs",]

plot_pop_gene <- ggplot(data = burden_pop, aes(x = NGENE_glm_OR, y = CNV_size)) +
  geom_effect(ggplot2::aes(xmin = NGENE_glm_lowerCI,xmax = NGENE_glm_upperCI,color = plot_color, filled =.filled_gene),position = ggstance::position_dodgev(height = 0.5)) +
  theme_classic() +
  geom_stripes(odd = "#33333333", even = "#00000000") +
  geom_vline(xintercept = 1, linetype=2) +
  xlab("OR per gene impacted by rCNV") +
  #ggtitle("Number of gene affected by CNVs") +
  ylab("rCNV Categories (kb)") +
  theme(axis.text.x = element_text(angle = 45, size=8, hjust=1)) +
  theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 10, face = "plain")) +
  theme(plot.tag = element_text(face='bold')) +
  scale_color_manual(values = cohort_colors, name="") +
  facet_grid(~CNV_type) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
  theme(legend.position="bottom") 

plot_pop_count <- ggplot(data = burden_pop, aes(x = NSEG_GENIC_glm_OR, y = CNV_size)) +
  geom_effect(ggplot2::aes(xmin = NSEG_GENIC_glm_lowerCI,xmax = NSEG_GENIC_glm_upperCI,color = plot_color, filled =.filled_count),position = ggstance::position_dodgev(height = 0.5)) +
  theme_classic() +
  geom_stripes(odd = "#33333333", even = "#00000000") +
  geom_vline(xintercept = 1, linetype=2) +
  xlab("OR per copy of rCNV") +
  #ggtitle("Number of gene affected by CNVs") +
  ylab("rCNV Categories (kb)") +
  theme(axis.text.x = element_text(angle = 45, size=8, hjust=1)) +
  theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 10, face = "plain")) +
  theme(plot.tag = element_text(face='bold')) +
  scale_color_manual(values = cohort_colors, name="") +
  facet_grid(~CNV_type) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
  theme(legend.position="bottom") 

plot_pop_kb <- ggplot(data = burden_pop, aes(x = KB_GENIC_glm_OR, y = CNV_size)) +
  geom_effect(ggplot2::aes(xmin = KB_GENIC_glm_lowerCI,xmax = KB_GENIC_glm_upperCI,color = plot_color, filled =.filled_KB),position = ggstance::position_dodgev(height = 0.5)) +
  theme_classic() +
  geom_stripes(odd = "#33333333", even = "#00000000") +
  geom_vline(xintercept = 1, linetype=2) +
  xlab("OR per 100kb of rCNV") +
  #ggtitle("Number of gene affected by CNVs") +
  ylab("rCNV Categories (kb)") +
  theme(axis.text.x = element_text(angle = 45, size=8, hjust=1)) +
  theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 10, face = "plain")) +
  theme(plot.tag = element_text(face='bold')) +
  scale_color_manual(values = cohort_colors, name="") +
  facet_grid(~CNV_type) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
  theme(legend.position="bottom") 

legend <- get_legend(
  # create some space to the left of the legend
  plot_pop_kb+
    theme(legend.position="bottom"))

cnv_event <- cowplot::plot_grid(plot_pop_kb+theme(legend.position="none"), plot_pop_count+theme(legend.position="none"), plot_pop_gene+theme(legend.position="none"),labels = 'auto', ncol=3, rel_widths = c(0.3,0.3,0.3))
cnv_event_final <- cowplot::plot_grid(cnv_event, legend, ncol=1, rel_heights = c(5,0.5))
ggsave("output/EF08_CNV_burden_population_comparison.pdf", cnv_event_final, width=10, height=6 ,units = "in")
