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

wd="/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/01main_figure/figure2_burden_test/"
setwd(wd)
loci_file_name <- list.files("input/",pattern="*novel_loci.*.burden$")
eur_file_name <- list.files("input/",pattern="*EUR_platform.*.burden$")
eas_file_name <- list.files("input/",pattern="*EAS_platform.*.burden$")

burden_loci <- read.table(paste0("input/",loci_file_name),header=T)
burden_eur <- read.table(paste0("input/",eur_file_name),header=T)
burden_eas <- read.table(paste0("input/",eas_file_name),header=T)
burden_loci <- rbind(burden_eas,burden_loci)

burden_eas$data_set <- gsub("combined","EAS",burden_eas$data_set)
idx <- intersect(names(burden_eur),names(burden_eas))
burden_pop <- rbind(burden_eur[burden_eur$data_set=="EUR",idx],burden_eas[burden_eas$data_set=="EAS",idx])

psignif <- 0.05
condition1 <- burden_loci$CNV_type == "allCNV" & (burden_loci$region_set == "allregions" | burden_loci$region_set == "novelregions")
condition2 <- burden_loci$CNV_type == "del" & (burden_loci$region_set == "allregions" | burden_loci$region_set == "novelregions")
condition3 <- burden_loci$CNV_type == "dup" & (burden_loci$region_set == "allregions" | burden_loci$region_set == "novelregions")
burden_loci <- burden_loci[condition1 | condition2 | condition3, ]

burden_loci$CNV_size <- factor(burden_loci$CNV_size,level=unique(burden_loci$CNV_size))
burden_loci$region_set <- factor(burden_loci$region_set,level=unique(burden_loci$region_set))
burden_loci$CNV_type <- gsub("allCNV","All rCNVs",burden_loci$CNV_type)
burden_loci$CNV_type <- gsub("del","Deletion",burden_loci$CNV_type)
burden_loci$CNV_type <- gsub("dup","Duplication",burden_loci$CNV_type)
burden_loci$region_set <- gsub("allregions","All Regions",burden_loci$region_set)
burden_loci$region_set <- gsub("novelregions.*","Novel Regions",burden_loci$region_set)
burden_loci$plot_color <- burden_loci$region_set

burden_loci <-
  burden_loci %>%
  dplyr::mutate(.filled_gene = NGENE_glm_pval< psignif)

burden_loci <-
  burden_loci %>%
  dplyr::mutate(.filled_count = NSEG_GENIC_glm_pval< psignif)

burden_loci <-
  burden_loci %>%
  dplyr::mutate(.filled_length = KB_glm_pval< psignif)

burden_loci$CNV_size <- gsub("kb","",burden_loci$CNV_size)
burden_loci$CNV_size <- factor(burden_loci$CNV_size,levels = unique(burden_loci$CNV_size))

#figure 2a
burden_combine <- burden_loci[burden_loci$data_set=="combined"&burden_loci$region_set=="All Regions",]
classify_gene <- reshape2::dcast(burden_combine, CNV_size + .filled_gene ~ CNV_type , value.var = "NGENE_glm_OR")
classify_gene <- melt(classify_gene)
cnv_event_2a_gene <- ggplot(classify_gene, aes(x = CNV_size, y = round(value,digits = 3),  group = variable,col = variable, fill = variable)) + 
  geom_line(stat = "identity") +  
  geom_point(aes(filled =.filled_gene)) +
  theme_minimal() +  
  scale_color_d3("category20c", name="") +
  scale_fill_d3("category20c", name="") +
  ylim(c(1,1.1)) +
  theme(legend.position="bottom") +
  labs(title = "", x = "", y = "OR per gene impacted by rCNV")  

classify_count <- reshape2::dcast(burden_combine, CNV_size + .filled_count ~ CNV_type, value.var = "NSEG_GENIC_glm_OR")
classify_count <- melt(classify_count)
cnv_event_2a_count <- ggplot(na.omit(classify_count), aes(x = CNV_size, y = value,  group = variable,col = variable, fill = variable)) + 
  geom_line(stat = "identity") +  
  geom_point(aes(filled =.filled_count)) +
  #geom_point(data=na.omit(classify_count[classify_count$.filled_count==T,])) +
  #geom_point(data=na.omit(classify_count[classify_count$.filled_count==F,]),aes(x = CNV_size, y = value, col = variable),shape=21,fill = "white" ) +
  theme_minimal() +  
  scale_color_d3("category20c", name="") +
  scale_fill_d3("category20c", name="") +
  ylim(c(0.99,2.1)) +
  theme(legend.position="bottom") +
  labs(title = "",x = "CNV Categories (kb)", y = "OR per copy of rCNV")  
 
classify_kb <- reshape2::dcast(burden_combine, CNV_size + .filled_length  ~ CNV_type , value.var = "KB_glm_OR")
classify_kb <- melt(classify_kb)
cnv_event_2a_kb <- ggplot(na.omit(classify_kb), aes(x = CNV_size, y = value,  group = variable,col = variable, fill = variable)) + 
  #geom_area(alpha = 0.4, position = 'identity') +  
  geom_line(stat = "identity") + 
  geom_point(aes(filled =.filled_length)) +
  #geom_point(data=na.omit(classify_kb[classify_kb$.filled_length==T,])) +
  #geom_point(data=na.omit(classify_kb[classify_kb$.filled_length==F,]),aes(x = CNV_size, y = value, col = variable),shape=21,fill = "white" ) +
  theme_minimal() +  
  scale_color_d3("category20c", name="") +
  scale_fill_d3("category20c", name="") +
  ylim(c(1,1.05)) +
  theme(legend.position="bottom") +
  labs(title = "", x = "" , y = "OR per 100kb of rCNV")

legend <- get_legend(
  # create some space to the left of the legend
  cnv_event_2a_count+
    theme(legend.position="bottom"))

cnv_event_2a <- cowplot::plot_grid(cnv_event_2a_kb+theme(legend.position="none"), cnv_event_2a_count+theme(legend.position="none"), cnv_event_2a_gene+theme(legend.position="none"),labels = NULL, ncol=3, rel_widths = c(0.3,0.3,0.3))
cnv_event_2a_final <- cowplot::plot_grid(cnv_event_2a, legend, ncol=1, rel_heights = c(5,0.5))

#figure 2b
burden_pop$plot_color <- factor(burden_pop$data_set,level=c("EAS","EUR"))
burden_pop$CNV_size <- factor(burden_pop$CNV_size,level=unique(burden_pop$CNV_size))
burden_pop$CNV_type <- gsub("allCNV","All rCNVs",burden_pop$CNV_type)
burden_pop$CNV_type <- gsub("del","Deletion",burden_pop$CNV_type)
burden_pop$CNV_type <- gsub("dup","Duplication",burden_pop$CNV_type)

burden_pop <-
  burden_pop %>%
  dplyr::mutate(.filled_gene = NGENE_glm_pval < psignif)

burden_pop <-
  burden_pop %>%
  dplyr::mutate(.filled_count = NSEG_GENIC_glm_pval < psignif)

burden_pop$CNV_size <- gsub("kb","",burden_pop$CNV_size)
burden_pop$CNV_size <- factor(burden_pop$CNV_size,levels = unique(burden_pop$CNV_size))

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
  theme(legend.position="bottom") 

#figure 2c
plot_loci_gene <- ggplot(data = burden_loci[burden_loci$data_set=="combined",], aes(x = NGENE_glm_OR, y = CNV_size)) +
  geom_effect(ggplot2::aes(xmin = NGENE_glm_lowerCI,xmax = NGENE_glm_upperCI,color = plot_color,filled =.filled_gene),position = ggstance::position_dodgev(height = 0.75)) +
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
  theme(legend.position="bottom") +
  facet_grid(~CNV_type)

plot_legend <- ggplot(burden_combine, aes(x = NSEG_GENIC_glm_pval, y = NSEG_GENIC_glm_pval)) + 
  #geom_hline(aes(yintercept = 1, linetype="Odds Ratio = 1",color = plot_color)) + 
  geom_hline(aes(yintercept = 1, linetype= "95% Confidence Interval")) + 
  geom_point(aes(shape="p-value <  0.05"))+
  geom_point(aes(shape="p-value >  0.05"))+
  #scale_color_d3("category20c", name="") + 
  scale_color_manual(values = cohort_colors, name="") +
  scale_linetype_manual(values=c("Odds Ratio = 1"="dashed","95% Confidence Interval"="solid"), name="") +
  scale_shape_manual(values=c("p-value <  0.05" = 16,"p-value >  0.05" = 21), name="") + theme_classic() 

legend_final <- get_legend(
  # create some space to the left of the legend
  plot_legend +
  theme(legend.position="bottom"))

cnv_event_1 <- cowplot::plot_grid(plot_pop_gene, plot_loci_gene, ncol=2, rel_heights = c(5,0.5), labels = c("b","c"))
cnv_event_2 <- cowplot::plot_grid(cnv_event_2a_final,cnv_event_1, nrow=2, rel_heights = c(5,5,0.5), labels = c("a",""))
plot_final <- cowplot::plot_grid(cnv_event_2, legend_final, ncol=1, rel_heights = c(5,0.5))
ggsave("output/Figure2_CNV_burden_combined.pdf", plot_final, width=10, height=6 ,units = "in")
