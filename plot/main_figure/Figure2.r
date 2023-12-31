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
  #EAS = "#0087FF",
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
