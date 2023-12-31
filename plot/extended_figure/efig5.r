setwd("/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/02extended_figure/extended_fig2_CNV_QC/pipeline_evaluation/")
rm(list=ls())
library("ggplot2")
library("cowplot")

data <- read.table("input/pipeline_comarison.txt",header=T)
data$cnv <- paste0(data$chr,"*",data$start,"*",data$end)
data <- data[!duplicated(data$cnv),]
data$Status <- ifelse(data$status=="detect","Detected","Not detected")

length_test <- wilcox.test(data[data$status=="detect","length"],data[data$status=="not_detect","length"])
p1 <- ggplot(data, aes(x = Status , y = as.numeric(length)/1000, fill = Status)) +
  geom_boxplot() +
  #geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  scale_y_log10() +
  labs(x = "rCNV type",
       y = "rCNV Length (kb)") +
  geom_text(aes(x=1.5, y=max(as.numeric(length)/1000) * 1.05, label=paste0("p-value = ",signif(length_test$p.value,2)))) +
  theme_minimal()

p2 <- ggplot(data, aes(x = status , y = SCORE, fill = Status)) +
  geom_violin() +
  #geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  #scale_y_log10() +
  labs(x = "rCNV type",
       y = "rCNV SCORE") +
  theme_minimal()
wilcox.test(data[data$status=="detect","SCORE"],data[data$status=="not_detect","SCORE"])

sites <- wilcox.test(data[data$status=="detect","SITES"],data[data$status=="not_detect","SITES"])
p3 <- ggplot(data, aes(x = Status , y = SITES, fill = Status)) +
  geom_boxplot() +
  #geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  scale_y_log10() +
  labs(x = "rCNV type",
       y = "rCNV SITES") +
  geom_text(aes(x=1.5, y=max(SITES) * 1.05, label=paste0("p-value = ",signif(sites$p.value,2)))) +
  theme_minimal()

p_final <- plot_grid(p1,p3,labels = "auto")
ggsave( "output/EF_05_pipeline_comarison.pdf",p_final, width=12, height=6)
