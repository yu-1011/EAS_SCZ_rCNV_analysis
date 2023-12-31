library(GenomicRanges)
wd="/Users/yu/Desktop/manuscript_ongoing/proj-CNV/figure/02extended_figure/extended_fig6_CNV_gene_assoc_pop_com/"
setwd(wd)
rm(list=ls())

# 读取文件
file1 <- read.table("input/EAS_CNV_region.txt", header = FALSE)
file2 <- read.table("input/EUR_CNV_region.txt", header = FALSE)

# 创建GenomicRanges对象

file1 <- file1[file1$V1 %in% c(1:22), ]
file2 <- file2[file2$V1 %in% c(1:22), ]

gr1 <- GRanges(seqnames = file1$V1, ranges = IRanges(start = file1$V2, end = file1$V3))
gr2 <- GRanges(seqnames = file2$V1, ranges = IRanges(start = file2$V2, end = file2$V3))

seqlevelsStyle(gr1) <- "UCSC"
seqlevelsStyle(gr2) <- "UCSC"

# 找到overlap的区域
overlaps <- findOverlaps(gr1, gr2)

# 手动计算并过滤至少有50% overlap的区域
hits <- subjectHits(overlaps)
widths <- width(pintersect(gr1[queryHits(overlaps)], gr2[hits]))

# 过滤结果
filtered_indices <- which(widths >= 0.5 * width(gr1[queryHits(overlaps)]))
filtered_overlaps <- overlaps[filtered_indices]

# 输出结果
filtered_overlaps

# 获取overlap的区域
overlap_regions <- pintersect(gr1[queryHits(filtered_overlaps)], gr2[subjectHits(filtered_overlaps)])
overlap_regions <- pintersect(gr1[queryHits(overlaps)], gr2[subjectHits(overlaps)])

# 合并重叠的区域
unique_overlap_regions <- reduce(overlap_regions)

# 转换为数据框
unique_overlap_df <- data.frame(
  chromosome = seqnames(unique_overlap_regions),
  start = start(unique_overlap_regions),
  end = end(unique_overlap_regions)
)

# 打印数据框
print(unique_overlap_df)

# 保存到文件
write.table(unique_overlap_df, "unique_overlap_all_regions.txt", sep = "\t", row.names = FALSE, quote = FALSE)
