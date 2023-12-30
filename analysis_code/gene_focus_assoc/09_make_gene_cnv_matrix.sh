#!/bin/bash
#$ -cwd
#$ -pe smp 4 -binding linear:4
#$ -l h_rt=12:00:00,h_vmem=15G,
#$ -N gene_assoc_gene_cnv_matrix

source /broad/software/scripts/useuse
use R-4.0

input=$1
phefile=$2
region_file_all=$3
output=$4


cat ${region_file_all} | awk 'NR=="'"${SGE_TASK_ID}"'"{print}' > /stanley/huang_lab/home/ychen/proj-CNV/07CNV_gene_association/01EAS_gene_assoc/tmp/${SGE_TASK_ID}_region_tmp.txt
region_file="/stanley/huang_lab/home/ychen/proj-CNV/07CNV_gene_association/01EAS_gene_assoc/tmp/${SGE_TASK_ID}_region_tmp.txt"

wdir="/stanley/huang_lab/home/ychen/proj-CNV/07CNV_gene_association/01EAS_gene_assoc"
gene_file="/stanley/huang_lab/home/ychen/proj-CNV/misc/ucsc.hg19.knowgene.exon.plink.txt"

Rscript /stanley/huang_lab/home/ychen/proj-CNV/07CNV_gene_association/01EAS_gene_assoc/scr/07_make_gene_cnv_matrix.r ${wdir} ${input} ${gene_file} ${phefile} ${region_file} ${output}
