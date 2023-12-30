#!/bin/bash
#$ -cwd
#$ -l h_vmem=1G
#$ -pe smp 10
#$ -binding linear:10
#$ -l h_rt=10:00:00
#$ -N meta_analysis

cd /stanley/huang_lab/home/ychen/proj-CNV/08CNV_meta_analysis/02gene_assoc_meta/
cnv_type=$1
#~/software/generic-metal/metal metaanalysis_script.txt

echo "################################################
SCHEME                 SAMPLESIZE
LOGPVALUE              OFF

## setting in metal
# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
MARKER GENE
ALLELE A0 	A1
EFFECT BETA
STDERR SE
PVALUE P 
WEIGHT N
PROCESS input/EAS.${cnv_type}.assoc


# === THE SECOND INPUT FILE HAS THE SAME FORMAT AND CAN BE PROCESSED IMMEDIATELY ===
PROCESS input/PGC.${cnv_type}.assoc

# === CARRY OUT AN INTERIM ANALYSIS OF THE FIRST FOUR FILES ===
OUTFILE output/META_PGC_EAS.${cnv_type}_ .txt
ANALYZE HETEROGENEITY
QUIT
" > tmp/tmp_${cnv_type}.txt

~/software/generic-metal/metal tmp/tmp_${cnv_type}.txt
