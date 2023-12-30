#!/bin/bash
#$ -l h_vmem=1G
#$ -pe smp 10
#$ -binding linear:10
#$ -l h_rt=10:00:00

source /broad/software/scripts/useuse
use PLINK2

WORKING_DIR=/stanley/huang_lab/home/ychen/proj-CNV/07CNV_gene_association/01EAS_gene_assoc/
phe_file=${WORKING_DIR}/input/cnv_bed/all_cohorts.reformatted.phe
cnv_type=$1

# plink2 \
#   --pfile ${WORKING_DIR}/regenie/all_genotype.pruned \
#   --remove ${WORKING_DIR}/regenie/dup_samples.ids \
#   --out ${WORKING_DIR}/regenie/all_genotype.pruned \
#   --make-pgen \
#   --geno 0.001 \
#   --extract ${WORKING_DIR}/regenie/all_genotype.pruned.prune.in \
#   --chr 1-22 \
#   --maf 0.05 \
#   --indep-pairwise 200 50 0.2 \

# use .regenie-2.2.4
# regenie \
#   --step 1 \
#   --bt \
#   --firth \
#   --pThresh 0.05 \
#   --pgen ${WORKING_DIR}/regenie/all_genotype.pruned \
#   --covarFile ${phe_file} \
#   --phenoFile ${phe_file} \
#   --bsize 100 \
#   --loocv \
#   --threads 10 \
#   --phenoColList AFF \
#   --covarColList SEX,C1,C2,C3,C4,C5 \
#   --catCovarList CNV_platform \
#   --out ${WORKING_DIR}/regenie/step1 \

use .regenie-3.2.2
regenie \
  --step 2 \
  --bt \
  --bed ${WORKING_DIR}/input/cnv_bed/cnv.${cnv_type} \
  --covarFile ${phe_file} \
  --phenoFile ${phe_file} \
  --bsize 2 \
  --firth \
  --pThresh 0.05 \
  --threads 10 \
  --phenoColList AFF \
  --covarColList SEX,C1,C2,C3,C4,C5 \
  --pred ${WORKING_DIR}/input/step1_pred.list \
  --out ${WORKING_DIR}/output/${cnv_type}.step2 \
  --catCovarList CNV_platform \
