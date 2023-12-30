#!/bin/bash
#$ -cwd
#$ -l h_vmem=1G
#$ -pe smp 10
#$ -binding linear:10
#$ -l h_rt=10:00:00
#$ -N input_format

permute=${SGE_TASK_ID}
cnv_type=$1
eas_file_wd="/stanley/huang_lab/home/ychen/proj-CNV/07CNV_gene_association/01EAS_gene_assoc/output/permute"
eur_file_wd="/stanley/huang_lab/home/ychen/proj-CNV/07CNV_gene_association/02EUR_gene_assoc/output/permute"
outpur_wd="/stanley/huang_lab/home/ychen/proj-CNV/08CNV_meta_analysis/02gene_assoc_meta/input/permute"

cat ${eas_file_wd}/${SGE_TASK_ID}_${cnv_type}.step2_AFF.regenie |awk 'NR==1{print "GENE BETA SE P"}NR>1{print $3,$9,$10,10^(-$12)}' |awk 'NR==1{print $0,"N A1 A0"}NR>1{print $0,44161,".","."}' >  ${outpur_wd}/EAS.${SGE_TASK_ID}_${cnv_type}.assoc
cat ${eur_file_wd}/${SGE_TASK_ID}_PGC.${cnv_type}.assoc|awk 'NR==1{print $0,"N A1 A0"}NR>1{print $0,34256,".","."}' >  ${outpur_wd}/PGC.${SGE_TASK_ID}_${cnv_type}.assoc

