#!/bin/bash
#$ -cwd
#$ -pe smp 4 -binding linear:4
#$ -l h_rt=12:00:00,h_vmem=15G,
#$ -N ped_to_bed

source /broad/software/scripts/useuse
use PLINK

WORKING_DIR=/stanley/huang_lab/home/ychen/proj-CNV/07CNV_gene_association/01EAS_gene_assoc/input/cnv_bed
cnv_type=$1

plink \
--file ${WORKING_DIR}/cnv.${cnv_type} \
--make-bed \
--out ${WORKING_DIR}/cnv.${cnv_type}
