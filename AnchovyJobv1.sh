#!/bin/bash -e
#$ -S /bin/bash
#$ -N threaded16_v1_anchovy
#$ -M Patrick.Dolan@nih.gov
#$ -m be
#$ -l h_vmem=15G
#$ -cwd
#$ -o anchOut/
#$ -pe threaded 16

module purge
module load python
indir='/hpcdata/lvd_qve/Sequencing_Data/QVEU_Seq_0010_Minion_ndas10xlib2/no_sample/20220727_2153_MC-113212_FAT27957_6eb71326/fastq_pass/'
for i in `ls ${indir}`
do
python ~/lab_share/anchovy/anchovy.0.1.py ${indir}/${i}/merge.sam ~/lab_share/anchovy/3M-february-2018.txt CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNTTTCTTATAT
done
