#!/bin/bash -e
#$ -S /bin/bash
#$ -N threaded16_anchovy
#$ -M Patrick.Dolan@nih.gov
#$ -m be
#$ -l h_vmem=15G
#$ -cwd
#$ -o anchOut/
#$ -pe threaded 16

module purge
module load python
for i in `ls ${1}`
do
python ~/lab_share/anchovy/anchovy.0.0.py ${1}/${i}/merge.sam ~/lab_share/anchovy/3M-february-2018.txt CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNTTTCTTATAT
done
