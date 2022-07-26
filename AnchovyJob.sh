#!/bin/bash -e
#$ -S /bin/bash
#$ -N threaded16_stickle
#$ -M Patrick.Dolan@nih.gov
#$ -m be
#$ -l h_vmem=32G
#$ -cwd
#$ -o anchOut/
#$ -pe threaded 16

module purge
module load python
for i in `ls /nethome/dolanpt/lab_share/Sequencing_Data/QVEU_Seq_0008_Minion_10xlib_7-12-2022/no_sample/20220712_2344_MC-113212_FAT05581_3e3d827f/fastq_pass/`
do
python ~/lab_share/anchovy/anchovy.0.0.py /nethome/dolanpt/lab_share/Sequencing_Data/QVEU_Seq_0008_Minion_10xlib_7-12-2022/no_sample/20220712_2344_MC-113212_FAT05581_3e3d827f/fastq_pass/${i}/merge.sam ~/lab_share/anchovy/3M-february-2018.txt CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNTTTCTTATAT
done
