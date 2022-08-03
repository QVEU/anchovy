#!/bin/bash -e
#$ -S /bin/bash
#$ -N threaded16_v1_anchovy
#$ -M Patrick.Dolan@nih.gov
#$ -m be
#$ -l h_vmem=15G
#$ -cwd
#$ -o anchOut/
#$ -pe threaded 16

pip install pysam
pip install pandas

template='/hpcdata/lvd_qve/QVEU_Code/sequencing/template_fastas/cvb3_rna.fa'
indir='/hpcdata/lvd_qve/Sequencing_Data/QVEU_Seq_0010_Minion_ndas10xlib2/no_sample/20220727_2153_MC-113212_FAT27957_6eb71326/fastq_pass/'

echo `ls -d ${indir}/*`

for BC in `ls -d ${indir}/*`
do

module purge
module load foss/2021a
module load python/3.7.3
module load Biopython/1.79-foss-2021a

#echo $indir
echo BARCODE: $BC
python ~/lab_share/anchovy/anchovy.0.1.py ${BC}/merge.sam /hpcdata/lvd_qve/anchovy/3M-february-2018.txt CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNTTTCTTATAT

#make fas for each cell with barcodes as headers
echo "Generating fastas for each cell..."
python ~/lab_share/anchovy/CellConsensus.py ${BC}
#Map reads in fastq files  and generate sam file per cell fa file.

echo "Mapping each fa..."
for i in `ls ${BC}/*.fa`
do
  module load minimap2
  echo $i
  echo "minimap..."
  minimap2 -ax map-ont $template $i > ${i/\.fa/}.sam
  echo ${i/\.fa/}.sam
  samtools view -S -b ${i/\.fa/}.sam > ${i/\.fa/}.bam # for Oxford Nanopore reads

echo "Changing Modules... $BC"
module purge
module load samtools
#sort and pile up reads. make bed cov file with bed files from template.
echo "Making cell pileups"
samtools sort ${i/\.fa/}.bam > ${i/\.fa/}_sort.bam
samtools mpileup ${i/\.fa/}_sort.bam > ${i/\.fa/}_sort_pile.pile

rm ${i/\.fa/}.sam
rm ${i/\.fa/}.bam
rm ${i/\.fa/}_sort.bam

echo "Done. "
done

done

module purge
