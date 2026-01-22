#!/bin/bash
#$ -S /bin/bash
#$ -N AnchovyJob_050323
#$ -M Patrick.Dolan@nih.gov
#$ -m be
#$ -l h_vmem=25G
#$ -cwd

module purge

whitelist=$1
indir=$2
outdir=$indir
template=$3
infile=$4

#echo `ls -d ${indir}/*`
if test -f "${indir}/${infile}"
  then
      echo "${indir}/${infile} exists, running anchovy..."
      module load python/3.9.5-GCCcore-10.3.0
      module load foss/2021a
      module load Biopython/1.79-foss-2021a
      pip install Levenshtein
      pip install pysam
      pip install pandas

      #echo $indir
      if test -f "${indir}/${infile/\.sam/}_anchovy_v2.csv"
      then
        echo "><> Anchovy already canned. ><>"
        echo "${indir}/${infile/\.sam/}_anchovy_v2.csv"
      else
        echo "Running Anchovy v2..."
        python /hpcdata/lvd_qve/QVEU_Code/anchovy/anchovy.py ${indir}/${infile} $whitelist CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNTTTCTTATAT
      fi
      #make fas for each cell with barcodes as headers
      echo "Generating fastas for each cell..."
      python /hpcdata/lvd_qve/QVEU_Code/anchovy/CBCtoFasta.py ${indir} ${infile/\.sam/}_anchovy_v2.csv #edited 5/4
      echo "Done."

      #Map reads in fastq files  and generate sam file per cell fa file.
      echo "Mapping each fa..."

      #for each barcode fasta...
      for i in $(ls ${indir}/*.fa)
      do
        module load minimap2
        echo "minimap... $i"
        minimap2 -ax map-ont $template $i > ${i/\.fa/}_BC.sam  # Map NANOPORE reads -PD 1/11/24
        echo "Prepared ${i/\.fa/}.sam."

        echo "Changing Modules... $2"
        module purge
        module load samtools

        #sort and pile up reads. make bed cov file with bed files from template.
        samtools view -S -b ${i/\.fa/}_BC.sam > ${i/\.fa/}_BC.bam # Convert sam to bam
        echo "Making cell pileups"
        samtools sort ${i/\.fa/}_BC.bam > ${i/\.fa/}_BC_sort.bam

        python /hpcdata/lvd_qve/QVEU_Code/sam2consensus/sam2consensus.py -c .5 -m 5 --outfolder $outdir -i ${i/\.fa/}_BC.sam #changed
        samtools depth -d 100000 ${i/\.fa/}_BC_sort.bam > ${i/\.fa/}_BC_sort.depth
        #Clean up all the files (there are so many files.)
        rm ${i/\.fa/}_BC.bam
        rm ${i/\.fa/}_BC.sam
        rm ${i/\.fa/}_BC_sort.bam
        rm $i #remove fa file
        echo "Done. "
      done
    cat ${indir}/*.fasta > ${indir}_allConsensus.fasta
    echo "Consensus files merged."

  else
    echo "${indir}/merge.sam does not exist."
fi
module purge
