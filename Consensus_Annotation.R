# Consensus_Plots_v3.R
# Analysis and annotation of genotypes from `anchovy`.
#
# This script analyzes consensus genotypes derived from single-cell viral sequencing.
# It annotates mutations relative to a provided reference consensus sequence,
# computes haplotype frequencies across cells, and can generate networks and plots.
#
# Usage (example):
# Rscript Consensus_Annotation_v3.R /path/to/reference_consensus.txt /path/to/filtConsensus.csv /path/to/output/prefix
#
# The script expects:
# - A plaintext reference consensus sequence file (single contiguous sequence).
# - A CSV of filtered consensus haplotypes (from ConsensusTools) with at least columns:
#   CBC_ID (cell barcode / cell identifier), genotype (mutations encoded like "12A_45T"), description (optional).
# - A NAME/prefix used to name output files.
#
# Output:
# - NAME_annot_v3.csv : annotated mutation table per haplotype
# - NAME_epistaticNetwork.csv and NAME_genotypeNetwork.csv (if network=TRUE)
# - NAME.pdf (if plothaps=TRUE) with haplotype/position visualizations
#
# NOTE: This file contains only inline comments and does not change behavior.

# Load libraries required for sequence handling, data manipulation, and plotting.
# DECIPHER is used for translation (codon->AA).
library(data.table,quietly = T)
library(Biostrings,quietly = T)
library(reshape2,quietly = T)
library(ggplot2,quietly = T)
library(DECIPHER,quietly = T)
library(cowplot,quietly = T)
library(purrr,quietly = T)
library(ggrepel,quietly = T)

# FUNCTIONS

# haploanalysis:
# Input: data.table `cons` with at least columns CBC_ID (cell identifier) and genotype (string of mutation tokens separated by "_").
# Purpose: Unroll genotype strings into individual mutation records, compute per-mutation counts and frequencies across cells,
#          and add an explicit "reference" row for cells that match the reference.
# Returns: data.table `mutantTable` with columns including mutants (token like "12A"), pos (position numeric),
#          base (mutant base), CBC_ID (cell), genotype, count, freq, total, and BCMutCount.
haploanalysis <- function(cons) {
  DEPTH = length(unique(cons$CBC_ID))
  print(paste("Number of Cells:", DEPTH))
  
  # Mark entries with an empty genotype as wildtype (WT == TRUE) per CBC_ID
  cons[genotype=="",WT:=T,by=CBC_ID]
  
  # Identify reference haplotypes (cells flagged WT)
  referenceHaplotypes<-unique(cons[WT==T,CBC_ID])
  
  # Split genotype strings on "_" to get individual mutation tokens (unrolled)
  unrolled = cons[, c(strsplit(genotype, "_")), by = CBC_ID]
  unrolledgeno = cons[, genotype, by = CBC_ID]
  unrolled <- merge(unrolled, unrolledgeno)
  
  colnames(unrolled) <- c("CBC_ID", "mutants", "genotype")
  
  # Extract base letter from token (last character of token) and convert to uppercase.
  bases = unrolled[, toupper(strsplit(mutants, "")[[1]][length(strsplit(mutants, "")[[1]])]), by = mutants]
  
  # Extract numeric position from token (all but last character(s) of token joined together).
  positions = unrolled[, as.integer(paste(collapse = "", strsplit(mutants, "")[[1]][1:(length(strsplit(mutants, "")[[1]]) - 1)])), by = mutants]
  
  colnames(positions) <- c("mutants", "pos")
  colnames(bases) <- c("mutants", "base")
  
  mutantKey = merge.data.table(positions, bases)
  
  # Remove tokens that indicate deletions represented by "-" base
  mutantKey<-mutantKey[base!="-"]
  
  # Merge to get per-cell-per-mutation table
  mutantTable = merge.data.table(mutantKey, unrolled, by = "mutants")
  
  # Compute number of mutations per CBC_ID (filter out cells with extremely many mutations)
  mutantTable[, BCMutCount := length(mutants), by = "CBC_ID"]
  mutantTable <- mutantTable[BCMutCount < 200]
  
  # Add explicit rows for reference cells (cells with no called mutations) so they appear in downstream summaries
  refTable<-data.table(mutants="",pos=NA,base=NA,CBC_ID=referenceHaplotypes,genotype="reference",BCMutCount=0)
  mutantTable<-rbindlist(fill = T,use.names = T,list(mutantTable,refTable))
  print(mutantTable)
  
  # Compute overall frequency of each mutation across cells
  mutantTable[, total := DEPTH]
  mutantTable[, count := nrow(.SD), by = mutants]
  mutantTable[, freq := count / total]
  print(mutantTable)
  return(mutantTable)
}

# hapNetworkGen:
# Input: haplocounts (output of haploanalysis) and NAME (prefix for saving CSV).
# Purpose: Generate binary presence/absence and count matrices for each genotype vs mutation token.
#          Compute pairwise genotype overlap (shared mutations) and export network CSVs that describe
#          genotype relationships and single-step epistatic edges.
# Returns: binaryMatrix (data.table with genotype rows and mutation columns; 1 if mutation present).
hapNetworkGen <- function(haplocounts, NAME) {
  # Create a binary matrix: rows=genotype, cols=mutants; 1 if the mutation is present in genotype
  binaryMatrix = dcast.data.table(
    haplocounts,
    genotype ~ mutants,
    fill = 0,
    value.var = "freq",
    fun.aggregate = function(X)
      ifelse(X > 0, 1, 0)
  )
  
  # Create a count matrix: how many cells (count) in each genotype contain each mutation
  countMatrix = dcast.data.table(haplocounts,
                                 genotype ~ mutants,
                                 fill = 0,
                                 value.var = "count")
  
  # Summarize counts into a per-genotype count value (used as metadata for edges)
  counts = data.table(countMatrix$genotype,
                      rowSums(countMatrix[, -1]) / rowSums(binaryMatrix[, -1]))
  colnames(counts) <- c("genotype", "count")
  
  allentries = NULL
  genoList <- as.factor(binaryMatrix$genotype)
  # Compare each genotype to every other genotype to record overlap in mutations
  for (i in as.factor(binaryMatrix$genotype)) {
    print(i)
    for (j in genoList) {
      Source = t(binaryMatrix[genotype == i])
      Target = t(binaryMatrix[genotype == j])
      comparison <- data.frame(Target, Source)
      # fcomp picks rows where both source and target have mutation present
      fcomp <- comparison[Target == 1 & Source == 1, ]
      if (nrow(fcomp) > 0) {
        entry = data.table(
          source = i,
          target = j,
          overlap = nrow(fcomp),
          mutNumSource = length(strsplit(i, "_")[[1]]),
          mutNumTarget = length(strsplit(j, "_")[[1]])
        )
        allentries <- rbindlist(list(allentries, entry))
      }
    }
    # Avoid duplicate comparisons by removing current genotype from future comparisons
    genoList <-
      genoList[genoList != i]
  }
  # Attach genotype counts metadata to the overlap table
  allentries <-
    merge.data.table(counts, allentries, by.x = "genotype", by.y = "source")
  
  # singleSteps selects pairs differing by exactly one mutation (adjacent in the genotype lattice)
  singleSteps = allentries[((mutNumSource == overlap) & (mutNumTarget == (overlap + 1))) |
                           ((mutNumTarget == overlap) & (mutNumSource == (overlap + 1)))]
  
  # selfSteps captures self comparisons (genotype==target)
  selfSteps = allentries[genotype == target]
  
  # Special-case: treat single mutation genotypes that match reference as edges to "reference"
  referenceEdges <- selfSteps[(`genotype` == `target`) &
                                mutNumTarget == 1]
  referenceEdges$target <- "reference"
  referenceEdges$mutNumTarget <- 0
  
  # Combine edges and export CSVs for downstream visualization/analysis
  singleSteps <- rbindlist(list(singleSteps, selfSteps, referenceEdges))
  write.csv(
    quote = F,
    row.names = F,
    singleSteps,
    file = paste(NAME, "_epistaticNetwork.csv", sep = "")
  )
  write.csv(
    quote = F,
    row.names = F,
    allentries,
    file = paste(NAME, "_genotypeNetwork.csv", sep = "")
  )
  return(binaryMatrix)
}

# runGenotypeAnalysis:
# Controls the full flow of analysis:
# - accepts an input table, a NAME prefix, a consensus sequence (string),
# - options to build a network and/or plot haplotypes.
# Steps:
# 1. run haploanalysis to extract per-mutation, per-cell counts and frequencies
# 2. annotate each variant (token) relative to the consensus reference
# 3. merge annotation information back with haplotype counts and compute per-genotype/haplotype frequencies
# 4. optionally build genotype networks and plot haplotype summaries
# Returns: haplocountsAnnot - a data.table with mutation annotation and per-genotype/haplotype frequency metrics.
runGenotypeAnalysis <- function(input, NAME = "", consensusSequence,network=F,plothaps=F) {
  consensus <- toupper(consensusSequence)
  # Compute haplotype counts / per-mutation frequencies across cells
  haplocounts = haploanalysis(input)
  
  # Break out genotype tokens for downstream annotation
  haplocountsVariants<-haplocounts[,(strsplit(genotype,"_")[[1]]),by=genotype]
  print("Haplotype analysis complete. Annotating mutations...")
  
  # Collect unique variant tokens (exclude "reference" sentinel)
  variants=unique(haplocounts[,strsplit(genotype,"_")[[1]],by=genotype]$V1)
  variants<-variants[variants!='reference']
  
  # Annotate each variant token using annotateMutation() and the consensus sequence
  MutationAnnotations<-rbindlist(map(.x = variants,annotateMutation,.progress = T,consensus))
  
  # Merge annotation (position, ref AA, mutant AA, etc.) back into haplotype counts
  haplocountsAnnot<-merge.data.table(haplocounts,MutationAnnotations,all.x = T,by = c("pos","base"))
  
  # Derive human-friendly names for genotypes and compute genotype/haplotype frequencies per total cells
  haplocountsAnnot[,genotypeName:=paste(unique(subName),collapse = "_"),by=genotype]
  haplocountsAnnot[,genoFreq:=length(unique(CBC_ID))/total,by=genotypeName]
  haplocountsAnnot[,haploFreq:=length(unique(CBC_ID))/total,by=genotype]
  
  # Optionally compute genotype similarity matrices and add basic metadata
  if(network){#if `network` option is TRUE.
    binaryMatrix = hapNetworkGen(haplocounts, NAME)
    
    similarityDist <- dist(binaryMatrix, method = "manhattan")
    DF = data.frame(cbind(target = rownames(as.matrix(similarityDist))), as.matrix(similarityDist))
    DF = merge(DF,
               input[, 1:3, ],
               by.x = "target",
               by.y = "CBC_ID",
               all.x = T)
  }

  # Optionally plot haplotype mutation frequency and per-cell mutation patterns to a multi-panel PDF
  if(plothaps==T){
    pdf(paste(NAME,".pdf",sep = ""), width = 8, height = 10)
    plot(
      cowplot::plot_grid(
        align = "v",
        nrow = 2,
        ggplot(haplocounts) + ylab("Mutation Frequency") +
          scale_color_brewer("Base", palette = "Paired", direction = 1) +
          xlab("Position") +
          geom_point(aes(
            as.integer(pos), freq, col = factor(base, levels = c("A", "T", "C", "G", "Y", "R", "S", "W", "-"))
          )),
        
        ggplot(haplocounts) +
          geom_point(aes(
            as.integer(pos),
            reorder(CBC_ID, pos, max),
            col = factor(base, levels = c("A", "T", "C", "G", "Y", "R", "S", "W", "-"))
          )) +
          theme(axis.text.y = element_text(size = 2)) +
          xlab("Position") +
          ylab("Cell") +
          scale_color_brewer("Base", palette = "Paired", direction = 1)
      )
    )
    dev.off()
  }
  
  # Persist the annotated haplotype table for downstream use
  fwrite(haplocountsAnnot,file = paste(NAME,"_annot_v3.csv",sep = ""))
  return(haplocountsAnnot)
}

# readReference:
# Reads an entire reference consensus sequence from a file into a single string.
# Expects a simple FASTA-like or raw sequence file; this function does not strip headers,
# so the reference file should contain only the sequence text (or the user should supply a pure sequence file).
readReference<-function(fileName){
  reference<-readChar(fileName, file.info(fileName)$size)
  return(reference)
}

# codon:
# Translate a DNA sequence into amino acids and extract codon information at a specific nucleotide position.
# Inputs:
# - SEQ: full nucleotide sequence (character)
# - pos: 1-based nucleotide position to examine
# - base: base being substituted (not used for translation itself, provided for context)
# Returns a data.table with:
# - codon: the 3-nt codon sequence containing pos (or "WT" if pos is NA)
# - resPos: amino-acid residue index (1-based)
# - AA: translated amino-acid letter at that residue (or "WT"/"X")
codon <- function(SEQ, pos, base) {
  if(is.na(pos)){
    return(data.table(codon = as.character("WT"), resPos=0, AA="WT"))
  }
  else{
    # Translate the sequence to AA; handle fuzzy codons where possible
    AA = try(translate(no.init.codon = T,if.fuzzy.codon = c("solve", "X"),DNAString(SEQ)))
    
    # Compute which amino-acid residue contains the nucleotide position
    AApos = function(pos) (((pos-1)-(pos-1)%%3)/3)+1
    OS=pos%%3
    if (is.na(as.character(AA))){
      returnAA<-"X"
    }
    else{
      returnAA=as.character(AA[AApos(pos)])}
    # Extract the three-nucleotide codon that contains pos.
    codonSeq = substring(as.character(SEQ),first = ifelse(OS==0,((pos - 3)+1),((pos - OS)+1)),last = ifelse(OS==0,pos,((pos - OS) + 3)))
    return(data.table(codon = as.character(codonSeq), resPos=AApos(pos), AA=returnAA))
  }
}

# annotateMutation:
# Input:
# - variant: token like "123A" or the literal string "reference"
# - SEQ: reference consensus sequence string
# Purpose:
# - For a nucleotide substitution token, compute:
#   - pos (numeric), base (mutant nt), reference codon/AA, mutant codon/AA
#   - subName: short substitution name like "K34R" (RefAA-resPos-MutAA)
#   - subClass: "Syn" if synonymous, "Non-Syn" if non-synonymous, "X" for ambiguous
# Returns: data.table with fields pos, base, ref (codon info), mut (codon info), subName, subClass
annotateMutation <- function(variant, SEQ){
  if(variant=="reference"){
    SubTable=data.table(pos=0, base=NA,ref = NA, mut = NA, subName="ref",subClass="WT")
    return(SubTable)
  }
  stringlist=strsplit(variant,"")[[1]]
  NT=stringlist[length(stringlist)]
  pos=as.integer(paste(stringlist[-length(stringlist)],collapse = ""))
  nt = toupper(NT)
  # compute reference codon and AA at pos
  ref = codon(SEQ, pos)
  # create mutated sequence (single-nucleotide substitution at pos)
  stringlist=strsplit(SEQ,"")[[1]]
  newSeq=paste(c(stringlist[1:(pos-1)],toupper(nt),stringlist[(pos+1):length(stringlist)]),collapse = "")
  # compute mutant codon and AA at pos
  mut = codon(newSeq, pos, nt)
  if (!is.na(mut$AA)){
    subClass=ifelse(ref$AA==mut$AA,"Syn","Non-Syn")
  }
  else{
    subClass="X"
  }
  # subName uses RefAA + residue position + MutAA, e.g., "K34R"
  SubTable=data.table(pos=pos, base=nt,ref = ref, mut = mut, subName=paste(sep="",ref$AA,mut$resPos,mut$AA),subClass=subClass)
  return(SubTable)
}

# extractCoverageStats:
# Parse description field from a consensus_table and extract coverage / length statistics.
# Assumes description contains tokens like "coverage:NN length:MM"
# Returns the table with new columns barcode, ref, cov, length, consthreshold (string-split results).
extractCoverageStats<-function(consensus_table){
  consensus_table[,c("barcode","ref","cov","length","consthreshold"):=tstrsplit(description," ")]
  consensus_table[,cov:=as.numeric(strsplit(cov,"coverage:")[[1]][2]),by=description]
  consensus_table[,length:=as.numeric(strsplit(cov,"length:")[[1]][2]),by=description]
  return(consensus_table)
}

# readAndAnnotData:
# Convenience wrapper to read reference and data files and run the analysis pipeline.
# Inputs:
# - refFileName: path to reference sequence file (string)
# - dataFile: path to filtered consensus CSV file of haplotypes
# - NAME: output file prefix
# - network, plothaps: booleans to enable network generation and plotting respectively
# Returns: result of runGenotypeAnalysis (annotated haplotype table)
readAndAnnotData<-function(refFileName,dataFile,NAME,network=F,plothaps=F){
  # READ IN DATA
  reference<-readReference(refFileName)
  dataTable<-fread(dataFile,header = T)
  return(runGenotypeAnalysis(dataTable,NAME,reference,network,plothaps))
}

# Main function / script entry point
# The script expects command-line arguments to be available in `allArgs`:
# allArgs[[1]] : reference sequence file
# allArgs[[2]] : filtered consensus CSV file
# allArgs[[3]] : NAME/prefix for outputs
#
# Example (shell):
# Rscript Consensus_Annotation_v3.R reference.txt filtConsensus.csv output_prefix
#
# Note: In some execution contexts `commandArgs(TRUE)` is used to get arguments into `allArgs`.
# Ensure `allArgs` is defined when invoking this script.
## Determined Consensus sequence.
ref <- allArgs[[1]]
## Filtered Consensus CSV from ConsensusTools. 
filteredConsensus=allArgs[[2]]
## Run Annotation Script, with network determination.
annot=readAndAnnotData(ref,filteredConsensus,NAME=allArgs[[3]],network=T,plothaps = T)
