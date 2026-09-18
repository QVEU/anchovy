# tests/annotate_analysis.R
#
# Analysis-only copy of Consensus_Annotation.R, used to FREEZE a golden reference
# for the Python port of the annotation stage.
#
# What's changed vs the original, and why it's still a faithful golden source:
#   - No package auto-install / no hardcoded loc_lib placeholder (just library()).
#   - Only the packages the ANALYSIS path uses are loaded (data.table, Biostrings,
#     DECIPHER, purrr) -- the ggplot/cowplot/ggrepel stack is plotting-only.
#   - plothaps = FALSE (we're dropping plotting); network = TRUE.
# The analysis FUNCTIONS (haploanalysis, annotateMutation, codon, hapNetworkGen,
# runGenotypeAnalysis) are copied verbatim from the original, so the CSV outputs
# are identical to what the original produces.
#
# Usage:
#   Rscript tests/annotate_analysis.R <reference.txt> <filtConsensus.csv> <out_prefix>

suppressMessages({
  library(data.table)
  library(Biostrings)
  library(DECIPHER)
  library(purrr)
})

allArgs <- commandArgs(trailingOnly = TRUE)

# ---- analysis functions (verbatim from Consensus_Annotation.R) --------------

haploanalysis <- function(cons) {
  DEPTH = length(unique(cons$CBC_ID))
  cons[genotype=="",WT:=T,by=CBC_ID]
  referenceHaplotypes<-unique(cons[WT==T,CBC_ID])
  unrolled = cons[, c(strsplit(genotype, "_")), by = CBC_ID]
  unrolledgeno = cons[, genotype, by = CBC_ID]
  unrolled <- merge(unrolled, unrolledgeno)
  colnames(unrolled) <- c("CBC_ID", "mutants", "genotype")
  bases = unrolled[, toupper(strsplit(mutants, "")[[1]][length(strsplit(mutants, "")[[1]])]), by = mutants]
  positions = unrolled[, as.integer(paste(collapse = "", strsplit(mutants, "")[[1]][1:(length(strsplit(mutants, "")[[1]]) - 1)])), by = mutants]
  colnames(positions) <- c("mutants", "pos")
  colnames(bases) <- c("mutants", "base")
  mutantKey = merge.data.table(positions, bases)
  mutantKey<-mutantKey[base!="-"]
  mutantTable = merge.data.table(mutantKey, unrolled, by = "mutants")
  mutantTable[, BCMutCount := length(mutants), by = "CBC_ID"]
  mutantTable <- mutantTable[BCMutCount < 200]
  refTable<-data.table(mutants="",pos=NA,base=NA,CBC_ID=referenceHaplotypes,genotype="reference",BCMutCount=0)
  mutantTable<-rbindlist(fill = T,use.names = T,list(mutantTable,refTable))
  mutantTable[, total := DEPTH]
  mutantTable[, count := nrow(.SD), by = mutants]
  mutantTable[, freq := count / total]
  return(mutantTable)
}

codon <- function(SEQ, pos, base) {
  if(is.na(pos)){
    return(data.table(codon = as.character("WT"), resPos=0, AA="WT"))
  } else {
    AA = try(translate(no.init.codon = T,if.fuzzy.codon = c("solve", "X"),DNAString(SEQ)))
    AApos = function(pos) (((pos-1)-(pos-1)%%3)/3)+1
    OS=pos%%3
    if (is.na(as.character(AA))){
      returnAA<-"X"
    } else {
      returnAA=as.character(AA[AApos(pos)])
    }
    codonSeq = substring(as.character(SEQ),first = ifelse(OS==0,((pos - 3)+1),((pos - OS)+1)),last = ifelse(OS==0,pos,((pos - OS) + 3)))
    return(data.table(codon = as.character(codonSeq), resPos=AApos(pos), AA=returnAA))
  }
}

annotateMutation <- function(variant, SEQ){
  if(variant=="reference"){
    SubTable=data.table(pos=0, base=NA,ref = NA, mut = NA, subName="ref",subClass="WT")
    return(SubTable)
  }
  stringlist=strsplit(variant,"")[[1]]
  NT=stringlist[length(stringlist)]
  pos=as.integer(paste(stringlist[-length(stringlist)],collapse = ""))
  nt = toupper(NT)
  ref = codon(SEQ, pos)
  stringlist=strsplit(SEQ,"")[[1]]
  newSeq=paste(c(stringlist[1:(pos-1)],toupper(nt),stringlist[(pos+1):length(stringlist)]),collapse = "")
  mut = codon(newSeq, pos, nt)
  if (!is.na(mut$AA)){
    subClass=ifelse(ref$AA==mut$AA,"Syn","Non-Syn")
  } else {
    subClass="X"
  }
  SubTable=data.table(pos=pos, base=nt,ref = ref, mut = mut, subName=paste(sep="",ref$AA,mut$resPos,mut$AA),subClass=subClass)
  return(SubTable)
}

hapNetworkGen <- function(haplocounts, NAME) {
  binaryMatrix = dcast.data.table(haplocounts, genotype ~ mutants, fill = 0,
    value.var = "freq", fun.aggregate = function(X) as.integer(any(X > 0)))
  countMatrix = dcast.data.table(haplocounts, genotype ~ mutants, fill = 0, value.var = "count")
  counts = data.table(countMatrix$genotype, rowSums(countMatrix[, -1]) / rowSums(binaryMatrix[, -1]))
  colnames(counts) <- c("genotype", "count")
  allentries = NULL
  genoList <- as.factor(binaryMatrix$genotype)
  for (i in as.factor(binaryMatrix$genotype)) {
    for (j in genoList) {
      Source = t(binaryMatrix[genotype == i])
      Target = t(binaryMatrix[genotype == j])
      comparison <- data.frame(Target, Source)
      fcomp <- comparison[Target == 1 & Source == 1, ]
      if (nrow(fcomp) > 0) {
        entry = data.table(source = i, target = j, overlap = nrow(fcomp),
          mutNumSource = length(strsplit(i, "_")[[1]]),
          mutNumTarget = length(strsplit(j, "_")[[1]]))
        allentries <- rbindlist(list(allentries, entry))
      }
    }
    genoList <- genoList[genoList != i]
  }
  allentries <- merge.data.table(counts, allentries, by.x = "genotype", by.y = "source")
  singleSteps = allentries[((mutNumSource == overlap) & (mutNumTarget == (overlap + 1))) |
                           ((mutNumTarget == overlap) & (mutNumSource == (overlap + 1)))]
  selfSteps = allentries[genotype == target]
  referenceEdges <- selfSteps[(`genotype` == `target`) & mutNumTarget == 1]
  referenceEdges$target <- "reference"
  referenceEdges$mutNumTarget <- 0
  singleSteps <- rbindlist(list(singleSteps, selfSteps, referenceEdges))
  write.csv(quote = F, row.names = F, singleSteps, file = paste(NAME, "_epistaticNetwork.csv", sep = ""))
  write.csv(quote = F, row.names = F, allentries, file = paste(NAME, "_genotypeNetwork.csv", sep = ""))
  return(binaryMatrix)
}

runGenotypeAnalysis <- function(input, NAME = "", consensusSequence, network=F, plothaps=F) {
  consensus <- toupper(consensusSequence)
  haplocounts = haploanalysis(input)
  variants=unique(haplocounts[,strsplit(genotype,"_")[[1]],by=genotype]$V1)
  variants<-variants[variants!='reference']
  MutationAnnotations<-rbindlist(map(.x = variants,annotateMutation,consensus))
  haplocountsAnnot<-merge.data.table(haplocounts,MutationAnnotations,all.x = T,by = c("pos","base"))
  haplocountsAnnot[,genotypeName:=paste(unique(subName),collapse = "_"),by=genotype]
  haplocountsAnnot[,genoFreq:=length(unique(CBC_ID))/total,by=genotypeName]
  haplocountsAnnot[,haploFreq:=length(unique(CBC_ID))/total,by=genotype]
  if(network){
    binaryMatrix = hapNetworkGen(haplocounts, NAME)
  }
  fwrite(haplocountsAnnot,file = paste(NAME,"_annot_v3.csv",sep = ""))
  return(haplocountsAnnot)
}

readReference<-function(fileName){
  readChar(fileName, file.info(fileName)$size)
}

# ---- entry point ------------------------------------------------------------
ref <- readReference(allArgs[[1]])
filteredConsensus <- fread(allArgs[[2]], header = TRUE)
invisible(runGenotypeAnalysis(filteredConsensus, NAME = allArgs[[3]],
                              consensusSequence = ref, network = TRUE, plothaps = FALSE))
cat("Analysis complete. Wrote", paste0(allArgs[[3]], "_annot_v3.csv"),
    "and network CSVs.\n")
