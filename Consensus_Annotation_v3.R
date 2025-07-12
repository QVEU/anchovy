#
# Consensus_Plots_v3.R
# This generates the network 
#
library(data.table)
library(Biostrings)
library(reshape2)
library(ggplot2)
library(DECIPHER)
library(cowplot)
library(purrr)
library(ggrepel)
#FUNCTIONS

haploanalysis <- function(cons) {
  DEPTH = length(unique(cons$CBC_ID))
  print(paste("Number of Cells:", DEPTH))
  
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
  print(mutantTable)
  mutantTable[, total := DEPTH]
  mutantTable[, count := nrow(.SD), by = mutants]
  mutantTable[, freq := count / total]
  print(mutantTable)
  return(mutantTable)
}

hapNetworkGen <- function(haplocounts, NAME) {
  #print(haplocounts)
  binaryMatrix = dcast.data.table(
    haplocounts,
    genotype ~ mutants,
    fill = 0,
    value.var = "freq",
    fun.aggregate = function(X)
      ifelse(X > 0, 1, 0)
  )
  
  binCov=cov(binaryMatrix)
  print(binCov[1:10])
  
  print(princomp(binaryMatrix))
  
  countMatrix = dcast.data.table(haplocounts,#reshape the haplotype counts into 
                                 genotype ~ mutants,
                                 fill = 0,
                                 value.var = "count")
  
  counts = data.table(countMatrix$genotype,
                      rowSums(countMatrix[, -1]) / rowSums(binaryMatrix[, -1]))
  
  colnames(counts) <- c("genotype", "count")
  #print(counts)
  allentries = NULL
  genoList <- as.factor(binaryMatrix$genotype)
  for (i in as.factor(binaryMatrix$genotype)) {
    #For each genotype,
    print(i)
    for (j in genoList) {
      #Compare to each genotype
      Source = t(binaryMatrix[genotype == i])
      Target = t(binaryMatrix[genotype == j])
      comparison <- data.frame(Target, Source)
      fcomp <- comparison[Target == 1 & Source == 1, ]
      #print(fcomp)
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
    genoList <-
      genoList[genoList != i]#remove genotype in future loops.
  }
  allentries <-
    merge.data.table(counts, allentries, by.x = "genotype", by.y = "source")
  
  singleSteps = allentries[((mutNumSource == overlap) & (mutNumTarget == (overlap + 1))) |
                           ((mutNumTarget == overlap) & (mutNumSource == (overlap + 1)))]
  
  selfSteps = allentries[genotype == target]
  referenceEdges <- selfSteps[(`genotype` == `target`) &
                                mutNumTarget == 1]
  referenceEdges$target <- "Consensus"
  referenceEdges$mutNumTarget <- 0
  #print(referenceEdges)
  singleSteps <- rbindlist(list(singleSteps, selfSteps, referenceEdges))
  #print(singleSteps)
  write.csv(
    quote = F,
    row.names = F,
    singleSteps,
    file = paste("~/Research/", NAME, "_epistaticNetwork.csv", sep = "")
  )
  write.csv(
    quote = F,
    row.names = F,
    allentries,
    file = paste("~/Research/", NAME, "_genotypeNetwork.csv", sep = "")
  )
  return(binaryMatrix)
}

# runGenotypeAnalysis(): function to control all steps in analysis of 
runGenotypeAnalysis <- function(input, NAME = "", consensusSequence,network=F,plothaps=F) {
  consensus <- toupper(consensusSequence)
  #print(consensus)
  haplocounts = haploanalysis(input)
  haplocountsVariants<-haplocounts[,(strsplit(genotype,"_")[[1]]),by=genotype]
  print("Haplotype analysis complete. Annotating mutations...")
  variants=unique(haplocounts[,strsplit(genotype,"_")[[1]],by=genotype]$V1)
  variants<-variants[variants!='reference']
  MutationAnnotations<-rbindlist(map(.x = variants,annotateMutation,.progress = T,consensus))
  
  haplocountsAnnot<-merge.data.table(haplocounts,MutationAnnotations,all.x = T,by = c("pos","base"))
  haplocountsAnnot[,genotypeName:=paste(unique(subName),collapse = "_"),by=genotype]
  haplocountsAnnot[,genoFreq:=length(unique(CBC_ID))/total,by=genotypeName]
  haplocountsAnnot[,haploFreq:=length(unique(CBC_ID))/total,by=genotype]
  
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

if(plothaps==T){#if haplotype plotting is selected (plothaps=T in )
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
  
  fwrite(haplocountsAnnot,file = paste(NAME,"_annot_v3.csv",sep = ""))
  return(haplocountsAnnot)
}

#readReference: read in reference sequence for comparison to reads
readReference<-function(fileName){
  reference<-readChar(fileName, file.info(fileName)$size)
  return(reference)
}

#codon: function to translate codons
codon <- function(SEQ, pos, base) {
  if(is.na(pos)){
    return(data.table(codon = as.character("WT"), resPos=0, AA="WT"))
  }
  else{
    AA = try(translate(no.init.codon = T,if.fuzzy.codon = c("solve", "X"),DNAString(SEQ)))
    
    AApos = function(pos) (((pos-1)-(pos-1)%%3)/3)+1
    OS=pos%%3
    if (is.na(as.character(AA))){
      returnAA<-"X"
    }
    else{
      returnAA=as.character(AA[AApos(pos)])}
    codonSeq = substring(as.character(SEQ),first = ifelse(OS==0,((pos - 3)+1),((pos - OS)+1)),last = ifelse(OS==0,pos,((pos - OS) + 3)))
    return(data.table(codon = as.character(codonSeq), resPos=AApos(pos), AA=returnAA))
  }
}

# annotateMutation: function to annotate SNV based on coding frame
annotateMutation <- function(variant, SEQ){
  if(variant=="reference"){
    SubTable=data.table(pos=0, base=NA,ref = NA, mut = NA, subName="ref",subClass="WT")
    return(SubTable)
  }
  stringlist=strsplit(variant,"")[[1]]
  NT=stringlist[length(stringlist)]
  pos=as.integer(paste(stringlist[-length(stringlist)],collapse = ""))
  nt = toupper(NT)
  #compute reference codon and AA
  ref = codon(SEQ, pos)
  #print(ref)
  stringlist=strsplit(SEQ,"")[[1]]
  newSeq=paste(c(stringlist[1:(pos-1)],toupper(nt),stringlist[(pos+1):length(stringlist)]),collapse = "")
  #compute mutant codon and AA
  mut = codon(newSeq, pos, nt)
  #print(mut)
  if (!is.na(mut$AA)){
    subClass=ifelse(ref$AA==mut$AA,"Syn","Non-Syn")
  }
  else{
    subClass="X"
  }
  SubTable=data.table(pos=pos, base=nt,ref = ref, mut = mut, subName=paste(sep="",ref$AA,mut$resPos,mut$AA),subClass=subClass)
  return(SubTable)
}

extractCoverageStats<-function(consensus_table){
  consensus_table[,c("barcode","ref","cov","length","consthreshold"):=tstrsplit(description," ")]
  consensus_table[,cov:=as.numeric(strsplit(cov,"coverage:")[[1]][2]),by=description]
  consensus_table[,length:=as.numeric(strsplit(cov,"length:")[[1]][2]),by=description]
  return(consensus_table)
}

readAndAnnotData<-function(refFileName,dataFile,NAME,network=F,plothaps=F){
#READ IN DATA
  reference<-readReference(refFileName)
  dataTable<-fread(dataFile,header = T)
  return(runGenotypeAnalysis(dataTable,NAME,reference,network,plothaps))
}


ref71 <- "/Volumes/LVD_qve/Projects/Freeman_Collab/Anchovy/Kinnex_Orgs_Consensus.txt"

EV71_Kinnex="/Volumes/LVD_qve/Projects/Freeman_Collab/Anchovy/filtConsensus_Kinnex.csv"
EV71_Kinnex_annot=readAndAnnotData(ref71,EV71_Kinnex,"/Volumes/LVD_qve/Projects/Freeman_Collab/Anchovy/Kinnex_EVA71")
