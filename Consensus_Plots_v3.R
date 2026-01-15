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
  
  fwrite(EV71_P1_annot,file = "EVA71_P1_annot_v3.csv")
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

ref71 <- "/Volumes/dolanpt/EV71_P1_consensus.txt"

EV71_6H_P1="/Volumes/LVD_QVE/Projects/PacBio_virus-inclusive_EVA71Passage/filtConsensus_EVA71ORF_6h_P1.csv"
EV71_P1_annot=readAndAnnotData(ref71,EV71_6H_P1,"A71_P1_6hpi")

EV71_6H_P3="/Volumes/lvd_qve/Projects/PacBio_virus-inclusive_EVA71Passage/filtConsensus_EVA71ORF_6h_P3.csv"
EV71_P3_annot<-readAndAnnotData(ref71,EV71_6H_P3,"A71_P3_6hpi")

EV71_6H_P5="/Volumes/lvd_qve/Projects/PacBio_virus-inclusive_EVA71Passage/filtConsensus_EVA71ORF_6h_P5.csv"
EV71_P5_annot<-readAndAnnotData(ref71, EV71_6H_P5,"A71_P5_6hpi")

EV71_P5_annot[,list(subName),by=c("genotypeName")]

#Passage Data organization
fwrite(EV71_P1_annot,file = "EVA71_P1_annot_v2.csv")
fwrite(EV71_P3_annot,file = "EVA71_P3_annot_v2.csv")
fwrite(EV71_P5_annot,file = "EVA71_P5_annot_v2.csv")

EV71_P1_annot<-fread(file = "EVA71_P1_annot_v2.csv")
EV71_P3_annot<-fread(file = "EVA71_P3_annot_v2.csv")
EV71_P5_annot<-fread(file = "EVA71_P5_annot_v2.csv")

EV71_P1_annot[,passage:="P1"]
EV71_P3_annot[,passage:="P3"]
EV71_P5_annot[,passage:="P5"]

EV71_P1_annot[,genoFreq:=length(unique(CBC_ID))/total,by=genotypeName]
EV71_P3_annot[,genoFreq:=length(unique(CBC_ID))/total,by=genotypeName]
EV71_P5_annot[,genoFreq:=length(unique(CBC_ID))/total,by=genotypeName]

EV71passages<-rbindlist(list(EV71_P1_annot,EV71_P3_annot,EV71_P5_annot))


#
write.csv(file = paste("genotypeKey.csv",sep=""),unique(EV71passages[,c("genotype","genotypeName")]),quote = F,row.names = F)

countsTable<-dcast.data.table(EV71passages,formula = genotypeName~passage,value.var = "CBC_ID",fun.aggregate = function(C) (length(unique(C))))
countsTable<-melt.data.table(countsTable,id.vars = 1,variable.name = "passage",value.name = "genotypeCount")
countsTablewRef<-rbindlist(list(countsTable,data.table(c('reference','reference','reference'),
                                                       c('P1','P3','P5'),
                                                       c(94,32,19))))

#Mueller Diagram
ggplot(countsTable)+
  geom_area(col="black",lwd=0.2,position='fill',stat = "identity",aes(passage,genotypeCount,group=genotypeName,fill=genotypeName),show.legend = F)+
  scale_fill_viridis_d(sample(length(countsTable$genotypeCount)))+theme_bw()

ggplot(countsTable)+
  geom_area(col="black",lwd=0.1,position='fill',stat = "identity",aes(passage,genotypeCount,group=reorder(genotypeName,-genotypeCount,mean),fill=genotypeName),show.legend = F)+
  scale_fill_viridis_d(sample(length(countsTable$genotypeCount)))+theme_bw()

THEME<-theme(axis.text.y = element_text(size=0.4))
THEME2<-theme(axis.text.y = element_text(size=7))

EV_Features<-fread("/Volumes/lvd_qve/QVEU_Code/bioinfo/EV71_4643_Features.csv")
EV68_Features<-fread("/Volumes/lvd_qve/QVEU_Code/bioinfo/EVD68_18947_MO-14_Features.csv")
CVB_Features<-fread("/Volumes/lvd_qve/QVEU_Code/bioinfo/EVD68_18947_MO-14_Features.csv")

pdf("EV71_passage_annotated.pdf",width=18,height=8)
cowplot::plot_grid(align = "hv",ncol = 3,nrow = 2,byrow = F,plotlist = 
                     list(
                       #Plot P1
                       ggplot(EV71_P1_annot)+theme_half_open()+
                         ylab("SNV Frequency")+
                         xlab("nt position (ORF)")+
                         geom_rect(data=EV_Features,aes(xmin=Start-746,xmax=End-746,ymin=0,ymax=1),alpha=rep(c(0,0.3,0,0.3,0,0.3,0,0.3,0,0.3,0)))+
                         geom_text(data=EV_Features,aes(x=(Start-746+End-746)/2,y=0.9,label=Product),size=3,angle=90,align="v")+
                         ggtitle("P1")+
                         geom_point(size=1,aes(pos,freq,col=subClass,))+
                         scale_color_manual("Sub. Class",breaks = c("Non-Syn","Syn","X"),values = c("darkred","blue","darkgrey"))+
                         THEME2,
                       
                        ggplot(EV71_P1_annot)+theme_half_open()+
                          ylab("Cell Haplotypes")+
                         xlab("nt position (ORF)")+
                          geom_point(size=0.3,aes(pos,reorder(CBC_ID,genotype),col=subClass))+
                          #geom_text(hjust=0,nudge_x=20,nudge_y=-20,size=3,aes(0,length(unique(CBC_ID)),label=paste(collapse = "","N:", P1Total,"\n",as.character(round(100-100*(length(unique(CBC_ID))/P1Total),digits = 1)),"% Master Sequence")))+
                          scale_color_manual("Sub. Class",breaks = c("Non-Syn","Syn","X"),values = c("darkred","blue","darkgrey"))+
                          THEME,
                        
                        #Plot P3
                        ggplot(EV71_P3_annot)+theme_half_open()+
                         ylab("SNV Frequency")+
                          geom_rect(data=EV_Features,aes(xmin=Start-746,xmax=End-746,ymin=0,ymax=1),alpha=rep(c(0,0.3,0,0.3,0,0.3,0,0.3,0,0.3,0)))+
                          geom_text(data=EV_Features,aes(x=(Start-746+End-746)/2,y=0.9,label=Product),size=3,angle=90)+
                          ggtitle("P3")+
                         xlab("nt position (ORF)")+
                          geom_point(size=1,aes(pos,col=subClass,freq))+
                          scale_color_manual("Sub. Class",breaks = c("Non-Syn","Syn","X"),values = c("darkred","blue","darkgrey"))+
                         THEME2,
                        
                        ggplot(EV71_P3_annot)+theme_half_open()+
                         ylab("Cell Haplotypes")+
                         xlab("nt position (ORF)")+
                          geom_point(size=.3,aes(pos,col=subClass,reorder(CBC_ID,genotype)))+
                          #geom_text(hjust=0,nudge_x=20,nudge_y=-20,size=3,aes(0,length(unique(CBC_ID)),label=paste(collapse = "","N:", P3Total,"\n",as.character(round(100-100*(length(unique(CBC_ID))/P3Total),digits = 1)),"% Master Sequence")))+
                          scale_color_manual("Sub. Class",breaks = c("Non-Syn","Syn","X"),values = c("darkred","blue","darkgrey"))+
                          THEME,
                        
                        #Plot P5
                       ggplot(EV71_P5_annot)+
                         theme_half_open()+
                         ggtitle("P5")+
                         xlab("nt position (ORF)")+
                         geom_rect(data=EV_Features,aes(xmin=Start-746,xmax=End-746,ymin=0,ymax=1),alpha=rep(c(0,0.3,0,0.3,0,0.3,0,0.3,0,0.3,0)))+
                         geom_text(data=EV_Features,aes(x=(Start-746+End-746)/2,y=0.9,label=Product),size=3,angle=90)+
                         geom_point(size=1,aes(pos,freq,col=subClass))+
                         ylab("SNV Frequency")+
                         scale_color_manual("Sub. Class",breaks = c("Non-Syn","Syn","X"),values = c("darkred","blue","darkgrey"))+
                         THEME2,
                      
                        ggplot(EV71_P5_annot)+
                         theme_half_open()+
                         xlab("nt position (ORF)")+
                         ylab("Cell Haplotypes")+
                          geom_point(size=0.3,aes(pos,reorder(CBC_ID,pos,mean),col=subClass))+
                          #geom_text(hjust=0,nudge_x=20,nudge_y=-40,size=3,aes(0,length(unique(CBC_ID)),label=paste(collapse = "","N:", P5Total,"\n",as.character(round(100-100*(length(unique(CBC_ID))/P5Total),digits = 1)),"% Master Sequence")))+
                          scale_color_manual("Sub. Class",breaks = c("Non-Syn","Syn","X"),values = c("darkred","blue","darkgrey"),na.value = "black")+
                          THEME
                     ))

dev.off()

#Passage Data organization

refCVB="/Volumes/dolanpt/CVB3_consensus.txt"
CVB="/Volumes/lvd_qve/Projects/PacBio_virus-inclusive/filtConsensus_CVB3ORF.csv"
CVB_annot<-readAndAnnotData(refCVB, CVB,"CV-B3")

EV71 <- "/Volumes/dolanpt/EV71_P1_consensus.txt"
EV71data="/Volumes/lvd_qve/Projects/PacBio_virus-inclusive/filtConsensus_EVA71ORF.csv"
EV71_annot<-readAndAnnotData(EV71, EV71data,"EV-A71")

EVd68 <- "/Volumes/dolanpt/EV68_consensus.txt"
EVd68data ="/Volumes/lvd_qve/Projects/PacBio_virus-inclusive/filtConsensus_EVd68ORF.csv"
EVd68_annot<-readAndAnnotData(EVd68, EVd68data,"EV-D68")

fwrite( CVB_annot,  file = "/Volumes/lvd_qve/Projects/PacBio_virus-inclusive/CVB3_annot.csv")
fwrite( EV71_annot,  file = "/Volumes/lvd_qve/Projects/PacBio_virus-inclusive/EV71_annot.csv")
fwrite( EVd68_annot, file = "/Volumes/lvd_qve/Projects/PacBio_virus-inclusive/EVD68_annot.csv")

CVB_annot<-fread(file = "/Volumes/lvd_qve/Projects/PacBio_virus-inclusive/CVB3_annot.csv")
EV71_annot<-fread(file = "/Volumes/lvd_qve/Projects/PacBio_virus-inclusive/CVB3_annot.csv")
EVd68_annot<-fread(file = "/Volumes/lvd_qve/Projects/PacBio_virus-inclusive/CVB3_annot.csv")

CVB_annot[,passage:="CVB3"]
EV71_annot[,passage:="EVA71"]
EVd68_annot[,passage:="EVD68"]

threeViruses<-rbindlist(list(CVB_annot,EV71_annot,EVd68_annot))

fwrite(threeViruses,file = "/Volumes/lvd_qve/Projects/PacBio_virus-inclusive/EV_3virus_annotated.csv")

pdf("/Volumes/lvd_qve/Projects/PacBio_virus-inclusive/EV_3virus_annotated.pdf",width=18,height=8)

cowplot::plot_grid(byrow = F,align = "hv",ncol = 3,plotlist = 
                     list(
                       #EV71
                       ggplot(EV71_annot)+
                         theme_half_open()+
                         ylab("SNV Frequency")+
                         geom_rect(data=EV_Features,aes(xmin=Start-746,xmax=End-746,ymin=0,ymax=0.5),alpha=rep(c(0,0.3,0,0.3,0,0.3,0,0.3,0,0.3,0)))+
                         geom_text(data=EV_Features,aes(x=(Start-746+End-746)/2,y=0.45,label=Product),size=3,angle=90)+
                         ggtitle("EV-A71")+
                         xlab("nt position (ORF)")+
                         geom_point(size=1,aes(pos,col=subClass,freq))+
                         scale_color_manual("Sub. Class",breaks = c("Non-Syn","Syn","X"),values = c("darkred","blue","darkgrey"))+
                         THEME2,
                       
                       ggplot(EV71_annot)+theme_half_open()+
                         ylab("Cell Haplotypes")+
                         xlab("nt position (ORF)")+
                         geom_point(size=.3,aes(pos,col=subClass,reorder(CBC_ID,genotype)))+
                         #geom_text(hjust=0,nudge_x=20,nudge_y=-20,size=3,aes(0,length(unique(CBC_ID)),label=paste(collapse = "","N:", 289,"\n",as.character(round(100-100*(length(unique(CBC_ID))/289),digits = 1)),"% Master Sequence")))+
                         scale_color_manual("Sub. Class",breaks = c("Non-Syn","Syn","X"),values = c("darkred","blue","darkgrey"))+
                         THEME,
                       
                       #Plot d68
                       ggplot(EVd68_annot)+
                         theme_half_open()+
                         ggtitle("EV-D68")+
                         xlab("nt position (ORF)")+
                         geom_rect(data=EV68_Features,aes(xmin=Start-697,xmax=End-697,ymin=0,ymax=0.5),alpha=rep(c(0,0.3,0,0.3,0,0.3,0,0.3,0,0.3,0)))+
                         geom_text(data=EV68_Features,aes(x=(Start-697+End-697)/2,y=0.45,label=Product),size=3,angle=90)+
                         geom_point(size=1,aes(pos,freq,col=subClass))+
                         ylab("SNV Frequency")+
                         scale_color_manual("Sub. Class",breaks = c("Non-Syn","Syn","X"),values = c("darkred","blue","darkgrey"))+
                         THEME2,
                       
                       ggplot(EVd68_annot)+
                         theme_half_open()+
                         xlab("nt position (ORF)")+
                         ylab("Cell Haplotypes")+
                         geom_point(size=0.3,aes(pos,reorder(CBC_ID,pos,mean),col=subClass))+
                         #geom_text(hjust=0,nudge_x=20,nudge_y=-40,size=3,aes(0,length(unique(CBC_ID)),label=paste(collapse = "","N:", 255,"\n",as.character(round(100-100*(length(unique(CBC_ID))/255),digits = 1)),"% Master Sequence")))+
                         scale_color_manual("Sub. Class",breaks = c("Non-Syn","Syn","X"),values = c("darkred","blue","darkgrey"),na.value = "black")+
                         THEME,
                       
                       #Plot CVB3
                       ggplot(CVB_annot)+
                         theme_half_open()+
                         ylab("SNV Frequency")+
                         xlab("nt position (ORF)")+
                         geom_rect(data=CVB_Features,aes(xmin=Start-696,xmax=End-696,ymin=0,ymax=0.5),alpha=rep(c(0,0.3,0,0.3,0,0.3,0,0.3,0,0.3,0)))+
                         geom_text(data=CVB_Features,aes(x=(Start-696+End-696)/2,y=0.45,label=Product),size=3,angle=90)+
                         ggtitle("CV-B3")+
                         geom_point(size=1,aes(pos,freq,col=subClass,))+
                         scale_color_manual("Sub. Class",breaks = c("Non-Syn","Syn","X"),values = c("darkred","blue","darkgrey"))+
                         THEME2,
                       
                       ggplot(CVB_annot)+theme_half_open()+
                         ylab("Cell Haplotypes")+
                         xlab("nt position (ORF)")+
                         geom_point(size=0.3,aes(pos,reorder(CBC_ID,genotype),col=subClass))+
                         #geom_text(hjust=0,nudge_x=20,nudge_y=-20,size=3,aes(0,length(unique(CBC_ID)),label=paste(collapse = "","N:", 404,"\n",as.character(round(100-100*(length(unique(CBC_ID))/404),digits = 1)),"% Master Sequence")))+
                         scale_color_manual("Sub. Class",breaks = c("Non-Syn","Syn","X"),values = c("darkred","blue","darkgrey"))+
                         THEME
                     ))
dev.off()

