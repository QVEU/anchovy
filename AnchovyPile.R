library(data.table)
inDir="/Volumes/lvd_qve/Sequencing_Data/QVEU_Seq_0010_Minion_ndas10xlib2/no_sample/20220727_2153_MC-113212_FAT27957_6eb71326/fastq_pass/"
for (barcode in list.dirs(inDir,full.names = T)){
  print(barcode)
  pileList<-list.files(barcode,pattern = ".pile")
  pileList<-pileList[pileList!="merge_sort_pile.pile"] 
  allPU<-NULL
  for (i in (pileList)){
    print(i)
    inputFile=fread(paste(barcode,i,sep = "/"))
    if (length(inputFile)>0){
      #print(inputFile)
      inputFile[,ID:=i]
      allPU=rbindlist(list(allPU,inputFile))
      print(allPU)  
    }
  }
  print(allPU)
  if (is.null(allPU)==F){
    BC_plots<-ggplot(allPU)+geom_line(aes(V2,V4,col=ID,group=ID),show.legend = F)+facet_wrap(~ID)+ggtitle(barcode)
    ggsave(filename = paste(barcode,".pdf",sep=""),plot = BC_plots)  
  }
  
}
