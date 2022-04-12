library(seqinr)
library("parallel")
library("foreach")
library(doSNOW)
library(progress)
library(Biostrings)
library(msa)
library("ggplot2")
library(nnet)
library(ggpubr)
library(writexl)
library(seqinr)

 cores<-10

results <-"/Noise/"   #Docker
temp <-paste(results,"rawnoise/",sep = "")   #Docker

bamfiles<-list.files(results, pattern = ".bam$", full.names = TRUE, recursive = TRUE)

if(length(bamfiles)>0){
  
  samp<-c(1:length(bamfiles))
  
  pb <- progress_bar$new(
    format = "Index: :samp.pb [:bar] :elapsed | eta: :eta",
    total = length(bamfiles),    # 100 
    width = 60)
  
  
  progress <- function(n){
    pb$tick(tokens = list(samp.pb = samp[n]))
  } 
  
  opts <- list(progress = progress)
  
  if(!dir.exists(temp)) dir.create(temp)
  cluster.cores<-makeCluster(cores)
  registerDoSNOW(cluster.cores)
  
  out.par<-foreach(i=1:length(bamfiles), .verbose=FALSE, .options.snow = opts) %dopar%{
    
    try(system(paste("cd /Noise && FINex -f ",bamfiles[i], " > ",
                     gsub(".*/",temp,gsub("\\.bam", "_NoisExtractorResult.tsv",bamfiles[i])), sep = "")))
    
  }
  stopCluster(cluster.cores)
  
  results.files<-list.files(temp, pattern = "_NoisExtractorResult\\.tsv$", full.names = TRUE)


 
  co<-0.1
  pb<-txtProgressBar(min = 0, max = length(results.files),initial = 0) #01072021 Nacho Garcia
  try(rm(extended.out))
  out.plots<-list()
  for (i in 1:length(results.files)) {
    
    try(rm(dummy), silent = TRUE)
    try(dummy<-read.csv(results.files[i],sep = "\t", header = FALSE), silent = TRUE)
   
    if(exists("dummy")){
      dummy<-dummy[,c(1:5)]
      colnames(dummy)<-c("Base","Noise","Reads","S1","S2")
      
      outlier.co<- max(median(dummy$Noise)+10*sd(dummy$Noise),0.10)
    
      genome.position<-as.data.frame(c(1:29903))
      colnames(genome.position)<-"Base"
      dummy<-merge(genome.position, dummy, by="Base", all=TRUE)
      dummy$Seq1[which(is.na(dummy$Reads))]<-"N"
      dummy$Reads[which(is.na(dummy$Reads))]<-0
      
      dummy$Outlier<-"NO"
      dummy$Outlier[dummy$Noise>=outlier.co]<-"YES"
      Overall.Noise<-sum(dummy$Noise[dummy$Outlier=="YES"])
      if(length(which(dummy$Noise>co))>5){
        
        if(!dir.exists(paste(results,"fasta/",sep = ""))) try(dir.create(paste(results,"fasta/",sep = "")))
        
        genome.position<-as.data.frame(c(1:29903))
        colnames(genome.position)<-"Base"
        dummy<-merge(dummy,genome.position,by="Base",all.y=TRUE)
        dummy$S1[which(is.na(dummy$S1))]<-"N"
        dummy$S2[which(is.na(dummy$S2))]<-"N"
        dummy$S2[which(dummy$Reads<10)]<-"N"
        dummy$S1[which(dummy$Reads<10)]<-"N"
        
        dummy$S2[which(dummy$Noise<co)]<-dummy$S1[which(dummy$Noise<co)]
        if(length(which(dummy$S1!=dummy$S2))>5 & length(which(dummy$S1=="N"))<1000){
          
          write.fasta(paste(dummy$S1[which(dummy$S1 %in% c("A","T","C","G","N"))], collapse = ""),
                      file.out = gsub("_NoisExtractorResult.tsv","_S1.fa", gsub(".*/",paste(results,"fasta/",sep = ""),results.files[i])), 
                      names = gsub("_NoisExtractorResult.tsv","_S1.fa", gsub(".*/","",results.files[i])))
          
          write.fasta(paste(dummy$S2[which(dummy$S2 %in% c("A","T","C","G","N"))], collapse = ""),
                      file.out = gsub("_NoisExtractorResult.tsv","_S2.fa", gsub(".*/",paste(results,"fasta/",sep = ""),results.files[i])), 
                      names = gsub("_NoisExtractorResult.tsv","_S2.fa", gsub(".*/","",results.files[i])))   
        }
        
      }
      
      names<-gsub("\\.sorted.*","",gsub("_S[0-9].*","",gsub("R[0-9].*","",gsub(".*/","",results.files[i]))))
      out.plots[[length(out.plots)+1]]<-ggplot(dummy)+
        geom_line(aes(Base, Noise))+
        geom_point(data=subset(dummy, Outlier=="YES"),aes(Base, Noise),col="red", alpha=0.3)+
        geom_point(data=subset(dummy, Reads<10),aes(Base, 0),col="blue", alpha=0.1)+
        ylim(0,1)+
        theme_minimal()+
        ggtitle(paste(names, " /Noise:", round(Overall.Noise, 2), sep=""))
    }
  }
  date<-gsub("-","",Sys.Date())
  if(length(list.files("/Noise/fasta/"))>0){
    system("cat /Noise/fasta/*.fa > /Noise/fasta/Coinfections_total.fa")
    system(paste("pangorunner.sh", "/Noise/fasta/Coinfections_total.fa"))
    df<- read.csv("/Noise/fasta/Coinfections_total.fa_pango.csv")
    
    df$Sample<-"S1"
    df$Sample[grep("_S2",df$taxon)]<-"S2"
    df$ID<-gsub( "_S.\\.fa","",df$taxon)
    
    
    S1.df<-df[which(df$Sample=="S1"),c("lineage","conflict","note","ID")]
    S2.df<-df[which(df$Sample=="S2"),c("lineage","conflict","note","ID")]
    
    colnames(S1.df)[-4]<-paste("S1_",colnames(S1.df)[-4], sep = "")
    colnames(S2.df)[-4]<-paste("S2_",colnames(S2.df)[-4], sep = "")
    
    total<-merge(S1.df, S2.df, by="ID")
    total<-total[,c("ID", "S1_lineage", "S2_lineage")]
    colnames(total)<-c("Sample", "Major", "Minor")
    
    if(length(which(total$Major==total$Minor))>0){
      total<-total[-which(total$Major==total$Minor),]
    }
    
    if(nrow(total)>0) write.csv(total, paste("/Noise/Coinfection_Results", date,".csv",sep=""), row.names = FALSE)
    if(nrow(total)==0) write.table("No Coinfections Found!",paste("/Noise/Coinfection_Results", date,".csv",sep=""), row.names = FALSE, quote = FALSE, col.names = FALSE)
    try(system("rm -rf /Noise/fasta"))
  }

  
  if(length(out.plots)<=40 ){
    if(length(out.plots)>0){
      ggarrange(plotlist =  out.plots[1:length(out.plots)], ncol = 4, nrow = 10)
      ggsave(paste(results,"ResultsNoisExtractor_",date,"_",".pdf", sep=""), width = 420, height = 600, units = "mm") #A4
    }  
  }else{
    plotting<-TRUE
    start<-1
    end<-40
    counter<-0
    
    while(plotting){
      if(end==length(out.plots)) plotting<-FALSE
      ggarrange(plotlist =  out.plots[start:end], ncol = 4, nrow = 10)
      start<-end+1
      end<-end+40
      if(end>=length(out.plots)) end<-length(out.plots)
      counter<-counter+1
      
      ggsave(paste(results,"ResultsNoisExtractor_",date,"_","_",counter,".pdf", sep=""), width = 420, height = 600, units = "mm") #A4
      
    }
  }
 
 library("pdftools")
  try(system("rm -rf /Noise/rawnoise"))
 pdf.list<-list.files(results, full.names = TRUE, pattern = ".*NoisExtractor.*\\.pdf")
 if(length(pdf.list)>1){
 pdf_combine(pdf.list, output = gsub("_.\\.pdf","_Merged.pdf",pdf.list[1]))
 file.remove(pdf.list)}
}else{
    print("No bam files found in the /Noise/ folder")
  }
  
