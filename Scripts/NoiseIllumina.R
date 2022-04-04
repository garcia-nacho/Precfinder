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
  
  dir.create(temp)
  cluster.cores<-makeCluster(cores)
  registerDoSNOW(cluster.cores)
  
  out.par<-foreach(i=1:length(bamfiles), .verbose=FALSE, .options.snow = opts) %dopar%{
    
    try(system(paste("FINex -f ",bamfiles[i], " > ",
                     gsub(".*/",temp,gsub("\\.bam", "_NoisExtractorResult.tsv",bamfiles[i])), sep = "")))
    
  }
  stopCluster(cluster.cores)
  
  results.files<-list.files(temp, pattern = "_NoisExtractorResult\\.tsv$", full.names = TRUE)

  out.plots<-list()
 
  summary<-as.data.frame(matrix(NA, nrow = length(results.files), ncol = 6))
  colnames(summary)<-c("Sample", "Noise", "NoisePositions", "MissingPositions", "NormNoisePositions", "OCO")
  pb<-txtProgressBar(min = 0, max = length(results.files),initial = 0) #01072021 Nacho Garcia
  try(rm(extended.out))
  out.plots<-list()
  for (i in 1:length(results.files)) {
    
    try(rm(dummy), silent = TRUE)
    try(dummy<-read.csv(results.files[i],sep = "\t", header = FALSE), silent = TRUE)
    summary$Sample[i]<-gsub(".*/","",results.files[i])
    if(exists("dummy")){
      colnames(dummy)<-c("Base","Noise","Reads","Seq1")
      
      outlier.co<- max(median(dummy$Noise)+10*sd(dummy$Noise),0.10)
      summary$OCO<-outlier.co
      genome.position<-as.data.frame(c(1:29903))
      colnames(genome.position)<-"Base"
      dummy<-merge(genome.position, dummy, by="Base", all=TRUE)
      dummy$Seq1[which(is.na(dummy$Reads))]<-"N"
      dummy$Reads[which(is.na(dummy$Reads))]<-0
      
      dummy$Outlier<-"NO"
      dummy$Outlier[dummy$Noise>=outlier.co]<-"YES"
      
      Overall.Noise<-sum(dummy$Noise[dummy$Outlier=="YES"])
      summary$Noise[i]<-Overall.Noise
      summary$NoisePositions[i]<-length(which(dummy$Outlier=="YES"))
      summary$MissingPositions[i]<-length(dummy$Base[which(dummy$Reads<5)])
      summary$NormNoisePositions[i]<- length(which(dummy$Outlier=="YES"))/(29903 - length(dummy$Base[which(dummy$Reads<5)]))
     
      if(is.infinite(summary$NormNoisePositions[i])) summary$NormNoisePositions[i]<-1 #Fix Nacho 20122021

      summary$meanNoise[i]<-mean(dummy$Noise, na.rm=TRUE)
      summary$SDNoise[i]<-sd(dummy$Noise, na.rm=TRUE)
      summary$meanCov[i]<-mean(dummy$Reads, na.rm=TRUE)
      summary$sdCov[i]<-sd(dummy$Reads, na.rm=TRUE)
      noise.hist<-hist(dummy$Base[dummy$Outlier=="YES"],breaks = seq(0,30000, by=3000), plot = FALSE)$counts
      miss.hist<-hist(dummy$Base[which(is.na(dummy$Seq1))],breaks = seq(0,30000, by=3000), plot = FALSE)$counts
      miss.noise.hist<-noise.hist*miss.hist
      miss.noise.sum<-sum(miss.noise.hist)
      total.hist<-c(miss.hist, noise.hist, miss.noise.hist, miss.noise.sum, summary$Sample[i])
      if(!exists("extended.out")){
        extended.out<-total.hist
      }else{
        extended.out<-rbind(extended.out, total.hist)
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
  #Classification
  

  date<-gsub("-","",Sys.Date())
  
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

 pdf.list<-list.files(results, full.names = TRUE, pattern = ".*NoisExtractor.*\\.pdf")
 if(length(pdf.list)>1){
 pdf_combine(pdf.list, output = gsub("_.\\.pdf","_Merged.pdf",pdf.list[1]))
 file.remove(pdf.list)}
}else{
    print("No bam files found in the /Noise/ folder")
  }
  
