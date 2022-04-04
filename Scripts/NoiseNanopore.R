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

  
# Plotting ----------------------------------------------------------------
  p.noise<-function(x){
    return(length(which(df$Noise[which(df$Base %in% stable.pos)]> x))/length(stable.pos))
  }
  p.minor<-function(x){
    return(length(which(df$FreqMinor[which(df$Base %in% stable.pos)]> x))/length(stable.pos))
  }
  
  
  results.files<-list.files(temp, pattern = "_NoisExtractorResult\\.tsv$", full.names = TRUE)

  
  out.plots<-list()
  for(i in 1:length(results.files)){

    df<-read.csv(results.files[i],sep = "\t",header = FALSE)
    
    df<-df[,c(1,2,3,6,7)]
    
    colnames(df)<-c("Base","Noise","Reads","FreqMinor","NoiseMinor")
    df$NonMR<-df$Reads*df$Noise
    df$MR<-df$Reads-df$NonMR
    
    #stable.pos<-stable.pos[which(stable.pos %in% df$Base)]
    
    
    genome.position<-as.data.frame(c(1:29903))
    colnames(genome.position)<-"Base"
    df<-merge(genome.position, df, by="Base", all=TRUE)
    df$Outlier<-"NO"
    
    lowerq = quantile(df$FreqMinor,na.rm=TRUE)[2]
    upperq = quantile(df$FreqMinor,na.rm=TRUE)[4]
    iqr = upperq - lowerq 
    co<-(iqr * 20) + upperq
    
    df$Outlier[which(df$FreqMinor>co)]<-"YES"
    
    df$Reads[which(is.na(df$Reads))]<-0

    df$NoiseNP<-0
    df$NoiseNP[which(df$Outlier=="YES")]<-df$FreqMinor[which(df$Outlier=="YES")]
    
    # ML test extraction --
    ml.noise.mean<-mean(df$Noise, na.rm=TRUE)
    ml.noise.sd<-sd(df$Noise, na.rm=TRUE)
    ml.minor.mean<-mean(df$FreqMinor, na.rm=TRUE)
    ml.minor.sd<-sd(df$FreqMinor, na.rm=TRUE)
    ml.noise.q3<-as.numeric(quantile(df$FreqMinor, na.rm=TRUE)[4])
    ml.minor.q3<-as.numeric(quantile(df$Noise, na.rm=TRUE)[4])
    ml.count.missing<- length(which(df$Reads<20))/29903
    ml.minor.count<- length(which(df$Outlier=="YES"))/29903
    ml.co <- as.numeric(co)
    
    to.ml<-as.data.frame(t(c(ml.noise.mean, 
                             ml.noise.sd,
                             ml.minor.mean,
                             ml.minor.sd,
                             ml.noise.q3,
                             ml.minor.q3,
                             ml.count.missing,
                             ml.minor.count,
                             ml.co)))
    # ML extraction end --
    
    
    names<-gsub("\\.sorted.*","",gsub("_S[0-9].*","",gsub("R[0-9].*","",gsub(".*/","",results.files[i]))))
    to.ml$Sample<-names
    if(!exists("ml.out")){
      ml.out<-to.ml
    }else{
      ml.out<-rbind(ml.out, to.ml)
    }
    
   
    #names<-gsub("\\.primertrimmed.*","",gsub("_S[0-9].*","",gsub(".*/","",results.files[i])))
    out.plots[[i]]<-ggplot(df)+
      geom_line(aes(Base, NoiseNP))+
      geom_point(data=subset(df, Outlier=="YES"),aes(Base, NoiseNP),col="red", alpha=0.3)+
      geom_point(data=subset(df, Reads<20),aes(Base, 0),col="blue", alpha=0.1)+
      ylim(0,0.501)+
      theme_minimal()+
      ggtitle(names)+
      ylab("Minor Frequency")

  }
  

  
  
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

 pdf.list<-list.files(results, full.names = TRUE, pattern = ".*NoisExtractor.*\\.pdf")
 library("pdftools")
 if(length(pdf.list)>1){
 pdf_combine(pdf.list, output = gsub("_.\\.pdf","_Merged.pdf",pdf.list[1]))
 file.remove(pdf.list)}
}else{
  print("No bam files found in the /Noise/ folder")
}
