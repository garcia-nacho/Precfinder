library("ggpubr")
library("parallel")
library("foreach")
library("doSNOW")
library("progress")
library("ggplot2")
library("readxl")
library("Biostrings")

# Nextclade & Pangolin running inside docker ---------------------------------

nc.pg.run <- function(collapse.lineage=FALSE, clean.bad=TRUE, freq.co=0.01, cov.co=0.95){
  fa.train <- list.files("/Training/", full.names = TRUE, pattern = ".*\\.fa.*")
  if(length(fa.train)>0){
  fa.train<-readDNAStringSet(fa.train)
  names(fa.train)<-gsub(" ","_",names(fa.train))
  
  #parsing large amount of sequences into chunks for nextstrain
  if(length(fa.train)<10000){
    writeXStringSet(fa.train, "/Training/Training.fasta")
    system("nextclade --input-fasta /Training/Training.fasta --output-csv /Training/Training_nc.csv")
    system(paste("pangorunner.sh /Training/Training.fasta"))
    df.nc<-read.csv("/Training/Training_nc.csv",sep = ";")
    df.pango<-read.csv("/Training/Training_pango.csv")
    colnames(df.pango)[1]<-"seqName"
    df<-merge(df.nc, df.pango, by="seqName")
    file.remove("/Training/Training_pango.csv")
    file.remove("/Training/Training_nc.csv")
    file.remove("/Training/Training.fasta")
  }else{
    continue<-TRUE

    start<-1
    end<-10000
    while (continue) {
      end<-min(end, length(fa.train))
      writeXStringSet(fa.train[c(start:end)], "/Training/Training.fasta")
      system("nextclade --input-fasta /Training/Training.fasta --output-csv /Training/Training_nc.csv")
      system(paste("pangorunner.sh /Training/Training.fasta"))
      df.nc<-read.csv("/Training/Training_nc.csv",sep = ";",stringsAsFactors = FALSE)
      df.pango<-read.csv("/Training/Training_pango.csv",stringsAsFactors = FALSE)
      colnames(df.pango)[1]<-"seqName"
      df.temp<-merge(df.nc, df.pango, by="seqName")
      file.remove("/Training/Training_pango.csv")
      file.remove("/Training/Training_nc.csv")
      file.remove("/Training/Training.fasta")
      
      if(!exists("df")){
        df<-df.temp
      }else{
        df<-rbind(df,df.temp)
      }
      
      if(end==length(fa.train)) continue<-FALSE
      start <- end+1
      end <- end+10000
    }
    if(length(which(df$lineage=="None"))>0) df<-df[-which(df$lineage=="None"),]
    if(length(which(is.na(df$lineage)))>0) df<-df[-which(is.na(df$lineage)),]
    
   #Collapse.
    if(collapse.lineage!=FALSE){
      collapse.lineage<-gsub(" ","",collapse.lineage)
      collapse.lineage<-unlist(base::strsplit(collapse.lineage,","))
      collapse.lineage<-paste(collapse.lineage, ".",sep = "")
      collapse.lineage<-gsub("\\.","\\\\.",collapse.lineage)
      for(i in 1:length(collapse.lineage)){
        if(length(grep(collapse.lineage[i], df$lineage))>0) df$lineage[grep(collapse.lineage[i], df$lineage)]<- paste(gsub("\\\\","",collapse.lineage[i]), "NN",sep = "")
      }
    }
    freqtable<-as.data.frame(table(df$lineage))
    freqtable$Freq<-freqtable$Freq/sum(freqtable$Freq)
    
    df<-df[which(df$lineage %in% freqtable$Var1[which(freqtable$Freq>=freq.co)]),]
    cov.co <- 29903 - (29903*cov.co)
    if(length(which(df$totalMissing>cov.co))>0)df<-df[-which(df$totalMissing>cov.co),]
    
  }
  write.csv(df, "/Training/Training_dataset.csv",row.names = FALSE)
  return(NULL)
  }else{
    print("Error! No sequences found")
    return(NULL)
}
  }



# Probability Calculator of mutations from Training -----------------------

P.calc<-function(df="/Training/Training_dataset.csv", cores=10){
  df<-read.csv(df)
  
  df$WHO.lineage<-df$lineage
  
  mutation.table<-as.data.frame(table(unlist(strsplit(df$substitutions, ","))))
  mutation.table<-mutation.table[-which(mutation.table$Freq<5),]
  
  unique(df$WHO.lineage)
  colnames(mutation.table)[1]<-"Mutation"
  
  add.lineages<-as.data.frame(matrix(data = NA, nrow = nrow(mutation.table),  ncol = (length(unique(df$WHO.lineage)))))
  colnames(add.lineages)<-unique(df$WHO.lineage)
  mutation.table<-cbind(mutation.table, add.lineages)
  mutation.table$Mutation<-as.character(mutation.table$Mutation)
  pb<-txtProgressBar(min = 1, max = nrow(mutation.table)*(ncol(mutation.table)-3), initial = 1)

  samp<-mutation.table$Mutation
  
  pb <- progress_bar$new(
    format = "Mutation: :samp.pb [:bar] :elapsed | eta: :eta",
    total = nrow(mutation.table),    # 100 
    width = 60)
  
  progress <- function(n){
    pb$tick(tokens = list(samp.pb = samp[n]))
  } 
  
  opts <- list(progress = progress)
  
  
  cluster.cores<-makeCluster(cores)
  registerDoSNOW(cluster.cores)
  
  
  out.par<-foreach(i=1:nrow(mutation.table), .verbose=FALSE, .options.snow = opts) %dopar%{
    #dummy<-mutation.table[i,j, drop=FALSE]
    for (j in 3:ncol(mutation.table)) {
      mutation.table[i,j]<-length(intersect(grep(mutation.table$Mutation[i], df$substitutions), which(df$WHO.lineage==colnames(mutation.table)[j]))  )  
      
    }
    return(mutation.table[i,])
  }
  stopCluster(cluster.cores)
  
  mutation.table<-do.call("rbind", out.par)
  rm(out.par)
  
  
  for (i in 1:ncol(add.lineages)) {
    add.lineages[,i]<-   mutation.table[,which(colnames(mutation.table)==colnames(add.lineages)[i])]/length(which(df$WHO.lineage==colnames(add.lineages)[i])) 
    
  }
  
  colnames(add.lineages)<-paste("P.",unique(df$WHO.lineage),sep = "")
  mutation.table<-cbind(mutation.table, add.lineages)
  
  add.lineages<-as.data.frame(matrix(data = NA, nrow = nrow(mutation.table),  ncol = (length(unique(df$WHO.lineage)))))
  colnames(add.lineages)<-unique(df$WHO.lineage)
  
  for (i in 1:nrow(mutation.table)) {
    for (j in 1:ncol(add.lineages)) {
      
      p.mut.lin <- mutation.table[i,which(colnames(mutation.table)==colnames(add.lineages)[j])]/
        length(which(df$WHO.lineage==colnames(add.lineages)[j]))
      p.mut <- mutation.table$Freq[i] / nrow(df)
      p.lin <- length(which(df$WHO.lineage==colnames(add.lineages)[j]))/nrow(df)
      #Bayes rule
      add.lineages[i,j]<- (p.mut.lin*p.lin)/p.mut 
    }
  }
  
  mutation.table$Lineage<-NA
  mutation.table$Lineage.P<-NA
  
  for (i in 1:nrow(mutation.table)) {
    mutation.table$Lineage[i]<- paste(colnames(add.lineages)[which(add.lineages[i,]==max(add.lineages[i,]))], collapse = "/")
    if(mutation.table$Lineage[i]!="")  mutation.table$Lineage.P[i]<- add.lineages[i,which(add.lineages[i,]==max(add.lineages[i,]))]
    if(mutation.table$Lineage[i]=="")mutation.table$Lineage[i]<-NA
    
  }
  
  mutation.table$Mutation<-as.character(mutation.table$Mutation)
  
  mutation.table$Lineage<-NA
  mutation.table$Lineage.P<-NA
  
  for (i in 1:nrow(mutation.table)) {
    mutation.table$Lineage[i]<- paste(colnames(add.lineages)[which(add.lineages[i,]==max(add.lineages[i,]))], collapse = "/")
    if(mutation.table$Lineage[i]!="")  mutation.table$Lineage.P[i]<- add.lineages[i,which(add.lineages[i,]==max(add.lineages[i,]))]
    if(mutation.table$Lineage[i]=="")mutation.table$Lineage[i]<-NA
    
  }
  
  mutation.table$Mutation<-as.character(mutation.table$Mutation)
  
  return(mutation.table)
}




# ML training -------------------------------------------------------------

ml.training <- function(){
  
} 
