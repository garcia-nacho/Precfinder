library(ggpubr)
library("parallel")
library("foreach")
library(doSNOW)
library(progress)
library(ggplot2)
library(readxl)

df<-read.csv("/home/nacho/Documents/Corona/Recombinants/Norway_Omicron/df_11022022.csv")


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

cores<-8

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

pb<-txtProgressBar(min = 0, max = nrow(df), initial = 1)

df<-read.csv("/home/nacho/Documents/Corona/Recombinants/Oppsett248/Run721.csv",sep = ";")


try(rm(out.list))
try(rm(plot.list))



for (i in 1:nrow(df)) {

  mut.raw<-mutation.table$Mutation[which(mutation.table$Mutation %in%  unlist(strsplit(df$substitutions[i], ",") ))]
  if(length(mut.raw)>0){
  mut.total<-mutation.table$Lineage[which(mutation.table$Mutation %in%  unlist(strsplit(df$substitutions[i], ",") ))]
  probs<-unlist(mutation.table$Lineage.P[which(mutation.table$Mutation %in%  unlist(strsplit(df$substitutions[i], ",") ))])
  
  

    tab.dum<-table(mut.total)

    
    dum.list <- list(mut.raw, mut.total, probs, table(mut.total)) 
    names(dum.list)<-c("Mutation", "Lineage", "Probability", "Overview")
    
    
    to.plot<-as.data.frame(mut.raw)
    colnames(to.plot)<-"Mutation"
    if(length(probs)==0){
      probs<-1
      mut.total<-"Unknown"
    }
    to.plot$Probability<- probs
    to.plot$Lineage<- mut.total
    to.plot$Position<-as.numeric(gsub("^[A-Z]","", gsub("[A-Z]$","",to.plot$Mutation)))
    missing.df<-df$missing[i]
    if(missing.df!=""){
    missing.df<-unlist(strsplit(df$missing[i], ",") )
    
    to.clean<-missing.df[grep("-", missing.df)]
    if(length(to.clean)>0){
    final.clean<-vector()
    for (tc in 1:length(to.clean)) {
      temp.clean<-as.numeric(unlist(strsplit(to.clean[tc], "-") ))
      final.clean<-c(final.clean, c(temp.clean[1]:temp.clean[2]))
    }
    if(length(missing.df[-grep("-", missing.df)])>0) final.clean<-c(final.clean, as.numeric(missing.df[-grep("-", missing.df)]))
    missing.df<-as.data.frame(final.clean)
    }else{
      missing.df<-as.data.frame(missing.df)
    }
    
    colnames(missing.df)<-"Missing"
    
    }else{
    missing.df<-as.data.frame(29904)
    colnames(missing.df)<-"Missing"
    }  
    dum.list$Missing<-missing.df$Missing
    etp<-0
    try(etp<-round(entropy::entropy(aggregate(Probability~Lineage, to.plot[which(to.plot$Probability>0.8),], sum)$Probability),3))
    
    dum.plot<- ggplot(to.plot)+
            geom_segment(aes(y=0, yend= Probability, x=Position, xend=Position), stat = "identity",alpha=0.3)+
            geom_point(aes(Position, Probability, col=Lineage),alpha=0.5)+
            scale_color_manual(values =rainbow(length(unique(to.plot$Lineage))))+
            geom_point(data=missing.df, aes(as.numeric(Missing),0),col="blue")+
            ylim(0,1.001)+
      
            theme_minimal()+
            xlim(0,29903)+
      ylab("Prob.")+
      xlab("Position")+
      ggtitle(paste(gsub("/.*","",df$seqName[i]), " /E:", 
     etp,sep = ""
                    ))
      
    

    
    if(!exists("out.list")){
      out.list<-list(dum.list)
      plot.list<-list(dum.plot)

    }else{
      out.list<-c(out.list,list(dum.list))
      plot.list<-c(plot.list, list(dum.plot))

    }
    
    names(out.list)[length(out.list)]<-df$seqName[i]
    
  }
  
}

ggarrange(plotlist = plot.list[c(1:4)], ncol = 4, nrow = 1)



par.extractor<-function(input, co=0){
  try(rm(results.out))
for (i in 1:length(out.list)) {
  
  muts<-out.list[[i]]$Mutation
  muts<-as.numeric(gsub(".$","",gsub("^.","",muts)))
  probs<-out.list[[i]]$Probability[order(muts)]
  lin<-out.list[[i]]$Lineage[order(muts)]
  muts<-muts[order(muts)]
  
  if(co>0){    
    lin<-lin[which(probs>co)]
    muts<-muts[which(probs>co)]
    probs<-probs[which(probs>co)]
  }
  
  bits<-0
  bit.l<-vector()
  bit.n<-vector()
  
  etp70<-etp80<-etp90<-0
  missing<-out.list[[i]]$Missing
  if(length(missing[which(missing==29904)])==1){ 
    missing<-length(which(missing==29904))-1
  }else{
    missing<-length(missing)/29903
  }
  
  if(length(unique(lin))>1){
    start.bit<-1
    start.n<-muts[1]
    for (l in 2:length(lin)) {
     if(lin[l]!= lin[l-1]) {
       bits<-bits+1
       bit.l<-c(bit.l, l-start.bit)
       bit.n<-c(bit.n, muts[l]-start.n)
       start.bit<-l
       start.n<-muts[l]
     }
    } 
    bit.l<-c(bit.l, length(lin)-start.bit)/length(lin)
    bit.n<-c(bit.n, muts[length(muts)]-start.n)/29903
    
    if(bit.l[length(bit.l)]==0)bit.l<-bit.l[-length(bit.l)]
    if(bit.n[length(bit.n)]==0)bit.n<-bit.n[-length(bit.n)]
    
    to.plot<-as.data.frame(probs[which(probs>0.8)])
    colnames(to.plot)<-"Probability"
    to.plot$Lineage<-lin[which(probs>0.8)] 
    try(etp80<- entropy::entropy(aggregate(Probability~Lineage, to.plot, sum)$Probability))
    
    to.plot<-as.data.frame(probs[which(probs>0.7)])
    colnames(to.plot)<-"Probability"
    to.plot$Lineage<-lin[which(probs>0.7)] 
    try(etp70<- entropy::entropy(aggregate(Probability~Lineage, to.plot, sum)$Probability))
    
    to.plot<-as.data.frame(probs[which(probs>0.9)])
    colnames(to.plot)<-"Probability"
    to.plot$Lineage<-lin[which(probs>0.9)] 
    try(etp90<- entropy::entropy(aggregate(Probability~Lineage, to.plot, sum)$Probability))
    
    
  }else{
    bit.l<-1
    bit.n<-1
   
  }
 missig<-
 results<-as.data.frame(t(c( names(out.list)[i],
   bits, min(bit.l), max(bit.l), mean(bit.l), sd(bit.l), etp70, etp80, etp90,missing,length(lin)/29903
 )
   
 ))
  colnames(results)<-c("Sample","Bits","MinBitSize","MaxBitSize","MeanBitSize", "SdBitSize", "E70", "E80","E90","Missing","Muts")
  if(!exists("results.out")){
    results.out<-results
  }else{
    results.out<-rbind(results.out, results)
  }

  }
results.out$SdBitSize[which(is.na(results.out$SdBitSize))]<-0
return(results.out)
}


out0<-par.extractor(out.list)
colnames(out0)[c(2,3,4,5,6,11)]<-paste("N0.",colnames(out0)[c(2,3,4,5,6,11)],sep = "")

out8<-par.extractor(out.list, co=0.8)
colnames(out8)[c(2,3,4,5,6,11)]<-paste("N8.",colnames(out8)[c(2,3,4,5,6,11)],sep = "")
total<-merge(out0, out8[,c(1,2,3,4,5,6,11)], by="Sample")

out9<-par.extractor(out.list, co=0.9)
colnames(out9)[c(2,3,4,5,6,11)]<-paste("N9.",colnames(out9)[c(2,3,4,5,6,11)],sep = "")
total<-merge(total, out9[,c(1,2,3,4,5,6,11)], by="Sample")
