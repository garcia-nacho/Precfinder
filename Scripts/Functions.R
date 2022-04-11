library("ggpubr")
library("parallel")
library("foreach")
library("doSNOW")
library("progress")
library("ggplot2")
library("readxl")
library("Biostrings")
library("keras")
library("entropy")
library("abind")

# Nextclade & Pangolin running inside docker ---------------------------------
nc.pg.run <- function(collapse.lineage="AY", clean.bad=TRUE, freq.co=0.01, cov.co=0.95, mode="Training"){
  if(mode!="Training") mode<-"Inference"
  
  fa.train <- list.files(paste("/",mode,"/",sep = ""), full.names = TRUE, pattern = ".*\\.fa.*")
  if(length(fa.train)>0){
  fa.train<-readDNAStringSet(fa.train)
  names(fa.train)<-gsub(" ","_",names(fa.train))
  f.temp<-tempfile(tmpdir = paste("/",mode,sep = ""))
  #parsing large amount of sequences into chunks for nextstrain
  
  if(length(fa.train)<10000){
    writeXStringSet(fa.train, f.temp)
    system(paste("nextclade --input-fasta ", f.temp ," --output-csv ", paste(f.temp, "_nc.csv", sep = ""),sep = ""))
    system(paste("pangorunner.sh", f.temp))
    df.nc<-read.csv(paste(f.temp, "_nc.csv", sep = ""),sep = ";")
    df.pango<-read.csv(paste(f.temp, "_pango.csv", sep = ""))
    
    colnames(df.pango)[1]<-"seqName"
    df<-merge(df.nc, df.pango, by="seqName", all.x=TRUE)
    file.remove(paste(f.temp, "_pango.csv", sep = ""))
    file.remove(paste(f.temp, "_nc.csv", sep = ""))
    file.remove(f.temp)
  }else{
    continue<-TRUE

    start<-1
    end<-10000
    while (continue) {
      end<-min(end, length(fa.train))
      writeXStringSet(fa.train[c(start:end)], f.temp)
      system(paste("nextclade --input-fasta ", f.temp ," --output-csv ", paste(f.temp, "_nc.csv", sep = ""),sep = ""))
      system(paste("pangorunner.sh", f.temp))
      df.nc<-read.csv(paste(f.temp, "_nc.csv", sep = ""),sep = ";")
      df.pango<-read.csv(paste(f.temp, "_pango.csv", sep = ""))
      
      colnames(df.pango)[1]<-"seqName"
      df.temp<-merge(df.nc, df.pango, by="seqName", all.x=TRUE)
      file.remove(paste(f.temp, "_pango.csv", sep = ""))
      file.remove(paste(f.temp, "_nc.csv", sep = ""))
      file.remove(f.temp)
      
      if(start==1){
        df<-df.temp
      }else{
        df<-rbind(df,df.temp)
      }
      
      if(end==length(fa.train)) continue<-FALSE
      start <- end+1
      end <- end+10000
    }
  }
  #Remove none and NA from pangolin, Only for Training
  if(mode=="Training"){
  if(length(which(df$lineage=="None"))>0) df<-df[-which(df$lineage=="None"),]
  if(length(which(is.na(df$lineage)))>0) df<-df[-which(is.na(df$lineage)),]
  }
  
  #Collapse lineages
  if(collapse.lineage!=FALSE & mode=="Training"){
    collapse.lineage<-gsub(" ","",collapse.lineage)
    collapse.lineage<-unlist(base::strsplit(collapse.lineage,","))
    collapse.lineage<-paste(collapse.lineage, ".",sep = "")
    collapse.lineage<-gsub("\\.","\\\\.",collapse.lineage)
    for(i in 1:length(collapse.lineage)){
      if(length(grep(collapse.lineage[i], df$lineage))>0) df$lineage[grep(collapse.lineage[i], df$lineage)]<- paste(gsub("\\\\","",collapse.lineage[i]), "NN",sep = "")
    }
  }
  
  #Cut-off for lineage ratio < co (1%) only for training
  if(mode=="Training"){
  freqtable<-as.data.frame(table(df$lineage))
  freqtable$Freq<-freqtable$Freq/sum(freqtable$Freq)
  df<-df[which(df$lineage %in% freqtable$Var1[which(freqtable$Freq>=freq.co)]),]
  }
  
  cov.co <- 29903 - (29903*cov.co)
  if(length(which(df$totalMissing>cov.co))>0 & mode=="Training")df<-df[-which(df$totalMissing>cov.co),]
 
  #Recode insertions and deletions
  for (d in 1:nrow(df)) {
    #Deletions
    if(length(which(is.na(df$deletions)) )>0) df$deletions[which(is.na(df$deletions)) ]<-""
    if(df$deletions[d]!=""){
      dummy.del<- unlist(base::strsplit(df$deletions[d], ",") )
      dummyvec<-vector()
      for (i in 1:length(dummy.del)) {
        if(length(grep("-",dummy.del[i]))>0 ){
          start<-gsub("-.*","",dummy.del[i])
          end<-gsub(".*-","",dummy.del[i])
          size<-as.numeric(end)-as.numeric(start)
          recoded<-paste("D",start, toupper(letters)[ min(size,26)], sep = "")
          
          dummyvec<-c(dummyvec,recoded)
          
        }else{
          recoded<-paste("D",dummy.del[i],"A", sep = "")
          dummyvec<-c(dummyvec,recoded)
        }
      }
      df$substitutions[d]<-paste(df$substitutions[d], paste(dummyvec, collapse = ","), sep = ",")
    }
    
    #Insertions
   if(length(which(is.na(df$insertions)) )>0) df$insertions[which(is.na(df$insertions)) ]<-""

    if(df$insertions[d]!=""){
      dummy.ins<- unlist(base::strsplit(df$insertions[d], ",") )
      dummyvec<-vector()
      for (i in 1:length(dummy.ins)) {
        start<-gsub(":.*","",dummy.ins[i]) 
        size<-nchar(gsub(".*:","",dummy.ins))
        recoded<-paste("I",start, toupper(letters)[ min(size,26)], sep = "")
        dummyvec<-c(dummyvec,recoded)
      }
      df$substitutions[d]<-paste(df$substitutions[d], paste(dummyvec, collapse = ","), sep = ",")
    }  
  }
  
  write.csv(df, paste(paste("/",mode,"/",mode,"_dataset.csv", sep = "")),row.names = FALSE)
  return(NULL)
  }else{
    print("Error! No sequences found")
    return(NULL)
}
  }


# Probability Calculator of mutations from Training -----------------------
table.generator<-function(df="/Training/Training_dataset.csv", cores=10){
  df<-read.csv(df, stringsAsFactors = FALSE)
  
  df$WHO.lineage<-df$lineage
  
  mutation.table<-as.data.frame(table(unlist(strsplit(df$substitutions, ","))))
  mutation.table<-mutation.table[-which(mutation.table$Freq<3),]
  
  colnames(mutation.table)[1]<-"Mutation"
  
  add.lineages<-as.data.frame(matrix(data = NA, nrow = nrow(mutation.table),  ncol = (length(unique(df$WHO.lineage)))))
  colnames(add.lineages)<-unique(df$WHO.lineage)
  mutation.table<-cbind(mutation.table, add.lineages)
  mutation.table$Mutation<-as.character(mutation.table$Mutation)

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
    mutation.table$Lineage[i]<- colnames(add.lineages)[which(add.lineages[i,]==max(add.lineages[i,]))][1]
    if(mutation.table$Lineage[i]!="")  mutation.table$Lineage.P[i]<- add.lineages[i,which(add.lineages[i,]==max(add.lineages[i,]))][1]
    if(mutation.table$Lineage[i]=="")mutation.table$Lineage[i]<-NA
    
  }
  
  mutation.table$Mutation<-as.character(mutation.table$Mutation)
  mutation.table$Lineage.P<-unlist(mutation.table$Lineage.P)
  return(mutation.table)
}

# Probability calculator -------------------------------------------------
P.calculator <- function(input.data, mutation.table){

  df<-input.data
  colnames(df)[1]<-"seqName"
  if(length(which(colnames(df)=="missing"))==0) df$missing=""

  try(rm(out.list))
  try(rm(plot.list))
  pb<-txtProgressBar(min = 0, max = nrow(df), initial = 0)
  for (i in 1:nrow(df)) {
    setTxtProgressBar(pb,i)
    
    mut.raw<-mutation.table$Mutation[which(mutation.table$Mutation %in%  unlist(strsplit(df$substitutions[i], ",") ))]
    if(length(mut.raw)>0){
      mut.total<-mutation.table$Lineage[which(mutation.table$Mutation %in%  unlist(strsplit(df$substitutions[i], ",") ))]
      probs<-unlist(mutation.table$Lineage.P[which(mutation.table$Mutation %in%  unlist(strsplit(df$substitutions[i], ",") ))])
      
      if(length(mut.raw[-which(mut.raw %in% mutation.table$Mutation)])>0){
        missing<-mut.raw[-which(mut.raw %in% mutation.table$Mutation)]
        mut.raw<-mut.raw[which(mut.raw %in% mutation.table$Mutation)]
        mut.raw<-c(mut.raw, missing)
        mut.total<-c(mut.total, rep("Unknown",length(missing)))
        probs<-c(probs, rep(1,length(missing)))
      }
      
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
      try(etp<-round(entropy::entropy(aggregate(Probability~Lineage, to.plot[which(to.plot$Probability>0.8),], sum)$Probability),3), silent = TRUE)
      
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
      names(plot.list)[length(plot.list)]<-df$seqName[i]
      
    }
    
  }
  close(pb)

  out.to.export<-list(out.list, plot.list)
  names(out.to.export)<-c("Probs", "Plots")

  return(out.to.export)

  
}

# Dataset preparation -------------------------------------------------------------
trainset.prep <- function(data, sample.mult=5, max.number=200000){
  pb<-txtProgressBar(min=1, max = (nrow(data)*sample.mult*2), initial = 1)
continue<-TRUE
counter<-0
  while (continue){
  counter<-counter+1  
  setTxtProgressBar(pb,counter)
  dumm<-data[1,,drop=FALSE]
  dumm[1,]<-NA
  index.a<-sample(nrow(data),1)
  seq.a <- unlist(base::strsplit(data$substitutions[index.a], ","))
  lin.a <- data$lineage[index.a]
  index.b <-sample(which(data$lineage!=lin.a),1) 
  seq.b <- unlist(base::strsplit(data$substitutions[index.b], ","))
  n.cuts <- sample(c(rep(c(1),90),rep(c(2),6),rep(c(3),3)),1)
  
  cuts<- sample(c(600:29000), n.cuts)
  try(rm(seq.c))
  if(length(cuts)==1 & 
     length(which(as.numeric(gsub("^.","",gsub(".$","",seq.a))) <= cuts)) &
     length(which(as.numeric(gsub("^.","",gsub(".$","",seq.b))) > cuts))){   seq.c<-c(seq.a[which(as.numeric(gsub("^.","",gsub(".$","",seq.a))) <= cuts)],
                                 seq.b[which(as.numeric(gsub("^.","",gsub(".$","",seq.b))) > cuts)])
  }
  
  if(length(cuts)==2){
    
    if(length(which(as.numeric(gsub("^.","",gsub(".$","",seq.a))) <= cuts[1]))>0 &
       length(which(as.numeric(gsub("^.","",gsub(".$","",seq.b))) > cuts[1] &  as.numeric(gsub("^.","",gsub(".$","",seq.b))) <= cuts[2]))>0 &
       length(which(as.numeric(gsub("^.","",gsub(".$","",seq.a))) > cuts[2]))>0){
    
    seq.c<-c(seq.a[which(as.numeric(gsub("^.","",gsub(".$","",seq.a))) <= cuts[1])],
             seq.b[which(as.numeric(gsub("^.","",gsub(".$","",seq.b))) > cuts[1] &  as.numeric(gsub("^.","",gsub(".$","",seq.b))) <= cuts[2])],
             seq.a[which(as.numeric(gsub("^.","",gsub(".$","",seq.a))) > cuts[2])])
    }
  }
  
  if(length(cuts)==3){
    
    if(length(which(as.numeric(gsub("^.","",gsub(".$","",seq.a))) <= cuts[1]))>0 &
       length(which(as.numeric(gsub("^.","",gsub(".$","",seq.b))) > cuts[1] &  as.numeric(gsub("^.","",gsub(".$","",seq.b))) <= cuts[2]))>0 &
       length(which(as.numeric(gsub("^.","",gsub(".$","",seq.a))) > cuts[2]& as.numeric(gsub("^.","",gsub(".$","",seq.a))) <= cuts[3] ))>0 &
       length(which(as.numeric(gsub("^.","",gsub(".$","",seq.b))) > cuts[3]))>0){
      
      seq.c<-c(seq.a[which(as.numeric(gsub("^.","",gsub(".$","",seq.a))) <= cuts[1])],
               seq.b[which(as.numeric(gsub("^.","",gsub(".$","",seq.b))) > cuts[1] &  as.numeric(gsub("^.","",gsub(".$","",seq.b))) <= cuts[2])],
               seq.a[which(as.numeric(gsub("^.","",gsub(".$","",seq.a))) > cuts[2] & as.numeric(gsub("^.","",gsub(".$","",seq.a))) <= cuts[3] )],
               seq.b[which(as.numeric(gsub("^.","",gsub(".$","",seq.b))) > cuts[3])])
    }
  }
if(exists("seq.c")){
  seq.c<-paste(seq.c, collapse = ",") 
  dumm$substitutions<-seq.c
  dumm$breakpoints<-n.cuts
  dumm$breaksites<-paste(cuts, collapse = "/")
    
  if(!exists("out.ml")){
    out.ml<-dumm
  }else{
    out.ml<-rbind(out.ml, dumm)
  }
  if(nrow(out.ml)==min((nrow(data)*sample.mult),max.number)) continue<-FALSE
}
  
  }
close(pb)
  out.ml$Class<-"Recombinant"
  data$Class<-"Non-Recombinant"
  data$breakpoints<-0
  data$breaksites<-0
  out.ml$seqName<-paste("Rec",rownames(out.ml),sep="_")
  data<-rbind(data, out.ml)
  
  data<-data[sample(1:nrow(data),nrow(data)),c("seqName","substitutions","Class","breakpoints","breaksites")]
  rownames(data)<-c(1:nrow(data))
  colnames(data)[1]<-"Sample"
  return(data)
  
} 

# MLCNN -------------------------------------------------------------------
cnn.training<- function(training.set, modelid){
  testCNN<-training.set$Probs
  mut.n<-unlist(lapply(testCNN, function(x) length(x$Mutation)))
  lin<-unique(mutation.table$Lineage)
  
  train.array<-array(0, dim = c(length(testCNN), max(mut.n), length(lin)))
  labelmatrix<-matrix(0, nrow = length(testCNN), ncol = 2)
  pb<-txtProgressBar(min = 1, max = length(testCNN), initial = 1)
  
  for (i in 1:length(testCNN)) {
    setTxtProgressBar(pb,i)
    dum.mut<-testCNN[[i]]$Mutation
    dum.lin<-testCNN[[i]]$Lineage
    dum.prb<-testCNN[[i]]$Probability
    if(length(grep("^Rec_",names(testCNN)[i]))>0){
      labelmatrix[i,1]<-1
    }else{
      labelmatrix[i,2]<-1
    }
    
    #Order
    dum.lin<-dum.lin[order(as.numeric(gsub("^.","",gsub(".$","",dum.mut))))]
    dum.prb<-dum.prb[order(as.numeric(gsub("^.","",gsub(".$","",dum.mut))))]
    dum.mut<-dum.mut[order(as.numeric(gsub("^.","",gsub(".$","",dum.mut))))]
    
    for (j in 1:length(lin)) {
      if(length(which(dum.lin==lin[j]))>0)
        train.array[i,which(dum.lin==lin[j]),j]<-dum.prb[which(dum.lin==lin[j])]
    }
  }
  
  #Augmentation 4X for lineages number > 3
  
  if(length(lin)>2){
    
    temp.array<-train.array[,,1]
    train.array.agg<-train.array 
    train.array.agg[,,1]<-train.array.agg[,,2]
    train.array.agg[,,2]<-temp.array
    
    train.array<- abind(train.array, train.array.agg, along = 1 )
    labelmatrix<-rbind(labelmatrix, labelmatrix)
    
    temp.array<-train.array[,,3]
    train.array.agg<-train.array 
    train.array.agg[,,3]<-train.array.agg[,,2]
    train.array.agg[,,2]<-temp.array
    
    train.array<- abind(train.array, train.array.agg, along = 1 )
    labelmatrix<-rbind(labelmatrix, labelmatrix)
    
  }
  
  
  
  size.dataset<-as.data.frame(c(dim(train.array)[2], lin))
  colnames(size.dataset)<-"Parameters"
  write.csv(size.dataset, paste("/Models/",modelid,"/",gsub("-","",Sys.Date()),"_InputSize.csv",sep = ""))
  #Split validation
  val.index<-sample(1:dim(train.array)[1], round(dim(train.array)[1]/10))
  val.array<-train.array[val.index,,]
  val.labs<-labelmatrix[val.index,]
  train.array<-train.array[-val.index,,]
  labelmatrix<-labelmatrix[-val.index,]
  
  model <- keras_model_sequential() %>% 
    layer_conv_1d(filters = 64, kernel_size = 5, activation = "relu", input_shape = c(dim(train.array)[-1])) %>% 
    layer_max_pooling_1d(pool_size =2) %>% 
    layer_conv_1d(filters = 32, kernel_size = 3, activation = "relu") %>% 
    layer_max_pooling_1d(pool_size =2) %>% 
    layer_conv_1d(filters = 12, kernel_size = 3, activation = "relu") %>% 
    layer_max_pooling_1d(pool_size =2) %>%
    layer_flatten() %>% 
    layer_dense(units = 24) %>% 
    layer_dense(units = 12) %>% 
    layer_dense(units = 2, activation = "softmax") 
  
  summary(model)
  
  model %>% compile(
    loss = 'binary_crossentropy',
    optimizer =optimizer_adagrad(),
    metrics = c('accuracy')
  )
  
  ES<-list(callback_early_stopping(
    monitor = "val_accuracy",
    min_delta = 0,
    patience = 8,
    verbose = 0,
    mode =  "max",
    restore_best_weights = TRUE
  ))
  
  
  history <- model %>% fit(
    train.array,labelmatrix , 
    epochs = 60, batch_size = 128, 
    #class_weight=list("0"=1, "1"=sum(labelmatrix[,2])/sum(labelmatrix[,1])),
    class_weight=list("0"=1, "1"=5),
    callbacks=ES,
    shuffle=TRUE,
    validation_split = 0.1
  )
  
  predictions<-model %>% predict(val.array)
  predictions<-as.data.frame(predictions)
  
  val.labs<-as.data.frame(val.labs)
  val.labs$Class<-"Non-Recombinant"
  val.labs$Class[which(val.labs$V1==1)]<-"Recombinant"
  
  predictions$Observed<-val.labs$Class
  FP <- length(which(predictions$V1>0.5 &  predictions$Observed=="Non-Recombinant"))
  FN <- length(which(predictions$V2>0.5 &  predictions$Observed=="Recombinant"))
  TP <- length(which(predictions$V1>0.5 &  predictions$Observed=="Recombinant"))
  TN <- length(which(predictions$V2>0.5 &  predictions$Observed=="Non-Recombinant"))
  TNR<-TN/(TN+FP)
  Recall<- TP/(TP+FN)
  Precision<- TP/(TP+FP)
  F1=2*((Precision*Recall)/(Precision+Recall))
  BA<-(TNR+Recall)/2
  
  cm<-as.data.frame(c(TP,FP,FN,TN))
  colnames(cm)<-"Value"
  cm$Observed<-c("Rec","NonRec","Rec","NonRec")
  cm$Predicted<-c("Rec","Rec","NonRec","NonRec")
  
  ggplot(cm)+
    geom_tile(aes(Observed, Predicted, fill=Value), show.legend=FALSE )+
    scale_fill_gradient(low="blue", high = "red")+
    geom_text(aes(Observed, Predicted, label=Value), col="white")+
    theme_minimal()+
    ggtitle(paste("Balanced Accuracy:", round(BA,3)))
  ggsave(paste("/Models/",modelid,"/",gsub("-","",Sys.Date()),"_Model_CM.pdf",sep = ""), width = 4, height = 4)
  
  plot(history)+
    theme_minimal()+
    ggtitle(paste("Training process Model",gsub("-","",Sys.Date()) ,sep=""))
  ggsave(paste("/Models/",modelid,"/",gsub("-","",Sys.Date()),"_Model_Training.pdf",sep = ""), width = 4, height = 6)
  
  
  #Save model
  
    path<-paste("/Models/",modelid,"/",gsub("-","",Sys.Date()),"_Recombinant_Model.hdf5",sep = "")
    save_model_hdf5(model,path)

  return(model)
  
}

# ML Inference ------------------------------------------------------------
ml.inference<-function(inference.set, model, model.id){
  
  size.file<-list.files(paste("/Models/", model.id,sep = ""), pattern = "_InputSize.csv", full.names = TRUE)
  size.dataset<-read.csv(size.file)
  
  labels.to.test<- names(inference.set$Probs)
  testCNN<-inference.set$Probs
  
  mut.n<-unlist(lapply(testCNN, function(x) length(x$Mutation)))
  lin<-size.dataset$Parameters[-1]
  
  train.array<-array(0, dim = c(length(testCNN), as.numeric(size.dataset$Parameters[1]), length(lin)))
  
  pb<-txtProgressBar(min = 1, max = length(testCNN), initial = 1)
  warn.list.size<-vector()
  warn.list.lin<-vector()
  for (i in 1:length(testCNN)) {
    setTxtProgressBar(pb,i)
    dum.mut<-testCNN[[i]]$Mutation
    dum.lin<-testCNN[[i]]$Lineage
    dum.prb<-testCNN[[i]]$Probability
    
    #Order
    dum.lin<-dum.lin[order(as.numeric(gsub("^.","",gsub(".$","",dum.mut))))]
    dum.prb<-dum.prb[order(as.numeric(gsub("^.","",gsub(".$","",dum.mut))))]
    dum.mut<-dum.mut[order(as.numeric(gsub("^.","",gsub(".$","",dum.mut))))]
    
    if(length(dum.lin)>as.numeric(size.dataset$Parameters[1])){
      dum.lin<-dum.lin[1, dim(train.array)[2]]
      dum.prb<-dum.prb[1, dim(train.array)[2]]
      dum.mut<-dum.mut[1, dim(train.array)[2]]
      warn.list.size<-c(warn.list.size, i)
    }
    if(length(unique(dum.lin)[-which(unique(dum.lin) %in% lin)])>0) warn.list.lin<-c(warn.list.lin, i)
    for (j in 1:length(lin)) {
      if(length(which(dum.lin==lin[j]))>0)
        train.array[i,which(dum.lin==lin[j]),j]<-dum.prb[which(dum.lin==lin[j])]
    }
  }
  
  
  inf.array<-train.array
  
  predictions<-model %>% predict(inf.array)
  predictions<-as.data.frame(predictions)
  colnames(predictions)<-c("Score.Recombinant", "Score Non-Recombinant")
  predictions$Sample<-labels.to.test
  predictions<-predictions[,c(3,1)]
  predictions$Class<-"Non Recombinant"
  predictions$Class[which(predictions$Score.Recombinant>0.5)]<-"Recombinant"
  predictions$Comments<-NA
  
  if(length(warn.list.lin)!=0 & length(warn.list.size)!=0){
  warn.full <- base::intersect(warn.list.size, warn.list.lin)
  if(length(warn.full)>0){
    warn.list.size<-warn.list.size[-warn.full]
    warn.list.lin<-warn.list.lin[-warn.full]
  }
  }
  if(length(warn.list.size)>0) predictions$Comments[warn.list.size]<-"Warning: Last mutations were non included"
  if(length(warn.list.lin)>0) predictions$Comments[warn.list.lin]<-"Warning: Lineages not seen on training"
  if(exists("warn.full")) predictions$Comments[warn.full]<-"Warning: Lineages not seen on training & last mutations were not included"
  
  return(predictions)
  
}




