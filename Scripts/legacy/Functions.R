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

# Parameter extractor -----------------------------------------------------
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
      try(etp80<- entropy::entropy(aggregate(Probability~Lineage, to.plot, sum)$Probability),silent = TRUE)
      
      to.plot<-as.data.frame(probs[which(probs>0.7)])
      colnames(to.plot)<-"Probability"
      to.plot$Lineage<-lin[which(probs>0.7)] 
      try(etp70<- entropy::entropy(aggregate(Probability~Lineage, to.plot, sum)$Probability),silent = TRUE)
      
      to.plot<-as.data.frame(probs[which(probs>0.9)])
      colnames(to.plot)<-"Probability"
      to.plot$Lineage<-lin[which(probs>0.9)] 
      try(etp90<- entropy::entropy(aggregate(Probability~Lineage, to.plot, sum)$Probability),silent = TRUE)
      
      
    }else{
      bit.l<-1
      bit.n<-1
      
    }
    
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

# Nextclade & Pangolin running inside docker ---------------------------------
nc.pg.run <- function(collapse.lineage=FALSE, clean.bad=TRUE, freq.co=0.01, cov.co=0.95, mode="Training"){
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
    df<-merge(df.nc, df.pango, by="seqName")
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
      df<-merge(df.nc, df.pango, by="seqName",all.x=TRUE)
      file.remove(paste(f.temp, "_pango.csv", sep = ""))
      file.remove(paste(f.temp, "_nc.csv", sep = ""))
      file.remove(f.temp)
      
      if(!exists("df")){
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
  if(clean.bad & mode=="Training"){
  if(length(which(df$lineage=="None"))>0) df<-df[-which(df$lineage=="None"),]
  if(length(which(is.na(df$lineage)))>0) df<-df[-which(is.na(df$lineage)),]
  }
  
  #Collapse lineages
  if(collapse.lineage!=FALSE){
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
  if(length(which(df$totalMissing>cov.co))>0)df<-df[-which(df$totalMissing>cov.co),]
  
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
P.calculator <- function(input.data, mutation.table, ML=TRUE){

  df<-input.data
  colnames(df)[1]<-"seqName"
  if(ML) df$missing=""
  try(rm(out.list))
  try(rm(plot.list))
  pb<-txtProgressBar(min = 1, max = nrow(df), initial = 1)
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
      
      
    }
    
  }
  close(pb)
  if(ML){
  environment(par.extractor) <- environment()
  out0<-par.extractor(out.list)
  colnames(out0)[c(2,3,4,5,6,11)]<-paste("N0.",colnames(out0)[c(2,3,4,5,6,11)],sep = "")
  
  out8<-par.extractor(out.list, co=0.8)
  colnames(out8)[c(2,3,4,5,6,11)]<-paste("N8.",colnames(out8)[c(2,3,4,5,6,11)],sep = "")
  total.ml<-merge(out0, out8[,c(1,2,3,4,5,6,11)], by="Sample")
  
  out9<-par.extractor(out.list, co=0.9)
  colnames(out9)[c(2,3,4,5,6,11)]<-paste("N9.",colnames(out9)[c(2,3,4,5,6,11)],sep = "")
  total.ml<-merge(total.ml, out9[,c(1,2,3,4,5,6,11)], by="Sample")
  
  total.ml<-as.matrix(total.ml)
  total.ml[which(is.na(total.ml))]<-0
  total.ml<-as.data.frame(total.ml)
  
  total.ml[,-1]<-apply(total.ml[,-1], 2, function(x) as.numeric(as.character(x)))
  
  total.ml[,grep("Bits",colnames(total.ml))]<-total.ml[,grep("Bits",colnames(total.ml))]/29903
  }
  out.to.export<-list(out.list, plot.list)
  names(out.to.export)<-c("Probs", "Plots")
  if(ML){
    return(total.ml)
  }else{
    return(out.to.export)
  }
  
}

# Dataset preparation -------------------------------------------------------------
trainset.prep <- function(data, sample.mult=5, max.number=80000){
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

# ML training -------------------------------------------------------------
ml.training<- function(ml.set, training.set){

  ml.set<-merge(ml.set,training.set,by="Sample")
  if(length(which(ml.set$N0.Bits==0 & ml.set$N8.Bits==0 & ml.set$N9.Bits==0))>0){
  ml.set<-ml.set[-which(ml.set$N0.Bits==0 & ml.set$N8.Bits==0 & ml.set$N9.Bits==0),]
  }
  
  training.index<-sample(1:nrow(ml.set), round(nrow(ml.set)*0.9))
  training <- as.matrix(ml.set[training.index,-c(1,10,24,25,26,27)])
  validation <- as.matrix(ml.set[-training.index,-c(1,10,24,25,26,27)])  
  
  ml.label<-as.data.frame(ml.set$Class)
  colnames(ml.label)<-"Class"
  ml.label$Recombinant<-0
  ml.label$Non.Recombinant<-0
  ml.label$Recombinant[which(ml.label$Class=="Recombinant")]<-1
  ml.label$Non.Recombinant[which(ml.label$Class!="Recombinant")]<-1
  ml.label$Class<-NULL
  
  label.training<-ml.label[training.index,]
  
  label.validation<-ml.label[-training.index,]
  

   #Model creation
  model <- keras_model_sequential() 
  model %>% 
    layer_dense(units = 20, activation = 'relu', input_shape = c(21)) %>% 
    layer_dropout(rate = 0.4) %>% 
    layer_dense(units = 20, activation = 'relu') %>%
    layer_dropout(rate = 0.3) %>%
    layer_dense(units = 10, activation = 'relu') %>%
    layer_dense(units = 2, activation = 'softmax')
  
  model %>% compile(
    loss = 'binary_crossentropy',
    optimizer =optimizer_adagrad(learning_rate = 0.001),
    
    metrics = c('accuracy')
  )
  
  
  ES<-list(callback_early_stopping(
    monitor = "val_accuracy",
    min_delta = 0,
    patience = 10,
    verbose = 0,
    mode =  "max",
    restore_best_weights = TRUE
  ))
  
  
  #Early stop
  history <- model %>% fit(
    as.matrix(training), as.matrix(label.training), 
    epochs = 30, batch_size = 128, 
    class_weight=list("0"=1, "1"=as.numeric(table(label.training$Recombinant))[2]/as.numeric(table(label.training$Recombinant))[1]),
    #callbacks=ES,
    validation_split = 0.1
    )

  predictions<-model %>% predict(validation)
  predictions<-as.data.frame(predictions)
  predictions$Observed<-ml.set$Class[-training.index]
  predictions$Mutations<-ml.set$substitutions[-training.index]
  
  #Confusion matrix
  #Positive -> Recombinant
  FP <- length(which(predictions$V1>0.5 &  predictions$Observed=="Non-Recombinant"))
  FN <- length(which(predictions$V2>0.5 &  predictions$Observed=="Recombinant"))
  TP <- length(which(predictions$V1>0.5 &  predictions$Observed=="Recombinant"))
  TN <- length(which(predictions$V2>0.5 &  predictions$Observed=="Non-Recombinant"))
  TNR<-TN/(TN+FP)
  Recall<- TP/(TP+FN)
  Precision<- TP/(TP+FP)
  F1=2*((Precision*Recall)/(Precision+Recall))
  BA<-(TNR+Recall)/2
  #print(paste("Balanced Accuracy:",BA))
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
  ggsave(paste("/Models/",gsub("-","",Sys.Date()),"_Model_CM.pdf",sep = ""), width = 4, height = 4)
  
  plot(history)+
    theme_minimal()+
    ggtitle(paste("Training process Model",gsub("-","",Sys.Date()) ,sep=""))
  ggsave(paste("/Models/",gsub("-","",Sys.Date()),"_Model_Training.pdf",sep = ""), width = 4, height = 6)
  
  
  if(length(c(which(predictions$V2>0.5 &  predictions$Observed=="Recombinant"),which(predictions$V1>0.5 &  predictions$Observed=="Non-Recombinant")))>0){
  
  errors<-predictions[c(which(predictions$V2>0.5 &  predictions$Observed=="Recombinant"),
                        which(predictions$V1>0.5 &  predictions$Observed=="Non-Recombinant")),]
  colnames(errors)[4]<-"substitutions"
  errors$missing<-""
  results.errors<-P.calculator(input.data = errors, mutation.table = mutation.table, ML=FALSE)
  ggarrange(plotlist = results.errors$Plots[c(1:min(40, length(results.errors$Plots)))], ncol = 4, nrow = 10)
  ggsave(paste("/Models/",gsub("-","",Sys.Date()),"_Model_ErrorExamples.pdf",sep = ""), width = 15, height = 30)
  }
  if(length(which(predictions$V1>0.5 &  predictions$Observed=="Recombinant"))>0){
  corrects<-predictions[which(predictions$V1>0.5 &  predictions$Observed=="Recombinant"),]
  colnames(corrects)[4]<-"substitutions"
  corrects$missing<-""
  results.corrects<-P.calculator(input.data = corrects, mutation.table = mutation.table, ML=FALSE)
  ggarrange(plotlist = results.corrects$Plots[c(1:min(40, length(results.corrects$Plots)))], ncol = 4, nrow = 10)
  ggsave(paste("/Models/",gsub("-","",Sys.Date()),"_Model_CorrectExamples.pdf",sep = ""), width = 15, height = 30)
  }
  
  #Save model
  if(BA > 0.6){
  path<-paste("/Models/",gsub("-","",Sys.Date()),"_Recombinant_ModelW.hdf5",sep = "")
  save_model_hdf5(model,path)
  }else{
    print("Warning: No meaningful model generated. Nothing was saved")
  }
  return(model)
}

# ML Inference ------------------------------------------------------------
ml.inference<-function(ml.set, model){

  labels<- ml.set$Sample
  
  if(length(which(ml.set$N0.Bits==0 & ml.set$N8.Bits==0 & ml.set$N9.Bits==0))>0){
  inference.to.test<-as.matrix(ml.set[-which(ml.set$N0.Bits==0 & ml.set$N8.Bits==0 & ml.set$N9.Bits==0),-c(1,10)])
  labels.clean<-labels[which(ml.set$N0.Bits==0 & ml.set$N8.Bits==0 & ml.set$N9.Bits==0)]
  labels.to.test<-labels[-which(ml.set$N0.Bits==0 & ml.set$N8.Bits==0 & ml.set$N9.Bits==0)]
  }else{
    inference.to.test<-as.matrix(ml.set)
    labels.to.test<-labels
  }
  
  predictions<-model %>% predict(inference.to.test)
  predictions<-as.data.frame(predictions)
  colnames(predictions)<-c("Score.Recombinant", "Score Non-Recombinant")
  predictions$Sample<-labels.to.test
  predictions<-predictions[,c(3,1)]
  predictions$Class<-"Non Recombinant"
  predictions$Class[which(predictions$Score.Recombinant>0.5)]<-"Recombinant"
  if(length(labels.clean)>0){
    labels.clean<-as.data.frame(labels.clean)
    colnames(labels.clean)<-"Sample"
    labels.clean$Score.Recombinant<-0
    labels.clean$Class<-"Non-Recombinant"
    predictions<-rbind(predictions, labels.clean)
  }
  
  
  return(predictions)
  
}


