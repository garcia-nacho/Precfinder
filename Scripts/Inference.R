source("/home/docker/Scripts/Functions.R")

args=commandArgs(TRUE)

if(dir.exists(paste("/Models/", args[1],sep = ""))){

nc.pg.run(mode="Inference")
mutation.table.file<-list.files(paste("/Models/", args[1],sep = ""), pattern = ".*_MutationTable.csv", full.names = TRUE)
mutation.table<- read.csv(mutation.table.file, stringsAsFactors = FALSE)

inference.raw<-read.csv("/Inference/Inference_dataset.csv")
resultsInference<-P.calculator(input.data = inference.raw, mutation.table = mutation.table)

model.file<-mutation.table.file<-list.files(paste("/Models/", args[1],sep = ""), pattern = ".*.hdf5", full.names = TRUE)
model<-load_model_hdf5(model.file)
results.out<-ml.inference(inference.set = resultsInference, model = model ,model.id=args[1])

library(writexl)
write_xlsx(results.out,"/Inference/Recombinant_Prediction.xlsx")
#Part for getting pdfs from inference
out.plots<-resultsInference$Plots
if(length(out.plots)<=40 ){
  if(length(out.plots)>0){
    ggarrange(plotlist =  out.plots[1:length(out.plots)], ncol = 4, nrow = 10)
    ggsave(paste("/Inference/",gsub("-","",Sys.Date()),"_RecombinantPlots.pdf",sep = ""), width = 420, height = 600, units = "mm") #A4
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
    ggsave(paste("/Inference/",gsub("-","",Sys.Date()),"_RecombinantPlots_",counter,".pdf",sep = ""), width = 420, height = 600, units = "mm") #A4
  }
}

library("pdftools")

pdf.list<-list.files("/Inference/", full.names = TRUE, pattern = ".*RecombinantPlots.*\\.pdf")
if(length(pdf.list)>1){
  pdf_combine(pdf.list, output = paste("/Inference/",gsub("-","",Sys.Date()),"_RecombinantPlots_MultiPage.pdf",sep = ""))
  file.remove(pdf.list)
}



# PlotOnly Recombinants -------------------------------------------------------------
pangonc<-read.csv("/Inference/Inference_dataset.csv")
colnames(pangonc)[1]<-"Sample"
results.out2<-merge(results.out, pangonc[,c("Sample","lineage")], by="Sample", all.x = TRUE)
write_xlsx(results.out2,"/Inference/Recombinant_Prediction_wPango.xlsx")

if(length(which(results.out$Class=="Recombinant"))>0){

out.plots<-out.plots[which(names(out.plots) %in% results.out$Sample[which(results.out$Class=="Recombinant")]) ]
if(length(out.plots)<=40 ){
  if(length(out.plots)>0){
    ggarrange(plotlist =  out.plots[1:length(out.plots)], ncol = 4, nrow = 10)
    ggsave(paste("/Inference/",gsub("-","",Sys.Date()),"_RecombinantSelected.pdf",sep = ""), width = 420, height = 600, units = "mm") #A4
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
    ggsave(paste("/Inference/",gsub("-","",Sys.Date()),"_RecombinantSelected_",counter,".pdf",sep = ""), width = 420, height = 600, units = "mm") #A4
  }
}


pdf.list<-list.files("/Inference/", full.names = TRUE, pattern = ".*RecombinantSelected.*\\.pdf")
if(length(pdf.list)>1){
  pdf_combine(pdf.list, output = paste("/Inference/",gsub("-","",Sys.Date()),"_RecombinantSelected_MultiPage.pdf",sep = ""))
  file.remove(pdf.list)
}

}
}else{
  print(paste("No model found. Available models: ", paste(gsub(".*/","",list.dirs("/Models/",recursive = FALSE)), collapse = " / "), sep = ""))
}