source("/home/docker/Scripts/Functions.R")


nc.pg.run(mode="Training")
nc.pg.run(mode="Inference")
mutation.table<-table.generator(df="/Training/Training_dataset.csv")

old.models<-list.files("/Models/", full.names = FALSE)
old.models<-old.models[-grep("legacy_models", old.models)]
if(length(old.models)>0){
  file.rename(paste("/Models/",old.models,sep = ""), paste("/Models/legacy_models/",old.models,sep = ""))
}

write.csv(mutation.table, paste("/Models/", gsub("-","",Sys.Date()),"_MutationTable.csv",sep = ""), row.names = FALSE)
data<-read.csv("/Training/Training_dataset.csv")
training.set<-trainset.prep(data = data)
ml.set<-P.calculator(input.data = training.set, mutation.table = mutation.table)
model<-cnn.training(training.set = ml.set)

inference.raw<-read.csv("/Inference/Inference_dataset.csv")
resultsInference<-P.calculator(input.data = inference.raw, mutation.table = mutation.table)

results.out<-ml.inference(inference.set = resultsInference, model = model )

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



