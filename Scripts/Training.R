source("/home/docker/Scripts/Functions.R")

nc.pg.run(mode="Training")
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
cnn.training(training.set = ml.set)
