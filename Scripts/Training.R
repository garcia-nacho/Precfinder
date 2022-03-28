source("/home/docker/Scripts/Functions.R")

args=commandArgs(TRUE)

nc.pg.run(mode="Training")
mutation.table<-table.generator(df="/Training/Training_dataset.csv")

write.csv(mutation.table, paste("/Models/", args[1], "/", gsub("-","",Sys.Date()),"_MutationTable.csv",sep = ""), row.names = FALSE)
data<-read.csv("/Training/Training_dataset.csv")
training.set<-trainset.prep(data = data)
ml.set<-P.calculator(input.data = training.set, mutation.table = mutation.table)
cnn.training(training.set = ml.set, modelid=args[1])
