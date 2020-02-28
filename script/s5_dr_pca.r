
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(Matrix))
suppressMessages(library(methods))
suppressMessages(library(scales))
suppressMessages(library(ggplot2))




args=commandArgs(T);
#sample=as.character(args[1]);
sample=args[1];


setwd("./")
path_out="./"

data1=read.table(paste("data_ref/expr.txt",sep="")) 
data2=read.table(paste("data_obj/",sample,".txt",sep="")) 
gene1=read.table(paste("data_ref/gene.txt",sep=""))


data3=rbind(data1[,1:ncol(data1)],data2[,1:ncol(data2)])
data4=t(data3)
cell_all=nrow(data3)
colnames(data4) <- paste("cell",c(1:cell_all),sep="")
rownames(data4) <- gene1[,2]



suppressWarnings(suppressMessages(bmmc <- CreateSeuratObject(counts = data4, project = "test")))
suppressWarnings(suppressMessages(bmmc <- ScaleData(object = bmmc)))

suppressWarnings(suppressMessages(bmmc <- RunPCA(bmmc, features = rownames(bmmc), npcs = 30,weight.by.var = TRUE)))


write.table(bmmc@reductions$pca@cell.embeddings[,1:30],file=paste(path_out,"temp_out/s5_embedding_",sample,".txt",sep=""),sep = "\t",row.names=F,col.names=F)


