library(ggplot2)
library(Matrix)
library(methods)
library(scales)


#args=commandArgs(T);
#sample=as.character(args[1]);


color_list<-hue_pal()(27)


setwd("./")
path_out="result/"


ident_list=as.matrix(read.table("data_skeleton/ident_list_coords.txt",sep="\t"))
ident=as.matrix(read.table("data_skeleton/ident.txt"))
coords=as.matrix(read.table("data_skeleton/coords.txt"))


df1<-data.frame(coordx=coords[,1],coordy=coords[,2],group=ident[,1])
df1$group<-factor(df1$group,levels=c(unique(ident_list[,1])), ordered=TRUE)
df2<-data.frame(coordx=as.numeric(ident_list[,2]),coordy=as.numeric(ident_list[,3]),label=ident_list[,1])
#df2$group<-factor(df2$group,levels=c(unique(ident_list[,1])), ordered=TRUE)

pdf(paste(path_out,"lineage_skeleton_ref.pdf",sep=""),width=6,height=6)
par(mai=c(0,0,0,0));
pl<-ggplot()+
geom_point(data=df1,aes(x=coordx,y=coordy,colour=group),size=0.5,alpha=1)+
geom_text(data=df2,aes(x=coordx,y=coordy,label=label),size=4,colour="black")+
    labs(x="tSNE1",y="tSNE2")+theme_bw()+
theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=16),
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.title=element_text(size=16),
        legend.position = "none",
        panel.background = element_blank()
)
pl<-pl+scale_colour_manual(values=color_list)
pl
dev.off()





