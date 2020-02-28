library(ggplot2)
library(ggforce)
library(Matrix)
library(methods)
library(scales) 
#library(igraph)
library(RColorBrewer)
#options(stringsAsFactors=FALSE)

args=commandArgs(T);
sample=as.character(args[1]);
workingDir = "./";
setwd(workingDir);


score1 <- read.table(paste("temp_out/s9_",sample,"_score.txt",sep=""))
color_list<-hue_pal()(2)

df1<-data.frame(score=score1[,3],coordx=score1[,1],coordy=score1[,2])
df1<-df1[with(df1,order(score)),]

pl<-ggplot()+
geom_point(data=df1,aes(x=coordx,y=coordy,colour=score),size=2,alpha=0.7)+
annotate(geom = 'text', label = sample, x = -Inf, y = Inf, hjust = -0, vjust = 1.1,size=9)+
#guides(fill=guide_legend(title="abundance"),title.position = "center")+
#guides(colour=guide_legend(title="abundance"),title.position = "center")+
#guides(colour=FALSE)+

labs(x="Component 1",y="Component 2")+theme_bw()+


theme(

        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
	legend.text=element_text(size=12), 
	legend.title=element_text(size=12), 
	#legend.position="top",
       # legend.position="none",

        # plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank())
#pl<-pl+scale_colour_gradient(name="score",low = "grey", high = "red")
#pl<-pl+scale_colour_gradientn(limits = c(0,1),colours=c("lightblue","blue","red"),breaks=c(0.2,0.4,0.6,0.8))
pl<-pl+scale_colour_gradientn(name="abundance",limits = c(0,1),colours=c("grey","red"),breaks=c(0.2,0.4,0.6,0.8))
pl<-pl+scale_fill_discrete(guide="none")


pdf(paste("result/projection_",sample,".pdf",sep=""),width=6,height=5)
pl
dev.off()

