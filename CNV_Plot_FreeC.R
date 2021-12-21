## Loading Library
library(ggplot2)
library(reshape2)
library(matrixStats)

### Loading Data
args <- commandArgs(TRUE)
data <- read.table(args[1],header=T)
head(data)
data <- data[,-2]
sample <- read.table(args[2],header=F)
filename <- args[3]
#outdir<-"/WPSnew/wangrui/Project/Human_Eight_Cell/PBAT/StatInfo/01.CNV/FreeC_Result/"
outdir<- args[4]
## Function Definition

do_CNV_Plot_FreeC <- function(x,sample=NULL){
  x <- as.data.frame(x)
  if(is.null(sample) == FALSE){
    colnames(x) <- c("Chr",as.character(sample[,1]))
  }
  x <- x[!x$Chr == "M",] #
  x <- x[!x$Chr == "Y",] #not includ Y chromosome
  x$Chr <- factor(x$Chr,levels = c(1:22,"X"),ordered=T)
  x <- x[order(factor(x$Chr,levels = c(1:22,"X"),ordered=T)),]
  x$Pos <- 1:nrow(x)

  chr_window_num <- as.data.frame(table(x$Chr))
  colnames(chr_window_num)<-c("Chr","Number")
  chr_window_num$position <-rep(0,23)
  for(i in 2:23){
    chr_window_num[i,2] <-chr_window_num[i,2]+chr_window_num[i-1,2]
  }
  chr_window_num[1,3] <- floor(chr_window_num[1,2]/2)
  for(i in 2:23){
    chr_window_num[i,3]<-floor(chr_window_num[i-1,2]+0.5*(chr_window_num[i,2]-chr_window_num[i-1,2]))
  }
  melt_data <- melt(x,id.vars = c("Chr","Pos"))
  colnames(melt_data)<-c("Chr","Pos","Sample","Copy_number")
  melt_data$Colors <-2
  melt_data[melt_data$Chr%in%c(seq(1,22,2),"X"),"Colors"] <-1
  melt_data$Colors <- factor(melt_data$Colors)
 P1<- ggplot(data=melt_data,aes(Pos,Copy_number,colour=Colors))+geom_point(size=2)+facet_grid(Sample ~.)+geom_vline(xintercept = chr_window_num[1:23,2],colour="black", linetype = "longdash")+geom_hline(yintercept=c(0.5,1,1.5,2),alpha=0.5,colour="grey")+scale_color_manual(values=c("blue", "red"))+scale_x_discrete(breaks=chr_window_num$position,labels=chr_window_num$Chr)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_text(angle=90,hjust=1))+ylim(0,2.3)

 print(P1)

  melt_data$CNV <- "Normal"
  melt_data[melt_data$Copy_number>1.3,"CNV"] <- "Gain"
  melt_data[melt_data$Copy_number<0.7,"CNV"] <- "Loss"
  melt_data$CNV <- factor(melt_data$CNV,levels=c("Normal","Gain","Loss"),ordered=T)
 P2<-  ggplot(data=melt_data,aes(Pos,Copy_number,colour=CNV))+geom_point(size=2)+facet_grid(Sample ~.)+geom_vline(xintercept = chr_window_num[1:23,2],colour="black", linetype = "longdash")+geom_hline(yintercept=c(0.5,1,1.5,2),alpha=0.5,colour="grey")+scale_color_manual(values=c("lightgrey","red","blue")) +scale_x_discrete(breaks=chr_window_num$position,labels=chr_window_num$Chr)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_text(angle=90,hjust=1))+ylim(0,2.3)

 print(P2)
}

pdf(file=paste(outdir,filename,sep="/"),width=17,height=4+nrow(sample)*1.5)
do_CNV_Plot_FreeC(data,sample=sample)
dev.off()
