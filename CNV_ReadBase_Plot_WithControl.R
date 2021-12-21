##Loading Library
library(ggplot2)
library(reshape2)
library(matrixStats)

##Loading Data
args <- commandArgs(TRUE)
data <- read.table(args[1],header=T)
sample <- read.table(args[2],header=F)
sample <- as.character(sample[,1])
filename <- args[3]
outdir <- as.character(args[4])

control <- c(grep("-N",sample,value=T),grep("_N",sample,value=T))

## Set control sample
#control_sample <- c("CAT_a","CAT_b","CAT_c")

print("Loading Data")
##Function Definition
do_CNV_Plot_Reads <- function(x,control=NULL,sample=NULL,remove=NULL,return_file=FALSE){
  x <- as.data.frame(x)
  if(is.null(sample) == FALSE){
    colnames(x) <- c("Chr","window",sample)
  }
  if(is.null(remove)==FALSE){
    x <- x[,!colnames(x)%in%remove]
  }

  x <- x[!x$Chr == "M",] #
  # Normalize data by seqence number
  x[,3:ncol(x)]<-t(t(x[,3:ncol(x)])/colSums(x[,3:ncol(x)]))
  #x <- x[!x$Chr == "Y",] #not includ Y chromosome
  x$Chr <- factor(x$Chr,levels = c(1:22,"X","Y"),ordered=T)
  x <- x[order(factor(x$Chr,levels = c(1:22,"X","Y"),ordered=T)),]

  if(is.null(control) == TRUE){
    control = x[,3:ncol(x)]
  }else{
#    control <- control[!control[,1] == "M",]
#    control[,1] <- factor(control[,1],levels = c(1:22,"X","Y"))
#    control <-control[order(factor(control[,1],levels = c(1:22,"X","Y"),ordered=T)),]
#    control = control[,3:ncol(control)]
         control = x[,colnames(x)%in%control]
  }
  #Normalize among sample by rowMeadians
  x[,3:ncol(x)]<- as.data.frame(x[,3:ncol(x)]/rowMedians(as.matrix(control)))

  x[is.na(x)]<- -1
  ## Set Big number to 2.5
  y <- x[,3:ncol(x)]
  y[y>2.5] <-2.5
  x[,3:ncol(x)]  <- y
#  write.table(file=paste(outdir,gsub(".pdf",".txt",filename),sep=""),quote=F)

  x$Pos <- 1:nrow(x)

  chr_window_num <- as.data.frame(table(x$Chr))
  colnames(chr_window_num)<-c("Chr","Number")
  chr_window_num$position <-rep(0,24)
  for(i in 2:24){
    chr_window_num[i,2] <-chr_window_num[i,2]+chr_window_num[i-1,2]
  }
  chr_window_num[1,3] <- floor(chr_window_num[1,2]/2)
  for(i in 2:24){
    chr_window_num[i,3]<-floor(chr_window_num[i-1,2]+0.5*(chr_window_num[i,2]-chr_window_num[i-1,2]))
  }
  melt_data <- melt(x,id.vars = c("Chr","window","Pos"))
  colnames(melt_data)<-c("Chr","window","Pos","Sample","Copy_number")
  melt_data$Colors <-2
  melt_data[melt_data$Chr%in%c(seq(1,22,2),"X"),"Colors"] <-1
  melt_data$Colors <- factor(melt_data$Colors)
  P1 <- ggplot(data=melt_data,aes(Pos,Copy_number,colour=Colors))+geom_point(size=1.5)+facet_grid(Sample ~.)+geom_vline(xintercept = chr_window_num[1:23,2],colour="black", linetype = "longdash")+geom_hline(yintercept=c(0.5,1,1.5,2,2.5),alpha=0.5,colour="grey")+geom_hline(yintercept=1,alpha=0.4,colour="darkgrey",size=0.9)+scale_color_manual(values=c("blue", "red"))+scale_x_discrete(breaks=chr_window_num$position,labels=chr_window_num$Chr)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_text(angle=90,hjust=1))+ylim(0,2.3)
 print(P1)

 P3 <- ggplot(data=melt_data,aes(Pos,2*Copy_number,colour=Colors))+geom_point(size=1.5)+facet_grid(Sample ~.)+geom_vline(xintercept = chr_window_num[1:23,2],colour="black", linetype = "longdash")+geom_hline(yintercept=c(1,2,3,4),alpha=0.5,colour="grey")+geom_hline(yintercept=1,alpha=0.4,colour="darkgrey",size=0.9)+scale_color_manual(values=c("blue", "red"))+scale_x_discrete(breaks=chr_window_num$position,labels=chr_window_num$Chr)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_text(angle=90,hjust=1))+ylim(0,4.5)
 print(P3)

# return(x)
#       ggplot(data=melt_data,aes(Pos,Copy_number,colour=Colors))+geom_point()+facet_grid(Sample ~.)+geom_vline(xintercept = chr_window_num[1:23,2],colour="black", linetype = "longdash")+geom_hline(yintercept=c(0.5,1,1.5,2,2.5),alpha=0.5,colour="grey")+geom_hline(yintercept=1,alpha=0.4,colour="darkgrey",size=0.9)+scale_color_manual(values=c("blue", "red"))+scale_x_discrete(breaks=chr_window_num$position,labels=chr_window_num$Chr)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_text(angle=90,hjust=1))+ylim(0,2.7)

 melt_data$CNV <- "Normal"
 melt_data[melt_data$Copy_number>1.3,"CNV"] <- "Gain"
 melt_data[melt_data$Copy_number<0.7,"CNV"] <- "Loss"
 melt_data$CNV <- factor(melt_data$CNV,levels=c("Normal","Gain","Loss"),ordered=T)
 P2 <- ggplot(data=melt_data,aes(Pos,Copy_number,colour=CNV))+geom_point(size=1.5)+facet_grid(Sample ~.)+geom_vline(xintercept = chr_window_num[1:23,2],colour="black", linetype = "longdash")+geom_hline(yintercept=c(0.5,1,1.5,2,2.5),alpha=0.5,colour="grey")+geom_hline(yintercept=1,alpha=0.4,colour="darkgrey",size=0.9)+scale_color_manual(values=c("lightgrey","red","blue"))+scale_x_discrete(breaks=chr_window_num$position,labels=chr_window_num$Chr)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_text(angle=90,hjust=1))+ylim(0,2.3)
 print(P2)
 if(return_file) {return(x)}
}

#nsample <- length(sample)-length(remove)
nsample <- length(sample)
pdf(file=paste(outdir,filename,sep="/"),width=22,height=4+nsample*1.5)
do_CNV_Plot_Reads(data,control=control,sample=sample,remove=NULL)
#tmp_matrix_Median<-do_CNV_Plot_Reads(data,control=NULL,sample=sample,remove=remove)
dev.off()

#pdf(file=paste(outdir,paste("WithControl",filename,sep="_"),sep=""),width=17,height=4+nsample*0.7)
tmp_matrix_Control<-do_CNV_Plot_Reads(data,control=control,sample=sample,return_file=TRUE)
#do_CNV_Plot_Reads(data,control=control_sample,sample=sample,remove=NULL)
#dev.off()
write.table(file=paste(outdir,paste("WithControl",gsub(".pdf",".txt",filename),sep="_"),sep=""),tmp_matrix_Control,quote=F)



#CNV_Summary_Matrix <- as.data.frame(matrix(rep(0,(ncol(tmp_matrix_Median)-2)*24),ncol=24))
#colnames(CNV_Summary_Matrix) <- paste("chr",c(1:22,"X","Y"),sep="")
#rownames(CNV_Summary_Matrix)<-colnames(tmp_matrix_Median)[-c(1:2)]
#head(CNV_Summary_Matrix)

chromosome <- c(1:22,"X","Y")
stat <-function(x){
        CNV_Summary_Matrix <- as.data.frame(matrix(rep(0,(ncol(x)-3)*24),ncol=24))
        colnames(CNV_Summary_Matrix) <- paste("chr",c(1:22,"X","Y"),sep="")
        rownames(CNV_Summary_Matrix)<-colnames(x)[-c(1:2,ncol(x))]
        for (i in 1:length(chromosome)){
                chr <-chromosome[i]
                sub_data<-x[x$Chr==chr,]
                sub_data<- sub_data[,-ncol(sub_data)]
                sub_data <- sub_data[,-c(1:2)]
                gain <- apply(sub_data,2,function(x){ratio<-length(x[x>1.3])/length(x);return(ratio)})
                lost <- apply(sub_data,2,function(x){ratio<-length(x[x<0.7 & x>=0])/length(x);return(ratio)})
                tmp <- rep(0,length(gain))
                for(k in 1:length(gain)){
                if(i<23){
                        if(gain[k]>0.1)tmp[k]=1
                        else if(lost[k]>0.1)tmp[k]=-1
                }else{
                        if(gain[k]>0.50)tmp[k]=1
            else if(lost[k]>0.50)tmp[k]=-1
                }
                }
        CNV_Summary_Matrix[,i]<-tmp
        }
        return(CNV_Summary_Matrix)
}

CNV_Summary_Matrix<-stat(tmp_matrix_Control)
#Sex_Summary_Matrix <- CNV_Summary_Matrix[,c(23,24)]
#Autosome_Summary_Matrix <-CNV_Summary_Matrix[,1:22]
#CNV_Summary_Matrix_Female <- stat(tmp_matrix_Female)
#CNV_Summary_Matrix_Male <- stat(tmp_matrix_Male)
write.table(file=paste(outdir,paste("CNV_Summary",gsub(".pdf",".txt",filename),sep="_"),sep="/"),CNV_Summary_Matrix,quote=F)
library("ComplexHeatmap")
library("methods")
pdf(paste(outdir,"/","Summary_",filename,sep=""))
Heatmap(CNV_Summary_Matrix[,1:22],name="CNV",column_title = "Chromosome(Human Eye)",column_title_gp = gpar(fontsize = 14, fontface = "bold"),row_title="Samples",row_title_gp = gpar(fontsize = 14, fontface = "bold"),cluster_rows=FALSE,cluster_columns=FALSE,col=c("1"="red","-1"="blue","0"="lightgrey"),rect_gp = gpar(col= "white"))

#Heatmap(CNV_Summary_Matrix[,1:22],name="CNV",column_title = "Autosome(Normal_Control)",column_title_gp = gpar(fontsize = 14, fontface = "bold"),row_title="Samples",row_title_gp = gpar(fontsize = 14, fontface = "bold"),cluster_rows=FALSE,cluster_columns=FALSE,col=c("1"="red","-1"="blue","0"="lightgrey"),rect_gp = gpar(col= "white")) + Heatmap(CNV_Summary_Matrix[,23:24],name="Sex",column_title = "Sex_C",show_row_names = FALSE,column_title_gp = gpar(fontsize = 7),cluster_rows=FALSE,cluster_columns=FALSE,col=c("1"="red","-1"="blue","0"="lightgrey"),rect_gp = gpar(col= "white")) #+ Heatmap(CNV_Summary_Matrix_Female[,23:24],name="Sex",column_title = "Sex_F",column_title_gp = gpar(fontsize = 7),cluster_rows=FALSE,cluster_columns=FALSE,col=c("1"="red","-1"="blue","0"="lightgrey"),rect_gp = gpar(col= "white"),show_row_names = FALSE) +Heatmap(CNV_Summary_Matrix_Male[,23:24],name="M",column_title = "Sex_M",column_title_gp = gpar(fontsize = 7),cluster_rows=FALSE,cluster_columns=FALSE,col=c("1"="red","-1"="blue","0"="lightgrey"),rect_gp = gpar(col= "white"))
dev.off()
