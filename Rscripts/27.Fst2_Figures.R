library(ggplot2)
library(ggpubr)
library(ggthemes)
library(plotrix)
library(grid)
library(purrr)
library(zoo)
library(tidyverse)
library(gplots)

source("Rscripts/baseRscript.R")

# read the files saved in Overview_output:
FstFiles<-list.files("Output_all/Diversity/",pattern="^Fst2.*.csv")

FstList<-list()
for (i in 1:length(FstFiles)){ 
        File<-read.csv(paste0("Output_all/Diversity/",FstFiles[i]),stringsAsFactors=FALSE)
        File<-File[,-1]
        FstList[[i]]<-File
        names(FstList)[i]<-substr(paste(FstFiles[i]),start=1,stop=7)
}

##############
#split the combination (xx vs. yy) names to each column
FstList2<-list()
for (j in 1:length(FstList)){
        data<-FstList[[j]]
        splitted<-strsplit(data$comb, split=" ")
        for (i in 1:nrow(data)){
                data$Sample1[i]<-splitted[[i]][1]
                data$Sample2[i]<-splitted[[i]][2]
        }
        FstList2[[j]]<-data
        names(FstList2)[j]<-names(FstList)[j]
}



geno<-c("1A","1B","3A")


colorRampBlue <- colorRampPalette(c("white", "steelblue1", "blue3"))

for (i in 1:length(FstList2)){
        dat<-FstList2[[i]]
        dat<-dat[,4:6]
        #for (j in 1:nrow(dat)){
        #        dat$Sample.1[j]<-paste0(substr(paste(dat$Sample1[j]), start=2,stop=2),".", samples$Sample2[samples$SampleID==dat$Sample1[j]])
        #        dat$Sample.2[j]<-paste0(substr(paste(dat$Sample2[j]), start=2,stop=2),".",samples$Sample2[samples$SampleID==dat$Sample2[j]])
        #        
        #}
        #s<-unique(c(dat$Sample.1,dat$Sample.2))
        #for (j in 1:length(s)){
        #        newrow<-c(0,"","",s[j],s[j])
        #        dat<-rbind(dat,newrow)
        #}
        #dat$Fst<-as.numeric(dat$Fst)
        dat$Fst<-as.numeric(format(round(dat$Fst,4),  nsmall=4))
        dat2<-dat[!is.na(dat$Fst),]
        dat2<-dat2[dat2$Sample1!="D75030"&dat2$Sample2!="D75030",]
        title<-paste0("Genotype ",substr(names(FstList2)[i], start=6, stop =7 ))
        p1<-ggplot(data=dat2, aes(Sample1,Sample2, fill=Fst)) +geom_tile(color="gray40")+ 
                scale_fill_gradientn(colors =colorRampBlue(62), limits=c(0,1))+
                theme_light()+
                theme(axis.title.x = element_blank(), axis.title.y = element_blank(),panel.grid.major = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      axis.text.x=element_text(size=9, angle=45,hjust=1))+ggtitle(title)
       
        ggsave(file=paste0("Output_all/Diversity/Fst_",names(FstList2)[i],".Plot.pdf"),plot=p1,width=9, height=9, units='in',device='pdf')
}




#### 
## with 1B
library(gplots)
library(heatmap3)

mat<-with(dat,tapply(Fst,list(Sample2,Sample1),"[[",1)) 
#heatmap.2(mat, col=colorRampPalette(c("white","blue"))(50) )
#heatmap3(mat, col=colorRampPalette(c("white","blue"))(50) )

name1<-dat$Sample1[1]
name2<-dat$Sample2[nrow(dat)]
        
mat2<-rbind(name1=rep(NA,times=nrow(mat)),mat)
mat3<-cbind(mat2,name2=rep(NA,times=nrow(mat2)))
diag(mat3)<-0
rownames(mat3)[1]<-name1


colorRampBlue <- colorRampPalette(c("white", "steelblue1", "blue3"))
mat3m<-melt(mat3, na.rm=T)
mat3m<-mat3m[!is.na(mat3m$value),]
ggplot(data=mat3m, aes(Var1,Var2, fill=value)) +geom_tile(color="white")+ 
        scale_fill_gradientn(colors=colorRampBlue(64),limit=c(0,1))+
        theme_minimal()

