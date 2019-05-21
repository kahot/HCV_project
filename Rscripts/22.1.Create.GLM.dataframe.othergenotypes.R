#GLM/Beta regression preparation:

library(tidyverse)
library(zoo)
library(purrr)
library(miscTools)


source("Rscripts/baseRscript.R")

# Read the summary data for different 
filtered<-list()
fnames<-c("1B","3A")
for (i in 1:length(fnames)){
        filename<-fnames[i]
        df<-read.csv(paste0("Output_all/Ts_summary_metadata.",filename,".csv"),stringsAsFactors = F)
        df<-df[,-1]
        filtered[[i]]<-df
        names(filtered)[i]<-filename
}


#### For preparing GLM data formatting

nucord <- c("a", "t", "c", "g")
glmform<-as.matrix(data.frame(a="",t="",c="",g="", Syn ="",Nonsyn="",Stop=""))
nuc<-as.matrix(data.frame(a="",t="",c="",g=""))
for (k in 1:length(filtered)){
        dat<-filtered[[k]]
        name<-fnames[k]
        dat1<-dat[,c("merged.pos",paste0("makesCpG.",name),paste0("bigAAChange.",name),"mean")]
        
        for (i in 1:nrow(dat)){
                atcg <- c(0,0,0,0)
                atcg[which(nucord == dat[i,paste0("nuc.",name)])] <- 1
                nuc<-insertRow(nuc,i,atcg)
                nonsyn <- as.numeric(regexpr("nonsyn",dat[i,paste0("Type.",name)]) > 0)
                stop <- as.numeric(regexpr("stop",dat[i,paste0("Type.",name)]) > 0)
                syn<-as.numeric(regexpr("^syn",dat[i,paste0("Type.",name)]) > 0)
                new<-c(atcg,syn,nonsyn,stop)
                glmform<-insertRow(glmform,i,new)
        }
        GlmData<-cbind(dat1$merged.pos,glmform[1:nrow(dat1),])
        GlmData<-cbind(GlmData,dat1[,2:4])
        colnames(GlmData)[1]<-"pos"
        colnames(GlmData)[9]<-"CpG"
        colnames(GlmData)[10]<-"bigAAChange"
        write.csv(GlmData, paste0("Output_all/GlmData_",name,".csv"))
}
        
