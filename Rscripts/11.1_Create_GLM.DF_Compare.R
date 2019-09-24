#GLM/Beta regression preparation for Q35, Q35_NA, Q35_zero:

library(tidyverse)
library(zoo)
library(purrr)
library(DataCombine)
library(miscTools)


source("Rscripts/baseRscript.R")

## Run 10.2.Filterout_lowfrew.R first ##
# Read the summary data files
filtered<-list()
fnames<-c("Ts", "Ts_NA", "Ts_zero")
for (i in 1:3){
        filename<-fnames[i]
        df<-read.csv(paste0("Output1A/Q35Compare/Summary_",filename,".Q35.csv"),stringsAsFactors = F)
        df<-df[,-1]
        filtered[[i]]<-df
        names(filtered)[i]<-filename
}

#### For preparing GLM data formatting

nucord <- c("a", "t", "c", "g")
glmform<-as.matrix(data.frame(a="",t="",c="",g="", Syn ="",Nonsyn="",Stop=""))
nuc<-as.matrix(data.frame(a="",t="",c="",g=""))
for (k in 1:3){
        dat<-filtered[[k]]
        dat1<-dat[,c("pos","makesCpG","bigAAChange","mean")]
        
        for (i in 1:nrow(dat)){
                atcg <- c(0,0,0,0)
                atcg[which(nucord == dat[i,]$ref)] <- 1
                nuc<-miscTools::insertRow(nuc,i,atcg)
                nonsyn <- as.numeric(regexpr("nonsyn",dat[i,]$Type) > 0)
                stop <- as.numeric(regexpr("stop",dat[i,]$Type) > 0)
                syn<-as.numeric(regexpr("^syn",dat[i,]$Type) > 0)
                new<-c(atcg,syn,nonsyn,stop)
                glmform<-insertRow(glmform,i,new)
        }
        GlmData<-cbind(dat1$pos,glmform[1:nrow(dat1),])
        GlmData<-cbind(GlmData,dat1[,2:4])
        colnames(GlmData)[9]<-"CpG"
        write.csv(GlmData, paste0("Output1A/GLM/GlmData_",names(filtered)[k],".csv"))
}
