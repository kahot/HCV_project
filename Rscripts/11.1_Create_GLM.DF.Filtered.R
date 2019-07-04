#GLM/Beta regression preparation:
# RUnnig only the filtered Q35 Ts data (no Ts_zero and Ts-NA)

library(tidyverse)
library(zoo)
library(purrr)
library(DataCombine)
library(miscTools)


source("Rscripts/baseRscript.R")

## Run 10.2.Filterout_lowfrew.R first ##

#### For preparing GLM data formatting

nucord <- c("a", "t", "c", "g")
glmform<-as.matrix(data.frame(a="",t="",c="",g="", Syn ="",Nonsyn="",Stop=""))
nuc<-as.matrix(data.frame(a="",t="",c="",g=""))

dat<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F, row.names=1)
dat<-dat[dat$pos>=342,]
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
write.csv(GlmData, "Output1A/GLM/GlmData_Ts_FilteredData.csv")

