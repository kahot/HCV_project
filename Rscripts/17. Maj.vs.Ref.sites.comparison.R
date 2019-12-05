library(reshape)
library(tidyverse)
library(zoo)
library(purrr)

source("Rscripts/baseRscript.R")


HCVFiles_overview2<-list.files("Output1A/Overview2/",pattern="overview2.csv")
MajRef<-data.frame(sample=substr(paste(HCVFiles_overview2),start=1,stop=6))

MajequalRef1<-list()
MajnotRef1<-list()
MajequalRef2<-list() #remove mut freq>0.2
MajnotRef2<-list()   #remove mut freq>0.2
for (i in 1:length(HCVFiles_overview2)){ 
        DF<-read.csv(paste0("Output1A/Overview2/",HCVFiles_overview2[i]),stringsAsFactors=FALSE, row.names = 1)
        sname<-substr(paste(HCVFiles_overview2[i]),start=1,stop=6)
        
        low_reads<-which(DF$TotalReads<1000) 
        DF[low_reads,c(8:42)]<-NA
        DF<-DF[!is.na(DF$freq.Ts),]
        
        mf_high<-which(DF$freq.Ts>=0.2|DF$freq.transv>=0.2) 
        DF2<-DF
        DF2[mf_high,c(8:19)]<-NA
        DF2<-DF2[!is.na(DF2$freq.Ts),]
        
        df1<-DF[DF$MajNt==DF$ref,]
        df2<-DF[DF$MajNt!=DF$ref,]
        
        df11<-DF2[DF2$MajNt==DF2$ref,]
        df22<-DF2[DF2$MajNt!=DF2$ref,]
        
        MajRef$Same1[i]<-nrow(df1)
        MajRef$Different1[i]<-nrow(df2)
        MajRef$Same2[i]<-nrow(df11)
        MajRef$Different2[i]<-nrow(df22)
        
        MajequalRef1[[i]]<-df1
        MajnotRef1[[i]]<-df2
        names(MajequalRef1)[i]<-sname
        names(MajnotRef1)[i]<-sname
        MajequalRef2[[i]]<-df11
        MajnotRef2[[i]]<-df22
        names(MajequalRef2)[i]<-sname
        names(MajnotRef2)[i]<-sname
}


MajRef$Percet.Diff1<-MajRef$Different1/(MajRef$Same1+MajRef$Different1)
MajRef$Percet.Diff2<-MajRef$Different2/(MajRef$Same2+MajRef$Different2)

mean(MajRef$Percet.Diff1)
#0.06643408
std.error(MajRef$Percet.Diff1)
# 0.001166803
mean(MajRef$Percet.Diff2)
#0.06275597
std.error(MajRef$Percet.Diff2)
#0.00122738


## Calculate mutation frequncies of maj==ref vs. maj!=ref
mf1<-data.frame(sample=substr(paste(HCVFiles_overview2),start=1,stop=6))
mf2<-data.frame(sample=substr(paste(HCVFiles_overview2),start=1,stop=6))

for (i in 1:length(MajequalRef1)){
        df1<-MajequalRef1[[i]]
        df2<-MajnotRef1[[i]]
        mf1$same.mean[i]<-mean(df1$freq.Ts, na.rm = T)
        mf1$diff.mean[i]<-mean(df2$freq.Ts, na.rm = T)
        
        df12<-MajequalRef2[[i]]
        df22<-MajnotRef2[[i]]
        mf2$same.mean[i]<-mean(df12$freq.Ts, na.rm = T)
        mf2$diff.mean[i]<-mean(df22$freq.Ts, na.rm = T)
}

mean(mf1$same.mean)  #0.0062938
mean(mf1$diff.mean, na.rm = T)  #0.03312996

mean(mf2$same.mean)  #0.004576401
mean(mf2$diff.mean, na.rm = T)  #0.01505648


result1<-wilcox.test(mf1$same.mean,mf1$diff.mean  , alternative = "less", paired = FALSE) 
result[3]
#p.value
# 4.266554e-48

result2<-wilcox.test(mf2$same.mean,mf2$diff.mean  , alternative = "less", paired = FALSE) 
result2[3]
#$p.value
#[1] 8.486596e-50


#where did they differ?

positions1<-c()
positions2<-c()

for (i in 1:length(MajequalRef1)){
        df2<-MajnotRef1[[i]]
        positions1<-c(positions1, df2[,"pos"])
        df22<-MajnotRef2[[i]]
        positions2<-c(positions2, df22[,"pos"])
 }

hist(positions1)
hist(positions2)
pos1<-as.data.frame(table(positions1))
pos2<-as.data.frame(table(positions2))
plot(pos1)
plot(pos2)

pos1$positions1<-as.numeric(as.character(pos1$positions1))
pos1.st<-pos1[pos1>341& pos1<=2579,]
pos1.nonst<-pos1[pos1>2579,]

sum(pos1.st$Freq)/2238 #12.25469
sum(pos1.nonst$Freq)/6070 #9.449259



#more structural regions than non-structural genes 