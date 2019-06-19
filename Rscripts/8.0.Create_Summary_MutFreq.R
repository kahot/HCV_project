library(purrr)
library(tidyverse)
library(zoo)
library(plotrix)
library(sfsmisc)

#Prep data
source("Rscripts/baseRscript.R")

HCVFiles_overview<-list.files("Output1A/Overview1/",pattern="overview.csv")

Overview_summary<-list()
for (i in 1:length(HCVFiles_overview)){ 
        overviews<-read.csv(paste0("Output1A/Overview1/",HCVFiles_overview[i]),stringsAsFactors=FALSE)
        overviews<-overviews[,-1]
        Overview_summary[[i]]<-overviews
        names(Overview_summary)[i]<-substr(paste(HCVFiles_overview[i]),start=1,stop=7)
}

## 

#source("Rscripts/MutationFreqSum.R")
MutFreq_Ts<-list()
MutFreq_tv1<-list()
MutFreq_tv2<-list()
MutFreq_tvs<-list()
MutFreq_all<-list()

MutFreq_Ts.ref<-list()
MutFreq_tv1.ref<-list()
MutFreq_tv2.ref<-list()
MutFreq_tvs.ref<-list()
MutFreq_all.ref<-list()
        
        
for (i in 1:length(Overview_summary)){
        dat<-Overview_summary[[i]]
        filename<-names(Overview_summary)[i]
        
        MutFreq_Ts[[i]]<-dat[,c("pos","freq.Ts")] 
        MutFreq_tv1[[i]]<-dat[,c("pos","freq.transv1")] 
        MutFreq_tv2[[i]]<-dat[,c("pos","freq.transv2")] 
        MutFreq_tvs[[i]]<-dat[,c("pos","freq.transv")] 
        MutFreq_all[[i]]<-dat[,c("pos","freq.mutations")]
        
        names(MutFreq_Ts)[i]<-filename
        names(MutFreq_tv1)[i]<-filename
        names(MutFreq_tv2)[i]<-filename
        names(MutFreq_tvs)[i]<-filename
        names(MutFreq_all)[i]<-filename
        
         MutFreq_Ts.ref[[i]]<-dat[,c("pos","freq.Ts.ref")] 
        MutFreq_tv1.ref[[i]]<-dat[,c("pos","freq.transv1.ref")] 
        MutFreq_tv2.ref[[i]]<-dat[,c("pos","freq.transv2.ref")] 
        MutFreq_tvs.ref[[i]]<-dat[,c("pos","freq.transv.ref")] 
        MutFreq_all.ref[[i]]<-dat[,c("pos","freq.mutations.ref")]
        
         names(MutFreq_Ts.ref)[i]<-filename
        names(MutFreq_tv1.ref)[i]<-filename
        names(MutFreq_tv2.ref)[i]<-filename
        names(MutFreq_tvs.ref)[i]<-filename
        names(MutFreq_all.ref)[i]<-filename
}        
        
#assign column names for the list
for (i in 1:length(MutFreq_Ts)) {
        colnames(MutFreq_Ts[[i]])<-c("pos",paste0(names(MutFreq_Ts[i])))
        colnames(MutFreq_tv1[[i]])<-c("pos",paste0(names(MutFreq_tv1[i])))
        colnames(MutFreq_tv2[[i]])<-c("pos",paste0(names(MutFreq_tv2[i])))
        colnames(MutFreq_tvs[[i]])<-c("pos",paste0(names(MutFreq_tvs[i])))
        colnames(MutFreq_all[[i]])<-c("pos",paste0(names(MutFreq_all[i])))
        
         colnames(MutFreq_Ts.ref[[i]])<-c("pos", paste0(names(MutFreq_Ts.ref[i])))
        colnames(MutFreq_tv1.ref[[i]])<-c("pos",paste0(names(MutFreq_tv1.ref[i])))
        colnames(MutFreq_tv2.ref[[i]])<-c("pos",paste0(names(MutFreq_tv2.ref[i])))
        colnames(MutFreq_tvs.ref[[i]])<-c("pos",paste0(names(MutFreq_tvs.ref[i])))
        colnames(MutFreq_all.ref[[i]])<-c("pos",paste0(names(MutFreq_all.ref[i])))
}

TsMutFreq<-MutFreq_Ts%>% purrr::reduce(full_join, by='pos')
Tv1.MutFreq<-MutFreq_tv1 %>% purrr::reduce(full_join, by='pos')
Tv2.MutFreq<-MutFreq_tv2 %>% purrr::reduce(full_join, by='pos')
Tvs.MutFreq<-MutFreq_tvs %>% purrr::reduce(full_join, by='pos')
AllMutFreq<-MutFreq_all %>% purrr::reduce(full_join, by='pos')

  TsMutFreq.ref<- MutFreq_Ts.ref %>% purrr::reduce(full_join, by='pos')
Tv1.MutFreq.ref<-MutFreq_tv1.ref %>% purrr::reduce(full_join, by='pos')
Tv2.MutFreq.ref<-MutFreq_tv2.ref %>% purrr::reduce(full_join, by='pos')
Tvs.MutFreq.ref<-MutFreq_tvs.ref %>% purrr::reduce(full_join, by='pos')
 AllMutFreq.ref<-MutFreq_all.ref %>% purrr::reduce(full_join, by='pos')

 
write.csv(Tv1.MutFreq, "Output1A/MutFreq/Tv1MutFreq_maj_summary.csv")
write.csv(Tv2.MutFreq, "Output1A/MutFreq/Tv2MutFreq_maj_summary.csv")
write.csv(TsMutFreq,    "Output1A/MutFreq/TsMutFreq_maj_summary.csv")
write.csv(Tvs.MutFreq, "Output1A/MutFreq/TvsMutFreq_maj_summary.csv")
write.csv(AllMutFreq,  "Output1A/MutFreq/AllMutFreq_maj_summary.csv")
 
write.csv(Tv1.MutFreq.ref, "Output1A/MutFreq/Tv1MutFreq_ref_summary.csv")
write.csv(Tv2.MutFreq.ref, "Output1A/MutFreq/Tv2MutFreq_ref_summary.csv")
write.csv(  TsMutFreq.ref, "Output1A/MutFreq/TsMutFreq_ref_summary.csv")
write.csv(Tvs.MutFreq.ref, "Output1A/MutFreq/TvsMutFreq_ref_summary.csv")
write.csv( AllMutFreq.ref, "Output1A/MutFreq/AllMutFreq_ref_summary.csv")


###### Check means

s<-length(Overview_summary)
mean(rowMeans(TsMutFreq[2:(s+1)],na.rm=T),na.rm=T) #[1] 0.008012215
mean(rowMeans(Tvs.MutFreq[2:(s+1)],na.rm=T),na.rm=T) # 0.001501265
mean(rowMeans(AllMutFreq[2:(s+1)],na.rm=T),na.rm=T) #  0.00951348

## Ref ones are not meaningful without filtering ##
mean(rowMeans(TsMutFreq.ref[2:(s+1)],na.rm=T),na.rm=T) # 0.06059215
mean(rowMeans(Tvs.MutFreq.ref[2:(s+1)],na.rm=T),na.rm=T) # 0.01262301
mean(rowMeans(AllMutFreq.ref[2:(s+1)],na.rm=T),na.rm=T) #  0.07321516

#############################