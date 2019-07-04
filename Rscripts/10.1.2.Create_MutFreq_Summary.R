library(plotrix)
library(reshape)
library(tidyverse)
library(zoo)
library(purrr)


#Create the filtered mut frequency table of all samples 

HCVFiles_overview3<-list.files("Output1A/Overview3/",pattern="overview3.csv")
FilteredOverview2<-list()
for (i in 1:length(HCVFiles_overview3)){ 
        overviews2<-read.csv(paste0("Output1A/Overview3/",HCVFiles_overview3[i]),stringsAsFactors=FALSE, row.names=1)
        FilteredOverview2[[i]]<-overviews2
        names(FilteredOverview2)[i]<-substr(paste(HCVFiles_overview3[i]),start=1,stop=7)
}

##################################
MutFreq_Ts<-list()
MutFreq_tv1<-list()
MutFreq_tv2<-list()
MutFreq_tvs<-list()
MutFreq_all<-list()

for (i in 1:length(FilteredOverview2)){
        dat<-FilteredOverview2[[i]]
        filename<-names(FilteredOverview2)[i]
        
        MutFreq_Ts[[i]]<-dat[,c("pos","freq.Ts.ref")] 
        MutFreq_tv1[[i]]<-dat[,c("pos","freq.transv1.ref")] 
        MutFreq_tv2[[i]]<-dat[,c("pos","freq.transv2.ref")] 
        MutFreq_tvs[[i]]<-dat[,c("pos","freq.transv.ref")] 
        MutFreq_all[[i]]<-dat[,c("pos","freq.mutations.ref")] 
        
        names(MutFreq_Ts)[i]<-filename
        names(MutFreq_tv1)[i]<-filename
        names(MutFreq_tv2)[i]<-filename
        names(MutFreq_tvs)[i]<-filename
        names(MutFreq_all)[i]<-filename
        
        
}
#assign column names for the list
for (i in 1:length(MutFreq_Ts)) {
        colnames(MutFreq_Ts[[i]])<-c("pos",paste0(names(MutFreq_Ts[i])))
        colnames(MutFreq_tv1[[i]])<-c("pos",paste0(names(MutFreq_tv1[i])))
        colnames(MutFreq_tv2[[i]])<-c("pos",paste0(names(MutFreq_tv2[i])))
        colnames(MutFreq_tvs[[i]])<-c("pos",paste0(names(MutFreq_tvs[i])))
        colnames(MutFreq_all[[i]])<-c("pos",paste0(names(MutFreq_all[i])))
}

Ts<-MutFreq_Ts%>% purrr::reduce(full_join, by='pos')
Tv1.MutFreq<-MutFreq_tv1 %>% purrr::reduce(full_join, by='pos')
Tv2.MutFreq<-MutFreq_tv2 %>% purrr::reduce(full_join, by='pos')
Tvs.MutFreq<-MutFreq_tvs %>% purrr::reduce(full_join, by='pos')
AllMutFreq<-MutFreq_all %>% purrr::reduce(full_join, by='pos')

### all mut. freq with metadata ###

files<-c("Ts","Tv1.MutFreq","Tv2.MutFreq","Tvs.MutFreq","AllMutFreq" )
cnames<-c("",".tv1",".tv2",".tvs","")
mf.files<-list()
s<-length(FilteredOverview2)
M<-FilteredOverview2[[3]]
colnames(M)[33:35]<-c("MutAA","MutAA.tv1","MutAA.tv2")
#muttypes1<-M[,c("pos","ref","Type","Type.tv1","Type.tv2","WTAA","MUTAA","TVS1_AA","TVS2_AA","makesCpG","makesCpG.tv1","makesCpG.tv2","bigAAChange","bigAAChange.tv1","bigAAChange.tv2")]
#muttypes<-M[,c(1,4,5,20,23,27,32,36,39)]

for (i in 1:5) {
        dat<-get(files[i])
        dat$mean<-rowMeans(dat[2:(s+1)],na.rm=T, dims=)
        if(i==1|i==2|i==3) muttypes<-M[,c("pos","ref", paste0("Type",cnames[i]),"WTAA",paste0("MutAA",cnames[i]), paste0("makesCpG",cnames[i]), paste0("bigAAChange",cnames[i]))]
        if(i==4|i==5) muttypes<-M[,c("pos","ref","Type","WTAA","MutAA","makesCpG","bigAAChange")]
        
        dat2<-merge(dat,muttypes,by="pos")
        mf.files[[i]]<-dat2
        names(mf.files)[i]<-files[i]
        dfname<-paste0(files[i],".F")
        assign(dfname, dat2)
        write.csv(dat2,paste0("Output1A/MutFreq.filtered/Filtered.",files[i],".Q35.csv"))
}


#AllMutFreq.F<-read.csv("Output1A/MutFreq.filtered/Filtered.AllMutFreq.Q35.csv",row.names = 1,stringsAsFactors = F)

##### Chceck the mean mut freq

mean(rowMeans(Ts.F[,2:(s+1)],na.rm=T),na.rm=T)   # 0.004804249 (Q35),  0.005549325 (Q30 reads<1000)
mean(rowMeans(Tvs.MutFreq.F[,2:(s+1)],na.rm=T),na.rm=T) # 0.0009677105 (Q35),   0.001138303 (Q30 reads<1000)

mean(rowMeans(AllMutFreq.F[,2:(s+1)],na.rm=T),na.rm=T) #0.005771959

std.error(AllMutFreq.F$mean) # 4.563061e-05
std.error(Ts.F$mean) # 4.241551e-05
std.error(Tvs.MutFreq.F$mean) #9.613703e-06
#coding regions only
TMutFreq2.F<-Ts.F[Ts.F$pos>=342,]
mean(TMutFreq2.F$mean, na.rm=T) #  0.004831417
#Syn vs NonSyn
mean(TMutFreq2.F$mean[TMutFreq2.F$Type=="syn"], na.rm=T) #  0.008093379
mean(TMutFreq2.F$mean[TMutFreq2.F$Type=="nonsyn"], na.rm=T) # 0.003351718




#Transition mutations
mean(Ts$mean[Ts$Type=="syn"], na.rm=T) #  0.0075788315 (Q35), #0.00911283 (Q30)
mean(Ts$mean[Ts$Type=="nonsyn"],na.rm=T) # 0.00329428(Q350,  #0.003973061 (Q30)
mean(Ts$mean[Ts$Type=="stop"],na.rm=T)  #0.001928567

mean(Ts$mean[Ts$Type=="syn"&Ts$makesCpG==0],na.rm=T)  #0.007185474
mean(Ts$mean[Ts$Type=="syn"&Ts$makesCpG==1],na.rm=T)  #0.01091608

mean(Ts$mean[Ts$Type=="nonsyn"&Ts$makesCpG==0],na.rm=T) #0.004844019
mean(Ts$mean[Ts$Type=="nonsyn"&Ts$makesCpG==1],na.rm=T) # 0.004355802


#compare the Q30 and Q35 mut freq
#Q30 mut freq
q30<-read.csv("Output/1Amut.freq.filtered/Summary_TsMutFreq.csv",stringsAsFactors = F)
r1<-wilcox.test(Ts$mean[Ts$Type=="syn"], q30$mean[q30$Type=="syn"], alternative = "less", paired = FALSE) 
#W = 2011800, p-value < 2.2e-16
r1[3]  # p.value =  5.283758e-24




