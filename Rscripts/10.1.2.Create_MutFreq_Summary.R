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
        if (i==1|i==5) muttypes<-M[,c("pos","ref", "Type","WTAA","MutAA","makesCpG","bigAAChange")]
        if(i==2|i==3) muttypes<-M[,c("pos","ref", paste0("Type",cnames[i]),"WTAA",paste0("MutAA",cnames[i]), paste0("makesCpG",cnames[i]), paste0("bigAAChange",cnames[i]))]
        if(i==4) muttypes<-M[,c("pos","ref","Type.tv1","Type.tv2","WTAA","MutAA.tv1","MutAA.tv2","makesCpG.tv1","makesCpG.tv1","bigAAChange.tv1","bigAAChange.tv2")]
        
        dat2<-merge(dat,muttypes,by="pos")
        mf.files[[i]]<-dat2
        names(mf.files)[i]<-files[i]
        dfname<-paste0(files[i],".F")
        assign(dfname, dat2)
        write.csv(dat2,paste0("Output1A/MutFreq.filtered/Filtered.",files[i],".Q35.csv"))
}

#####
#read in files
Ts<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",row.names = 1,stringsAsFactors = F)
Tv1.MutFreq<-read.csv("Output1A/MutFreq.filtered/Filtered.Tv1.MutFreq.Q35.csv",row.names = 1,stringsAsFactors = F)
Tv2.MutFreq<-read.csv("Output1A/MutFreq.filtered/Filtered.Tv2.MutFreq.Q35.csv",row.names = 1,stringsAsFactors = F)
Tvs.MutFreq<-read.csv("Output1A/MutFreq.filtered/Filtered.Tvs.MutFreq.Q35.csv",row.names = 1,stringsAsFactors = F)
AllMutFreq<-read.csv("Output1A/MutFreq.filtered/Filtered.AllMutFreq.Q35.csv",row.names = 1,stringsAsFactors = F)
###



##### Chceck the mean mut freq
#Summarize the mean and se for all types of mutations

tb<-data.frame(type=c("Ts","Ts.syn","Ts.ns","Ts.stop", "Tv1","Tv1.syn","Tv1.ns","Tv1.stop","Tv2","Tv2.syn","Tv2.ns","Tv2.stop","Tvs","All" ))
files<-c("Ts","Tv1.MutFreq","Tv2.MutFreq","Tvs.MutFreq","AllMutFreq" )

for (i in 1:length(files)){
        dt<-mf.files[[i]]
        #coding region only
        dt<-dt[dt$pos>341,]
        
        if (i<=3){
                colnames(dt)[199]<-"Type"
                k<-(i-1)*4+1
                tb$mean[k]<-mean(dt$mean,na.rm=T)
                tb$se[k]<-std.error(dt$mean,na.rm=T)
                k<-k+1
                tb$mean[k]<-mean(dt$mean[dt$Type=="syn"],na.rm=T)
                tb$se[k]<-std.error(dt$mean[dt$Type=="syn"],na.rm=T)
                k<-k+1
                tb$mean[k]<-mean(dt$mean[dt$Type=="nonsyn"],na.rm=T)
                tb$se[k]<-std.error(dt$mean[dt$Type=="nonsyn"],na.rm=T)
                k<-k+1
                tb$mean[k]<-mean(dt$mean[dt$Type=="stop"],na.rm=T)
                tb$se[k]<-std.error(dt$mean[dt$Type=="stop"],na.rm=T)
        }
        
        if (i==4) {k=13
                tb$mean[k]<-mean(dt$mean,na.rm=T)
                tb$se[k]<-std.error(dt$mean,na.rm=T)}
        if (i==5) {k=14
                tb$mean[k]<-mean(dt$mean,na.rm=T)
                tb$se[k]<-std.error(dt$mean,na.rm=T)}
}

tb$type<-as.character(tb$type)

addtb<-data.frame(type=c("tvs.syn","tvs.ns","tvs.stop", "all.syn","all.ns","all.stop"),
                    mean=c(mean(tb$mean[6],tb$mean[10]),mean(tb$mean[7],tb$mean[11]),mean(tb$mean[8],tb$mean[12]),
                           mean(tb$mean[2],tb$mean[6],tb$mean[10]),mean(tb$mean[3],tb$mean[7],tb$mean[11]),mean(tb$mean[4],tb$mean[8],tb$mean[12])),
                    se  =c(mean(tb$se[6],tb$se[10]),mean(tb$se[7],tb$se[11]),mean(tb$se[8],tb$se[12]),
                           mean(tb$se[2],tb$se[6],tb$se[10]),mean(tb$se[3],tb$se[7],tb$se[11]),mean(tb$se[4],tb$se[8],tb$se[12])))

tb2<-rbind(tb, addtb)

write.csv(tb2, "Output1A/MutFreq.filtered/MF.Mean.SE.summary.csv")



#####
#compare the Q30 and Q35 mut freq
#Q30 mut freq
#q30<-read.csv("Output1A/mut.freq.filtered/Summary_TsMutFreq.csv",stringsAsFactors = F)
#r1<-wilcox.test(Ts$mean[Ts$Type=="syn"], q30$mean[q30$Type=="syn"], alternative = "less", paired = FALSE) 
#W = 2011800, p-value < 2.2e-16
#r1[3]  # p.value =  5.283758e-24




