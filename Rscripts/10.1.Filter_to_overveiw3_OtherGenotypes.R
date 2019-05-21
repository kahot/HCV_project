library(reshape)
library(tidyverse)
library(zoo)
library(purrr)

source("Rscripts/baseRscript.R")

#Specify the Genotype
geno<-"1B"
geno<-"3A"

#Load the overview files (summarized using the ref (H77))

HCVFiles_overview2<-list.files(paste0("Output",geno,"/Overview2/"),pattern="overview2.csv")
Overview_sum_ref<-list()
for (i in 1:length(HCVFiles_overview2)){ 
        overviews<-read.csv(paste0("Output",geno,"/Overview2/",HCVFiles_overview2[i]),stringsAsFactors=FALSE)
        overviews<-overviews[,-1]
        Overview_sum_ref[[i]]<-overviews
        names(Overview_sum_ref)[i]<-substr(paste(HCVFiles_overview2[i]),start=1,stop=7)
}

#remove D75213R (read depth not enough)



s<-length(HCVFiles_overview2)
################################################

# 1. Filter mutation frequency >0.2 (this will remove the sites MajNT != ref) and coverage <1000
# put NA for the sites with total-reads <1000 (3/5/2019)
FilteredOverview1<-list()
for (i in 1:length(Overview_sum_ref)){
        dat<-Overview_sum_ref[[i]]
        filename<-names(Overview_sum_ref)[i]
        low_reads<-which(dat$TotalReads<1000) 
        dat[low_reads,c(8:42)]<-NA

        mf_high<-which(dat$freq.Ts.ref>=0.2|dat$freq.transv.ref>=0.2) # row numbers of high mut freq
        dat[mf_high,c(8:42)]<-NA
        FilteredOverview1[[i]]<-dat
        names(FilteredOverview1)[i]<-filename
        #write.csv(dat,paste0("Output",geno,"/Overview3/",filename,"_overview3.csv"))
        
        
}

###################
T_Freq_all<-list()
for (i in 1:length(FilteredOverview1)){
        dat<-FilteredOverview1[[i]]
        filename<-names(FilteredOverview1)[i]
        T_Freq_all[[i]]<-dat[,c("pos","freq.Ts.ref")] 
        names(T_Freq_all)[i]<-filename
}
#assign column names for the list
for (i in 1:length(T_Freq_all)) {
        colnames(T_Freq_all[[i]])<-c("pos",paste0(names(T_Freq_all[i])))
}

TMutFreq<-T_Freq_all %>% purrr::reduce(full_join, by='pos') #8537 sites

write.csv(TMutFreq, paste0("Output",geno,"/MutFreq/TsMutFreqSummary_",geno,".all.csv"))


### Calculate mean transition mut. freq
mean(rowMeans(TMutFreq[2:(s+1)],na.rm=T),na.rm=T)  #1B:0.005223944  1A:0.005135379 3A:0.00443436
TMutFreq$mean<-rowMeans(TMutFreq[2:(s+1)],na.rm=T)
mean(TMutFreq$mean[TMutFreq$pos>=342],na.rm=T) # 1B: 0.005252472    1A: 0.005157235  3A: 0.004449846

dat<-FilteredOverview1[[1]]
dat<-dat[,c(1,4,5,20,23,27,31,32,33,36,39)]
TMutFreq<-merge(TMutFreq,dat, by="pos")


#Syn vs NonSyn
mean(TMutFreq$mean[TMutFreq$Type=="syn" & TMutFreq$pos>=342], na.rm=T) #1B:0.008920395  1A:0.007561845    3A: 0.006943268
mean(TMutFreq$mean[TMutFreq$Type=="nonsyn"& TMutFreq$pos>=342], na.rm=T) #1B: 0.003489323  1A:0.003308032   3A: 0.003266465


# Remove the sites with >30% NA in mutation frequency in FilteredOverview1 files, 
#count the # of non NA samples
TMutFreq$sum.nonNA<-apply(TMutFreq[2:(s+1)],1,function(x) sum(!is.na(x)))
table(TMutFreq$sum.nonNA) #241 sites are all NA
#TMutFreq<-TMutFreq[TMutFreq$sum!=0,] 

TMutFreq$keep0.5<-((TMutFreq$sum)/s)>=0.5
TMutFreq$keep0.3<-((TMutFreq$sum)/s)>=(1/3)

#how many sites?
sum(TMutFreq$keep0.5==T) #7392 (3A)     7733(1B)   7481
sum(TMutFreq$keep0.3==T) #7879 (3A)     7935 (1B)  8025

#what % are REMOVED?
1-sum(TMutFreq$keep0.5==T)/nrow(TMutFreq) #3A:0.1343249    1B: 0.09417828        0.1236968
1-sum(TMutFreq$keep0.3==T)/nrow(TMutFreq) #3A:0.07729242   1B: 0.07051657        0.05880286

#create a vector of positions to keep  : 1/3 for now (3/12/1) 
Keep<-data.frame(TMutFreq$pos[TMutFreq$keep0.3==T])
colnames(Keep)<-"pos"

# Retain only 'Keep'sites
FilteredOverview2<-list()
for (i in 1:length(FilteredOverview1)){
        dat<-FilteredOverview1[[i]]
        filename<-names(FilteredOverview1)[i]
        dat<-merge(Keep, dat, by="pos") #trim down to pos='Keep' 
        FilteredOverview2[[i]]<-dat
        names(FilteredOverview2)[i]<-filename
        write.csv(dat,paste0("Output",geno,"/Overview3/",filename,"_overview3.csv"))
        
}

##################################
MutFreq_Ts<-list()
#MutFreq_tv1<-list()
#MutFreq_tv2<-list()
#MutFreq_tvs<-list()
#MutFreq_all<-list()

for (i in 1:length(FilteredOverview2)){
        dat<-FilteredOverview2[[i]]
        filename<-names(FilteredOverview2)[i]
        
        MutFreq_Ts[[i]]<-dat[,c("pos","freq.Ts.ref")] 
        #MutFreq_tv1[[i]]<-dat[,c("pos","freq.transv1.ref")] 
        #MutFreq_tv2[[i]]<-dat[,c("pos","freq.transv2.ref")] 
        #MutFreq_tvs[[i]]<-dat[,c("pos","freq.transv.ref")] 
        #MutFreq_all[[i]]<-dat[,c("pos","freq.mutations.ref")] 
        
        names(MutFreq_Ts)[i]<-filename
        #names(MutFreq_tv1)[i]<-filename
        #names(MutFreq_tv2)[i]<-filename
        #names(MutFreq_tvs)[i]<-filename
        #names(MutFreq_all)[i]<-filename
        
        
}
#assign column names for the list
for (i in 1:length(MutFreq_Ts)) {
        colnames(MutFreq_Ts[[i]])<-c("pos",paste0(names(MutFreq_Ts[i])))
        #colnames(MutFreq_tv1[[i]])<-c("pos",paste0(names(MutFreq_tv1[i])))
        #colnames(MutFreq_tv2[[i]])<-c("pos",paste0(names(MutFreq_tv2[i])))
        #colnames(MutFreq_tvs[[i]])<-c("pos",paste0(names(MutFreq_tvs[i])))
        #colnames(MutFreq_all[[i]])<-c("pos",paste0(names(MutFreq_all[i])))
}

Ts<-MutFreq_Ts%>% purrr::reduce(full_join, by='pos')
#Tv1.MutFreq<-MutFreq_tv1 %>% purrr::reduce(full_join, by='pos')
#Tv2.MutFreq<-MutFreq_tv2 %>% purrr::reduce(full_join, by='pos')
#Tvs.MutFreq<-MutFreq_tvs %>% purrr::reduce(full_join, by='pos')
#AllMutFreq<-MutFreq_all %>% purrr::reduce(full_join, by='pos')

### Transition only ####

M<-FilteredOverview2[[3]]
muttypes<-M[,c(1,4,5,20,23,27,31,32,33,36,39)]
Ts$mean<-rowMeans(Ts[2:s+1],na.rm=T)
mean(Ts$mean)  #0.004790245

Ts2<-merge(Ts,muttypes,by="pos")
write.csv(Ts2,paste0("Output", geno,"/MutFreq/Filetered.Summary.Ts.",geno,".csv"))







### Other files  ## work on later ## 3.25.19
files<-c("Ts","Tv1.MutFreq","Tv2.MutFreq","Tvs.MutFreq","AllMutFreq" )
mf.files<-list()
s<-length(FilteredOverview2)
M<-FilteredOverview2[[3]]
muttypes<-M[,c("pos","ref","Type","Type.tv1","Type.tv2","WTAA","MUTAA","TVS1_AA","TVS2_AA","makesCpG","makesCpG.tv1","makesCpG.tv2","bigAAChange","bigAAChange.tv1","bigAAChange.tv2")]
#muttypes<-M[,c(1,4,5,20,23,27,32,36,39)]

for (i in 2:5){
        dat<-get(files[i])
        dat$mean<-rowMeans(dat[2:s+1],na.rm=T)
        dat2<-merge(dat,muttypes,by="pos")
        mf.files[[i]]<-dat2
        names(mf.files)[i]<-files[i]
        dfname<-paste0(files[i],".F")
        assign(dfname, dat2)
        write.csv(dat2,paste0("Output/MutFreqQ35/Filetered.Summary.",files[i],".Q35.csv"))
}


##### Chceck the mean mut freq

mean(rowMeans(TMutFreq2.F[,2:(s+1)],na.rm=T),na.rm=T)   # 0.004804249 (Q35),  0.005549325 (Q30 reads<1000)
mean(rowMeans(Tvs.MutFreq.F[,2:(s+1)],na.rm=T),na.rm=T) # 0.0009677105 (Q35),   0.001138303 (Q30 reads<1000)


TMutFreq2.F<-TMutFreq2.F[TMutFreq2.F$pos>=342,]
mean(TMutFreq2.F$mean, na.rm=T) #  0.004831417
#Syn vs NonSyn
mean(TMutFreq2.F$mean[TMutFreq2.F$Type=="syn"], na.rm=T) #  0.007578831
mean(TMutFreq2.F$mean[TMutFreq2.F$Type=="nonsyn"], na.rm=T) #  0.00329428

write.csv(TMutFreq2.F, "Output/TMutFreq2.F.csv")







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
q30<-read.csv("Output/mut.freq.filtered/Summary_TsMutFreq.csv",stringsAsFactors = F)
r1<-wilcox.test(Ts$mean[Ts$Type=="syn"], q30$mean[q30$Type=="syn"], alternative = "less", paired = FALSE) 
#W = 2011800, p-value < 2.2e-16
r1[3]  # p.value =  5.283758e-24




