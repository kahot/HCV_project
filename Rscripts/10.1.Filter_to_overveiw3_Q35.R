library(reshape)
library(tidyverse)
library(zoo)
library(purrr)
source("Rscripts/baseRscript.R")


#Load the overview files (summarized using the ref (H77))

HCVFiles_overview2<-list.files("Output1A/Overview2/",pattern="overview2.csv")
Overview_sum_ref<-list()
for (i in 1:length(HCVFiles_overview2)){ 
        overviews<-read.csv(paste0("Output1A/Overview2/",HCVFiles_overview2[i]),stringsAsFactors=FALSE)
        overviews<-overviews[,-1]
        Overview_sum_ref[[i]]<-overviews
        names(Overview_sum_ref)[i]<-substr(paste(HCVFiles_overview2[i]),start=1,stop=7)
}

################################################

# 1. Filter mutation frequency >0.2 (this will remove the sites MajNT != ref) and coverage <1000
# put NA for the sites with total-reads <1000 (3/5/2019)
FilteredOverview1<-list()
for (i in 1:length(Overview_sum_ref)){
        dat<-Overview_sum_ref[[i]]
        filename<-names(Overview_sum_ref)[i]
        low_reads<-which(dat$TotalReads<1000) 
        dat[low_reads,c(8:19)]<-NA

        mf_high<-which(dat$freq.Ts.ref>=0.2|dat$freq.transv.ref>=0.2) # row numbers of high mut freq
        dat[mf_high,c(8:19)]<-NA
        FilteredOverview1[[i]]<-dat
        names(FilteredOverview1)[i]<-filename
}

####
## Create the summary of filtered mutation frequnecy
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

#write.csv(TMutFreq, "Output1A/MutFreq.filtered/Ts_MutFreq_1A_all.csv")


##### Check the results   ####
#       ### Calculate mean transition mut. freq
#       mean(rowMeans(TMutFreq[2:196],na.rm=T),na.rm=T) 
#       TMutFreq$mean<-rowMeans(TMutFreq[2:196],na.rm=T)
#       mean(TMutFreq$mean[TMutFreq$pos>=342],na.rm=T) # 0.005157235
#       
#       dat<-FilteredOverview1[[1]]
#       dat<-dat[,c(1,4,5,20,23,27,31,32,33,36,39)]
#       TMutFreq<-merge(TMutFreq,dat, by="pos")        
#       #Syn vs NonSyn
#       mean(TMutFreq$mean[TMutFreq$Type=="syn" & TMutFreq$pos>=342], na.rm=T) # 0.007561845
#       mean(TMutFreq$mean[TMutFreq$Type=="nonsyn"& TMutFreq$pos>=342], na.rm=T) # 0.003308032

###  Remove the sites with >50% NA in mutation frequency in FilteredOverview1 files, 
s<-length(FilteredOverview1)

#count the # of non NA samples
TMutFreq$sum.nonNA<-apply(TMutFreq[2:(s+1)],1,function(x) sum(!is.na(x)))
table(TMutFreq$sum.nonNA) #169 sites are all NA

TMutFreq$keep0.5<-((TMutFreq$sum)/s)>=0.5
TMutFreq$keep0.3<-((TMutFreq$sum)/s)>=(1/3)

        #how many sites?
        sum(TMutFreq$keep0.5==T) #7481
        sum(TMutFreq$keep0.3==T) #8035
        
        #what % of sites are REMOVED?
        1-sum(TMutFreq$keep0.5==T)/nrow(TMutFreq) #0.1236968
        1-sum(TMutFreq$keep0.3==T)/nrow(TMutFreq) #0.05880286
        
#create a vector of positions to keep  : >1/3 for now
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
        write.csv(dat,paste0("Output1A/Overview3/",filename,"_overview3.csv"))
        
}


