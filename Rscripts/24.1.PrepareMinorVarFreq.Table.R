library(stringr)
library(tidyverse)
library(zoo)

#############
# find the sites that are most different in mutation frequencies between the genomes

# create non-filtered summary file

MergedPositions<-read.csv("Data/MergedPositionInfo.csv", stringsAsFactors = F)
geno<-c("1A","1B","3A")

HCV1A<-list.files(paste0("Output/Overview2/"),pattern="overview2.csv")
HCV3A<-list.files(paste0("Output3A/Overview2/"),pattern="overview2.csv")
HCV1B<-list.files(paste0("Output1B/Overview2/"),pattern="overview2.csv")

for (g in 1:3){
        flist<-get(paste0("HCV",geno[g]))
        if (g==1) dir<-"Output/Overview2/"
        if (g==2) dir<-"Output1B/Overview2/"
        if (g==3) dir<-"Output3A/Overview2/"
        
        Positions<-MergedPositions[,c("merged.pos",paste0("org.pos.",geno[g]))]
        pos<-data.frame(merged.pos=Positions$merged.pos)
        for (i in 1:length(flist)){
                dat<-read.csv(paste0(dir,flist[i]), stringsAsFactors = F)
                dat<-dat[,-1]
                colnames(dat)[1]<-paste0("org.pos.",geno[g])
                
                dat2<-merge(Positions,dat, by=paste0("org.pos.",geno[g]),all.x=T )
                dat3<-merge(pos, dat2, by="merged.pos", all.x=T)
                
                fname<-substr(paste(flist[i]),start=1,stop=7)
                write.csv(dat3,paste0("Output_all/Overview",geno[g],"/Overview2.2/", fname, "_overview2.2.csv"))
        }
}

### create a transition mutation summary of Overview2.2 (unfiltered) ###

merged.meta<-read.csv("Output_all/Ts_summary_metadata.1A.csv", stringsAsFactors = F)
merged.meta<-merged.meta[,c(2:33)]


for (f in 1:3){
        flist<-list.files(paste0("Output_all/Overview",geno[f],"/Overview2.2/"),pattern="overview2.2.csv")
        
        MutFreq_Ts<-list()
        for (i in 1:length(flist)){
                dat<-read.csv(paste0("Output_all/Overview",geno[f],"/Overview2.2/",flist[i]), stringsAsFactors = F)
                dat<-dat[,-1]
                #filter out the sites with read depth <1000
                low_reads<-which(dat$TotalReads<1000) 
                dat[low_reads,c(7:ncol(dat3))]<-NA
                
                #get the minor variant frequency
                MutFreq_Ts[[i]]<-dat[,c("merged.pos","freq.Ts")] 
                
                filename<-substr(paste0(flist[i]),start=1,stop=7)
                names(MutFreq_Ts)[i]<-filename
        }
        for (i in 1:length(MutFreq_Ts)) {
                colnames(MutFreq_Ts[[i]])<-c("merged.pos",paste0(names(MutFreq_Ts[i])))
        }
        
        Ts<-MutFreq_Ts%>% purrr::reduce(full_join, by='merged.pos')
        Ts$mean<-rowMeans(Ts[2:ncol(Ts)],na.rm=T)
        write.csv(Ts,paste0("Output_all/Unfiltered/Ts_",geno[f],".csv"))
        s<-length(flist)
        Ts2<-Ts[,c(1,s+2)]
        colnames(Ts2)[2]<-paste0("mean.",geno[f]) 
        Ts2<-merge(merged.meta,Ts2,by="merged.pos")
        write.csv(Ts2,paste0("Output_all/Unfiltered/Ts_Mean_metadata.",geno[f],".csv"))
}

#Combine the mean mut freq of 3 genotypes into 1 table.

Summary<-read.csv(paste0("Output_all/Unfiltered/Ts_Mean_metadata.1A.csv"),stringsAsFactors = F)
Summary<-Summary[,-1]

for (f in 2:3){
        d<-read.csv(paste0("Output_all/Unfiltered/Ts_Mean_metadata.",geno[f],".csv"),stringsAsFactors = F)
        d<-d[,c(2,ncol(d))]
        Summary<-merge(Summary,d,by="merged.pos")
}


write.csv(Summary,"Output_all/Unfiltered/Ts_Mean_3genotypes.csv")
