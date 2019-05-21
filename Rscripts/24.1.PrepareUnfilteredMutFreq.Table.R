library(stringr)
library(tidyverse)
library(zoo)

#############
# find the sites that are most different in mutation frequencies between the genomes

# create non-filtered summary file

#MergedPositions<-read.csv("Data/MergedPositionInfo.csv", stringsAsFactors = F)
geno<-c("1A","1B","3A")
#
#HCV1A<-list.files(paste0("Output/Overview2/"),pattern="overview2.csv")
#HCV3A<-list.files(paste0("Output3A/Overview2/"),pattern="overview2.csv")
#HCV1B<-list.files(paste0("Output1B/Overview2/"),pattern="overview2.csv")
#
#for (g in 1:3){
#        flist<-get(paste0("HCV",geno[g]))
#        if (g==1) dir<-"Output/Overview2/"
#        if (g==2) dir<-"Output1B/Overview2/"
#        if (g==3) dir<-"Output3A/Overview2/"
#        
#        Positions<-MergedPositions[,c("merged.pos",paste0("org.pos.",geno[g]))]
#        pos<-data.frame(merged.pos=Positions$merged.pos)
#        for (i in 1:length(flist)){
#                dat<-read.csv(paste0(dir,flist[i]), stringsAsFactors = F)
#                dat<-dat[,-1]
#                colnames(dat)[1]<-paste0("org.pos.",geno[g])
#                
#                dat2<-merge(Positions,dat, by=paste0("org.pos.",geno[g]),all.x=T )
#                dat3<-merge(pos, dat2, by="merged.pos", all.x=T)
#                
#                fname<-substr(paste(flist[i]),start=1,stop=7)
#                write.csv(dat3,paste0("Output_all/Overview",geno[g],"/Overview2.2/", fname, "_overview2.2.csv"))
#        }
#}
#


#### create a transition mutation summary (ref) of Overview2.2 (unfiltered) ###

#create merged metadata since some are missing from "Output_all/Ts_summary_metadata.1A.csv"

merged.meta1<-read.csv(paste0("Output_all/Overview1A/Overview2.2/D75002-_overview2.2.csv"), stringsAsFactors = F)
merged.meta1<-merged.meta1[,c(2,3,6,7,22,38,41)]
colnames(merged.meta1)[3:7]<-paste0(colnames(merged.meta1)[3:7],".1A")
merged.meta2<-read.csv(paste0("Output_all/Overview1B/Overview2.2/D75046-_overview2.2.csv"), stringsAsFactors = F)
merged.meta2<-merged.meta2[,c(2,3,6,7,22,38,41)]
colnames(merged.meta2)[3:7]<-paste0(colnames(merged.meta2)[3:7],".1B")
merged.meta3<-read.csv(paste0("Output_all/Overview3A/Overview2.2/D75007-_overview2.2.csv"), stringsAsFactors = F)
merged.meta3<-merged.meta3[,c(2,3,6,7,22,38,41)]
colnames(merged.meta3)[3:7]<-paste0(colnames(merged.meta3)[3:7],".3A")

merged.meta<-cbind(merged.meta1,merged.meta2[2:7],merged.meta3[2:7])

genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}

genetable<-data.frame("merged.pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector
merged.meta<-merge(genetable,merged.meta, by="merged.pos")
#add codon positions
codon<-rep(1:3, times=nrow(merged.meta)/3)
codon<-c(2,3,codon)

merged.meta<-cbind(codon[1:nrow(merged.meta)],merged.meta)
merged.meta<-merged.meta[c(2,1,3:ncol(merged.meta))]
colnames(merged.meta)[2]<-"codon"

for (f in 1:3){
        flist<-list.files(paste0("Output_all/Overview",geno[f],"/Overview2.2/"),pattern="overview2.2.csv")
        
        MutFreq_Ts<-list()
        for (i in 1:length(flist)){
                dat<-read.csv(paste0("Output_all/Overview",geno[f],"/Overview2.2/",flist[i]), stringsAsFactors = F)
                dat<-dat[,-1]
                #filter out the sites with read depth <1000
                low_reads<-which(dat$TotalReads<1000) 
                dat[low_reads,c(7:ncol(dat))]<-NA
                
                #get the Ts. mutation frequency
                MutFreq_Ts[[i]]<-dat[,c("merged.pos","freq.Ts.ref")] 
                
                filename<-substr(paste0(flist[i]),start=1,stop=7)
                names(MutFreq_Ts)[i]<-filename
        }
        for (i in 1:length(MutFreq_Ts)) {
                colnames(MutFreq_Ts[[i]])<-c("merged.pos",paste0(names(MutFreq_Ts[i])))
        }
        
        Ts<-MutFreq_Ts%>% purrr::reduce(full_join, by='merged.pos')
        Ts$mean<-rowMeans(Ts[2:ncol(Ts)],na.rm=T)
        #write.csv(Ts,paste0("Output_all/Unfiltered/Ts.ref_",geno[f],".csv"))
        s<-length(flist)
        Ts2<-Ts[,c(1,s+2)]
        colnames(Ts2)[2]<-paste0("mean.",geno[f]) 
        Ts2<-merge(merged.meta,Ts2,by="merged.pos")
        write.csv(Ts2,paste0("Output_all/Unfiltered/Ts.ref_Mean_metadata.",geno[f],".csv"))
}

#Combine the mean mut freq of 3 genotypes into 1 table.

Summary<-read.csv(paste0("Output_all/Unfiltered/Ts.ref_Mean_metadata.1A.csv"),stringsAsFactors = F)
Summary<-Summary[,-1]

for (f in 2:3){
        d<-read.csv(paste0("Output_all/Unfiltered/Ts.ref_Mean_metadata.",geno[f],".csv"),stringsAsFactors = F)
        d<-d[,c(2,ncol(d))]
        Summary<-merge(Summary,d,by="merged.pos")
}


write.csv(Summary,"Output_all/Unfiltered/Ts.ref_Mean_3genotypes.csv")



range(Summary$mean.1A, na.rm=T)
range(Summary$mean.3A, na.rm=T)



#################################
# separate the sites that are Maj!=Ref

for (f in 1:3){
        flist<-list.files(paste0("Output_all/Overview",geno[f],"/Overview2.2/"),pattern="overview2.2.csv")
        
        MutFreq_Ts.same<-list()
        MutFreq_Ts.diff<-list()
        pos<-data.frame(merged.pos=c(261:8699))
        for (i in 1:length(flist)){
                dat<-read.csv(paste0("Output_all/Overview",geno[f],"/Overview2.2/",flist[i]), stringsAsFactors = F)
                dat<-dat[,-1]
                #filter out the sites with read depth <1000
                low_reads<-which(dat$TotalReads<1000) 
                dat[low_reads,c(9:ncol(dat))]<-NA
                
                filename<-substr(paste0(flist[i]),start=1,stop=7)
                
                dat.same<-dat[dat$MajNt==paste0(dat$ref),]
                dat.same<-dat.same[!is.na(dat.same$merged.pos),]
                dat.s<-merge(pos,dat.same, by="merged.pos", all.x = T)
                
                MutFreq_Ts.same[[i]]<-dat.s[,c("merged.pos","freq.Ts.ref")] 
                names(MutFreq_Ts.same)[i]<-filename
                
                
                dat.diff<-dat[dat$MajNt!=dat$ref,]
                dat.diff<-dat.diff[!is.na(dat.diff$merged.pos),]
                dat.d<-merge(pos,dat.diff, by="merged.pos", all.x = T)
                
                MutFreq_Ts.diff[[i]]<-dat[,c("merged.pos","freq.Ts")] 
                names(MutFreq_Ts.diff)[i]<-filename
        }
        
        for (i in 1:length(MutFreq_Ts.same)) {
                colnames(MutFreq_Ts.same[[i]])<-c("merged.pos",paste0(names(MutFreq_Ts.same[i])))
                colnames(MutFreq_Ts.diff[[i]])<-c("merged.pos",paste0(names(MutFreq_Ts.diff[i])))
        }
        
        Ts.same<-MutFreq_Ts.same %>% purrr::reduce(full_join, by='merged.pos')
        Ts.same$mean<-rowMeans(Ts.same[2:ncol(Ts.same)],na.rm=T)
        write.csv(Ts.same,paste0("Output_all/Unfiltered/Ts.same_",geno[f],".csv"))
        
        Ts.diff<-MutFreq_Ts.diff %>% purrr::reduce(full_join, by='merged.pos')
        Ts.diff$mean<-rowMeans(Ts.diff[2:ncol(Ts.diff)],na.rm=T)
        write.csv(Ts.diff,paste0("Output_all/Unfiltered/Ts.diff_",geno[f],".csv"))

        s<-length(flist)
        Ts2.same<-Ts.same
        Ts2.same$sum<-apply(Ts2.same[2:(s+1)],1,function(x) sum(!is.na(x)))
        Ts2.same$keep0.3<-(Ts2.same$sum/s)>=1/3
        
        #what % are REMOVED?
        x<-1-sum(Ts2.same$keep0.3==T)/nrow(Ts2.same) 
        cat("% of sites removed by removing sites at least 1/3 samples was ", x) 
        Ts2.same$mean[Ts2.same$keep0.3==F]<-NA
        
        Ts2.same<-Ts2.same[,c("merged.pos","mean")]
        
        colnames(Ts2.same)[2]<-paste0("mean.",geno[f]) 
        Ts2.same<-merge(merged.meta,Ts2.same,by="merged.pos")
        write.csv(Ts2.same,paste0("Output_all/Unfiltered/Ts.same_Mean_metadata.",geno[f],".csv"))


        
        Ts2.diff<-Ts.diff
        Ts2.diff$sum<-apply(Ts2.diff[2:(s+1)],1,function(x) sum(!is.na(x)))
        Ts2.diff$keep0.3<-(Ts2.diff$sum/s)>=1/3

        Ts2.diff<-Ts.diff[,c("merged.pos","mean")]
        colnames(Ts2.diff)[2]<-paste0("mean.",geno[f]) 
        Ts2.diff<-merge(merged.meta,Ts2.diff,by="merged.pos")
        write.csv(Ts2.diff,paste0("Output_all/Unfiltered/Ts.diff_Mean_metadata.",geno[f],".csv"))
        
}

#proportion of sites removed by removing sites at least 1/3 samples was  0.04621401 : 1A
#proportion of sites removed by removing sites at least 1/3 samples was  0.05853774  : 1B
#proportion of sites removed by removing sites at least 1/3 samples was  0.0652921   :3A



#Combine the mean mut freq of 3 genotypes into 1 table.

Summary.same<-read.csv(paste0("Output_all/Unfiltered/Ts.same_Mean_metadata.1A.csv"),stringsAsFactors = F)
Summary.same<-Summary.same[,-1]

for (f in 2:3){
        d<-read.csv(paste0("Output_all/Unfiltered/Ts.same_Mean_metadata.",geno[f],".csv"),stringsAsFactors = F)
        d<-d[,c(2,ncol(d))]
        Summary.same<-merge(Summary.same,d,by="merged.pos")
}
write.csv(Summary.same,"Output_all/Unfiltered/Ts.Same_Mean_3genotypes.csv")




Summary.diff<-read.csv(paste0("Output_all/Unfiltered/Ts.diff_Mean_metadata.1A.csv"),stringsAsFactors = F)
Summary.diff<-Summary.diff[,-1]

for (f in 2:3){
        d<-read.csv(paste0("Output_all/Unfiltered/Ts.diff_Mean_metadata.",geno[f],".csv"),stringsAsFactors = F)
        d<-d[,c(2,ncol(d))]
        Summary.diff<-merge(Summary.diff,d,by="merged.pos")
}


write.csv(Summary.diff,"Output_all/Unfiltered/Ts.Diff_Mean_3genotypes.csv")

