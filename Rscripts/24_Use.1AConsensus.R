# 24.1 Continuation: use Consensus seq as a reference for 1A

merged.meta<-read.csv("Output_all/Unfiltered/merged.metadata.csv",stringsAsFactors = F,row.names = 1)

cons1A<-read.dna("Data/HCV1A_CONSENSUS.fasta", format="fasta")
cons1A<-cons1A[1:8800]
con1A<-data.frame(ref.1A.con=paste0(cons1A), org.pos.1A=1:8800)
con1A$ref.1A.con<-as.character(con1A$ref.1A.con)
con1A<-merge(con1A, merged.meta, by="org.pos.1A", all.y = T)
con1A<-con1A[,c("merged.pos", "ref.1A.con")]


f=1
flist<-list.files(paste0("Output_all/Overview",geno[f],"/Overview2.2/"),pattern="overview2.2.csv")

MutFreq_Ts.same<-list()
MutFreq_Ts.diff<-list()
MutFreq_Tvs.same<-list()
MutFreq_all.same<-list()
pos<-data.frame(merged.pos=c(261:8699))
for (i in 1:length(flist)){
        dat<-read.csv(paste0("Output_all/Overview",geno[f],"/Overview2.2/",flist[i]), stringsAsFactors = F)
        dat<-dat[,-1]
        #filter out the sites with read depth <1000
        low_reads<-which(dat$TotalReads<1000) 
        dat[low_reads,c(9:ncol(dat))]<-NA
        
        filename<-substr(paste0(flist[i]),start=1,stop=7)
        #for 1A only
        dat<-merge(dat, con1A, by="merged.pos")
        
        dat.same<-dat[dat$MajNt==dat$ref.1A.con,]
        dat.same<-dat.same[!is.na(dat.same$merged.pos),]
        dat.s<-merge(pos,dat.same, by="merged.pos", all.x = T)
        
        MutFreq_Ts.same[[i]]<-dat.s[,c("merged.pos","freq.Ts.ref")] 
        names(MutFreq_Ts.same)[i]<-filename
        
        MutFreq_Tvs.same[[i]]<-dat.s[,c("merged.pos","freq.transv.ref")] 
        MutFreq_all.same[[i]]<-dat.s[,c("merged.pos","freq.mutations.ref")] 
        names(MutFreq_Tvs.same)[i]<-filename
        names(MutFreq_all.same)[i]<-filename
        
        dat.diff<-dat[dat$MajNt!=dat$ref.1A.con,]
        dat.diff<-dat.diff[!is.na(dat.diff$merged.pos),]
        dat.d<-merge(pos,dat.diff, by="merged.pos", all.x = T)
        
        MutFreq_Ts.diff[[i]]<-dat.d[,c("merged.pos","freq.Ts")] 
        names(MutFreq_Ts.diff)[i]<-filename
}

for (i in 1:length(MutFreq_Ts.same)) {
        colnames(MutFreq_Ts.same[[i]])<-c("merged.pos",paste0(names(MutFreq_Ts.same[i])))
        colnames(MutFreq_Ts.diff[[i]])<-c("merged.pos",paste0(names(MutFreq_Ts.diff[i])))
        colnames(MutFreq_Tvs.same[[i]])<-c("merged.pos",paste0(names(MutFreq_Tvs.same[i])))
        colnames(MutFreq_all.same[[i]])<-c("merged.pos",paste0(names(MutFreq_all.same[i])))}

Ts.sameC<-MutFreq_Ts.same %>% purrr::reduce(full_join, by='merged.pos')
Ts.sameC$mean<-rowMeans(Ts.sameC[2:ncol(Ts.sameC)],na.rm=T)
write.csv(Ts.sameC,paste0("Output_all/Unfiltered/Ts.sameC_",geno[f],".csv"))

Ts.diffC<-MutFreq_Ts.diff %>% purrr::reduce(full_join, by='merged.pos')
Ts.diffC$mean<-rowMeans(Ts.diffC[2:ncol(Ts.diffC)],na.rm=T)
write.csv(Ts.diffC,paste0("Output_all/Unfiltered/Ts.diffC_",geno[f],".csv"))

Tvs.sameC<-MutFreq_Tvs.same %>% purrr::reduce(full_join, by='merged.pos')
Tvs.sameC$mean<-rowMeans(Tvs.sameC[2:ncol(Tvs.sameC)],na.rm=T)
write.csv(Tvs.sameC,paste0("Output_all/Unfiltered/Tvs.sameC_",geno[f],".csv"))

All.sameC<-MutFreq_all.same %>% purrr::reduce(full_join, by='merged.pos')
All.sameC$mean<-rowMeans(All.sameC[2:ncol(All.sameC)],na.rm=T)
write.csv(Ts.same,paste0("Output_all/Unfiltered/All.sameC_",geno[f],".csv"))


s<-length(flist)
Ts2.same<-Ts.sameC
Ts2.same$sum<-apply(Ts2.same[2:(s+1)],1,function(x) sum(!is.na(x)))
Ts2.same$keep0.3<-(Ts2.same$sum/s)>=1/3
Ts2.same<-Ts2.same[Ts2.same$keep0.3==T,]

Tvs2.same<-Tvs.sameC
Tvs2.same$sum<-apply(Tvs2.same[2:(s+1)],1,function(x) sum(!is.na(x)))
Tvs2.same$keep0.3<-(Tvs2.same$sum/s)>=1/3
Tvs2.same<-Tvs2.same[Tvs2.same$keep0.3==T,]

All2.same<-All.sameC
All2.same$sum<-apply(All2.same[2:(s+1)],1,function(x) sum(!is.na(x)))
All2.same$keep0.3<-(All2.same$sum/s)>=1/3
All2.same<-All2.same[All2.same$keep0.3==T,]

#what % are REMOVED?
x<-1-sum(Ts2.same$keep0.3==T)/nrow(Ts.sameC) 
cat("% removed by removing sites at least 1/3 samples in ", geno[f], " was ", x) 
Ts2.same$mean[Ts2.same$keep0.3==F]<-NA
Tvs2.same$mean[Tvs2.same$keep0.3==F]<-NA
All2.same$mean[All2.same$keep0.3==F]<-NA

Ts2.same<-Ts2.same[,c("merged.pos","mean")]
Tvs2.same<-Tvs2.same[,c("merged.pos","mean")]
All2.same<-All2.same[,c("merged.pos","mean")]

colnames(Ts2.same)[2]<-paste0("mean.",geno[f]) 
colnames(Tvs2.same)[2]<-paste0("mean.",geno[f]) 
colnames(All2.same)[2]<-paste0("mean.",geno[f]) 

Ts2.same<-merge(merged.meta,Ts2.same,by="merged.pos", all.x=T)
write.csv(Ts2.same,paste0("Output_all/Unfiltered/Ts.same.Consensus_Mean_metadata.",geno[f],".csv"))

Tvs2.same<-merge(merged.meta,Tvs2.same,by="merged.pos", all.x=T)
write.csv(Tvs2.same,paste0("Output_all/Unfiltered/Tvs.same.Consensus_Mean_metadata.",geno[f],".csv"))
All2.same<-merge(merged.meta,All2.same,by="merged.pos", all.x=T)
write.csv(All2.same,paste0("Output_all/Unfiltered/All.same.Consensus_Mean_metadata.",geno[f],".csv"))


Ts2.diff<-Ts.diff[,c("merged.pos","mean")]
colnames(Ts2.diff)[2]<-paste0("mean.",geno[f]) 
Ts2.diff<-merge(merged.meta,Ts2.diff,by="merged.pos", all.x=T)
write.csv(Ts2.diff,paste0("Output_all/Unfiltered/Ts.diff.Consensus_Mean_metadata.",geno[f],".csv"))


nrow(Ts2.same[!is.na(Ts2.same$mean.1A),])
