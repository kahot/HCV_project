library(stringr)
library(ape)
library(seqinr)
library(e1071)
library(msa)
library(car)
library(readtext)
library(purrr)

geno<-c("1A","1B","3A")
# 1) Create merged overview files for 3 genotyes (.overview4.csv)

refs<-read.dna("Data/HCV_3GenoRef.fasta", format = "fasta",as.character=TRUE)

Seq3A<-as.character(refs[[1]])
Seq1B<-as.character(refs[[3]])
Seq1A<-as.character(refs[[2]])

for (g in 1:length(geno)){
        nuc<-get(paste0("Seq",geno[g]))
        fname<-paste0("seq",geno[g])
        data<-data.frame("nuc"=nuc)
        data$merged.pos<-1:nrow(data)
        data$org.pos<-""
        n=1
        for (i in 1:nrow(data)){
                if (data$nuc[i]!="-") {
                        data$org.pos[i]<-n
                        n<-n+1
                }
                else  data$org.pos[i]<-NA
        }
        data$org.pos<-as.integer(data$org.pos)
        assign(fname,data)
}

colnames(seq1A)[3]<-"org.pos.1A"
colnames(seq1B)[3]<-"org.pos.1B"
colnames(seq3A)[3]<-"org.pos.3A"
MergedPositions<-cbind(seq1A,seq1B$org.pos.1B)
MergedPositions<-cbind(MergedPositions,seq3A$org.pos.3A)
MergedPositions<-MergedPositions[,-1]
colnames(MergedPositions)[3:4]<-c("org.pos.1B","org.pos.3A")
write.csv(MergedPositions,"Data/MergedPositionInfo.csv")

for (f in 1:3){
        flist<-list.files(paste0("Output",geno[f],"/Overview3/"),pattern="overview3.csv")
        Positions<-get(paste0("seq",geno[f]))
        for (i in 1:length(flist)){
                dat<-read.csv(paste0("Output",geno[f],"/Overview3/",flist[i]), stringsAsFactors = F)
                dat<-dat[,-1]
                colnames(dat)[1]<-paste0("org.pos.",geno[f])
                
                pos<-data.frame(merged.pos=Positions$merged.pos)
                dat2<-merge(Positions,dat, by=paste0("org.pos.",geno[f]),all.x=T )
                dat3<-merge(pos, dat2, by="merged.pos", all.x=T)
                
                fname<-substr(paste(flist[i]),start=1,stop=7)
                write.csv(dat3,paste0("Output_all/Overview",geno[f],"/", fname, "_overview4.csv"))
        }
}


### 2) create a transition mutation summary (filtered data (mut req<0.2)) ###
M<-read.csv("Output_all/merged.metadata.csv", stringsAsFactors = F)
M<-M[,-1]
for (f in 1:3){
    flist<-list.files(paste0("Output_all/Overview",geno[f],"/"),pattern="overview4.csv")
    
    MutFreq_Ts<-list()
    for (i in 1:length(flist)){
        dat<-read.csv(paste0("Output_all/Overview",geno[f],"/",flist[i]), stringsAsFactors = F)
        filename<-substr(paste0(flist[i]),start=1,stop=7)
        MutFreq_Ts[[i]]<-dat[,c("merged.pos","freq.Ts.ref")] 
        names(MutFreq_Ts)[i]<-filename
    }
    for (i in 1:length(MutFreq_Ts)) {
            colnames(MutFreq_Ts[[i]])<-c("merged.pos",paste0(names(MutFreq_Ts[i])))
    }
    
    Ts<-MutFreq_Ts%>% purrr::reduce(full_join, by='merged.pos')
    Ts$mean<-rowMeans(Ts[2:ncol(Ts)],na.rm=T)
    write.csv(Ts,paste0("Output_all/Filtered/Ts_summary_",geno[f],".csv"))
    s<-length(flist)
    Ts2<-Ts[,c(1,s+2)]
    Ts2<-merge(M,Ts2,by="merged.pos")
    write.csv(Ts2,paste0("Output_all/Filtered/Ts_summary_metadata.",geno[f],".csv"))
}

