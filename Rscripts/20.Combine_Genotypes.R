library(stringr)
library(ape)
library(seqinr)
library(e1071)
library(msa)
library(car)
library(readtext)


refs<-read.dna(paste0("Data/HCV_3refseq2.fasta"), format = "fasta",as.character=TRUE)
#refs<-read.dna(paste0("Data/References_aligned.fasta"), format = "fasta",as.character=TRUE)

# 
Seq3A<-as.character(refs[[1]])
Seq1B<-as.character(refs[[2]])
Seq1A<-as.character(refs[[3]])


geno<-c("1A","1B","3A")
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

library(Hmisc)
HCV1A<-list.files(paste0("Output/Overview3/"),pattern="overview3.csv")
HCV3A<-list.files(paste0("Output3A/Overview3/"),pattern="overview3.csv")
HCV1B<-list.files(paste0("Output1B/Overview3/"),pattern="overview3.csv")

for (f in 1:3){
        flist<-get(paste0("HCV",geno[f]))
        if (f==1) dir<-"Output3A/Overview3/"
        if (f==2) dir<-"Output/Overview3/"
        if (f==3) dir<-"Output1B/Overview3/"
        
        
        for (i in 1:length(flist)){
                dat<-read.csv(paste0(dir,flist[i]), stringsAsFactors = F)
                dat<-dat[,-1]
                colnames(dat)[1]<-"org.pos"
                
                Positions<-get(paste0("seq",geno[f]))
                pos<-data.frame(merged.pos=Positions$merged.pos)
                dat2<-merge(Positions,dat, by="org.pos",all.x=T )
                dat3<-merge(pos, dat2, by="merged.pos", all.x=T)
                
                fname<-substr(paste(flist[i]),start=1,stop=7)
                write.csv(dat3,paste0("Output_all/Overview",geno[f],"/", fname, "_overview4.csv"))
        }

        
        
}

### create a transition mutation summary ###


M1<-read.csv("Output_all/Overview3A/D75003-_overview4.csv", stringsAsFactors = F)
M2<-read.csv("Output_all/Overview1A/D75002-_overview4.csv", stringsAsFactors = F)
M3<-read.csv("Output_all/Overview1B/D75046-_overview4.csv", stringsAsFactors = F)

for (f in 1:3){
        m<-get(paste0("M",f))
        muttypes1<-m[,c(2,3,4,23,26,30,34,35,36,39,42)]
        columnnames1<-paste0(colnames(muttypes1),".",geno[f])
        colnames(muttypes1)[2:11]<-columnnames1[2:11]
        filename<-paste0("muttypes",geno[f])
        assign(filename,muttypes1)
        
}
M<-merge(muttypes3A, muttypes1A, by="merged.pos")
M<-merge(M, muttypes1B, by="merged.pos")

genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}

genetable<-data.frame("merged.pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector
M<-merge(genetable,M, by="merged.pos")




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
    write.csv(Ts,paste0("Output_all/Ts_summary_",geno[f],".csv"))
    s<-length(flist)
    Ts2<-Ts[,c(1,s+2)]
    Ts2<-merge(M,Ts2,by="merged.pos")
    write.csv(Ts2,paste0("Output_all/Ts_summary_metadata.",geno[f],".csv"))
}

