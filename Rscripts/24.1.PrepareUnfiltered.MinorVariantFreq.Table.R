library(stringr)
library(tidyverse)
library(zoo)
library(reshape2)
library(colorspace)
library(RColorBrewer)
library(sfsmisc)

source("Rscripts/baseRscript.R")

#Geonotype colors 
colors2<-qualitative_hcl(6, palette="Dark3")
# 1A=1, 1B=3, 3A=5, c("#E16A86","#50A315","#009ADE")

geno<-c("1A","1B","3A")

#If run previously, skip to ## **START** ##

#############
# find the sites that are most different in mutation frequencies between the genomes

#1. Create non-filtered Overview files

MergedPositions<-read.csv("Data/MergedPositionInfo.csv", stringsAsFactors = F)
geno<-c("1A","1B","3A")

for (i in 1:3){
        lname<-paste0("HCV",geno[i])
        lS<-list.files(paste0("Output",geno[i],"/Overview2/"), pattern="overview2.csv")
        assign(lname,lS)
}

for (g in 1:3){
        flist<-get(paste0("HCV",geno[g]))
        
        Positions<-MergedPositions[,c("merged.pos",paste0("org.pos.",geno[g]))]
        pos<-data.frame(merged.pos=Positions$merged.pos)
        for (i in 1:length(flist)){
                dat<-read.csv(paste0("Output",geno[g],"/Overview2/",flist[i]), stringsAsFactors = F, row.names = 1)
                colnames(dat)[1]<-paste0("org.pos.",geno[g])
                
                dat2<-merge(Positions,dat, by=paste0("org.pos.",geno[g]),all.x=T )
                dat3<-merge(pos, dat2, by="merged.pos", all.x=T)
                
                fname<-substr(paste(flist[i]),start=1,stop=7)
                write.csv(dat3,paste0("Output_all/Overview",geno[g],"/Overview2.2/", fname, "_overview2.2.csv"))
        }
}



#2. Create transition mutation summary (ref) files of Overview2.2 (unfiltered) 
#2.1 Create merged metadata since some are missing from "Output_all/Ts_summary_metadata.1A.csv"

filename<-c("D75002","D75046", "D75007")
Metadata<-list()
for (i in 1:3){
        meta<-read.csv(paste0("Output_all/Overview",geno[i],"/Overview2.2/",filename[i],"-_overview2.2.csv"), stringsAsFactors = )
        meta<-meta[,c(2,3,7,9,22,34,35,38,41)]
        colnames(meta)[3:9]<-paste0(colnames(meta)[3:9],".",geno[i])
        Metadata[[i]]<-meta
        }

merged.meta<-do.call(cbind, Metadata)
merged.meta<-merged.meta[,-c(10,19)]

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

write.csv(merged.meta, "Output_all/Unfiltered/merged.metadata.csv")




#############################################################
merged.meta<-read.csv("Output_all/Unfiltered/merged.metadata.csv",stringsAsFactors = F,row.names = 1)
#############################################################

## 1) combine all "minor variant" frequency tables regardless of Maj==Ref or not
for (f in 1:3){
        flist<-list.files(paste0("Output_all/Overview",geno[f],"/Overview2.2/"),pattern="overview2.2.csv")
        
        MutFreq_Ts<-list()
        MutFreq_all<-list()
        MutFreq_Tvs<-list()
        
        for (i in 1:length(flist)){
                dat<-read.csv(paste0("Output_all/Overview",geno[f],"/Overview2.2/",flist[i]), stringsAsFactors = F, row.names=1)
                #filter out the sites with read depth <1000
                low_reads<-which(dat$TotalReads<1000) 
                dat[low_reads,c(7:ncol(dat))]<-NA
                
                #get the Ts. mutation frequency
                MutFreq_Ts[[i]]<-dat[,c("merged.pos","freq.Ts")] 
                MutFreq_Tvs[[i]]<-dat[,c("merged.pos","freq.transv")] 
                MutFreq_all[[i]]<-dat[,c("merged.pos","freq.mutations")] 
                
                filename<-substr(paste0(flist[i]),start=1,stop=7)
                names(MutFreq_Ts)[i]<-filename
                names(MutFreq_all)[i]<-filename
                names(MutFreq_Tvs)[i]<-filename
        }
        for (i in 1:length(MutFreq_Tvs)) {
                colnames(MutFreq_Ts[[i]])<-c("merged.pos",paste0(names(MutFreq_Ts[i])))
                colnames(MutFreq_Tvs[[i]])<-c("merged.pos",paste0(names(MutFreq_Tvs[i])))
                colnames(MutFreq_all[[i]])<-c("merged.pos",paste0(names(MutFreq_all[i])))
        }
        
        Ts<-MutFreq_Ts%>% purrr::reduce(full_join, by='merged.pos')
        Ts$mean<-rowMeans(Ts[2:ncol(Ts)],na.rm=T)
        
        Tvs<-MutFreq_Tvs%>% purrr::reduce(full_join, by='merged.pos')
        Tvs$mean<-rowMeans(Tvs[2:ncol(Tvs)],na.rm=T)
        
        all<-MutFreq_all%>% purrr::reduce(full_join, by='merged.pos')
        all$mean<-rowMeans(all[2:ncol(all)],na.rm=T)
        
        write.csv(Ts,paste0("Output_all/Unfiltered/Ts.MinorVarient_",geno[f],".csv"))
        write.csv(Tvs,paste0("Output_all/Unfiltered/Tvs.MinorVarient_",geno[f],".csv"))
        write.csv(all,paste0("Output_all/Unfiltered/All.MinorVarient_",geno[f],".csv"))
        
        
        s<-length(flist)
        Ts2<-Ts
        Ts2$sum<-apply(Ts2[2:(s+1)],1,function(x) sum(!is.na(x)))
        Ts2$keep0.3<-(Ts2$sum/s)>=1/3
        Ts2<-Ts2[Ts2$keep0.3==T,]
        
        Tvs2<-Tvs
        Tvs2$sum<-apply(Tvs2[2:(s+1)],1,function(x) sum(!is.na(x)))
        Tvs2$keep0.3 <-(Tvs2$sum/s)>=1/3
        Tvs2<-Tvs2[Tvs2$keep0.3==T,]
        
        All2<-all
        All2$sum<-apply(All2[2:(s+1)],1,function(x) sum(!is.na(x)))
        All2$keep0.3 <-(All2$sum/s)>=1/3
        All2<-All2[All2$keep0.3==T,]
        
        #what % are REMOVED?
        x<-1-sum(Ts2$keep0.3==T)/nrow(Ts) 
        cat("% removed by removing sites at least 1/3 samples in ", geno[f], " was ", x,"\n") 
         Ts2$mean[ Ts2$keep0.3==F]<-NA
        Tvs2$mean[Tvs2$keep0.3==F]<-NA
        All2$mean[All2$keep0.3==F]<-NA

         Ts2<- Ts2[,c("merged.pos","mean")]
        Tvs2<-Tvs2[,c("merged.pos","mean")]
        All2<-All2[,c("merged.pos","mean")]
        
        colnames( Ts2)[2]<-paste0("mean.",geno[f]) 
        colnames(Tvs2)[2]<-paste0("mean.",geno[f]) 
        colnames(All2)[2]<-paste0("mean.",geno[f]) 
        
        Ts2<-merge(merged.meta,Ts2,by="merged.pos", all.x = T)
        write.csv(Ts2,paste0("Output_all/Unfiltered/Ts.MinorVariant_Mean_metadata.",geno[f],".csv"))
        
        Tvs2<-merge(merged.meta,Tvs2,by="merged.pos",all.x = T)
        write.csv(Tvs2,paste0("Output_all/Unfiltered/Tvs.MinorVariant_Mean_metadata.",geno[f],".csv"))
        
        All2<-merge(merged.meta,All2,by="merged.pos",all.x = T)
        write.csv(All2,paste0("Output_all/Unfiltered/All.MinorVariant_Mean_metadata.",geno[f],".csv"))
        
}

#% removed by removing sites at least 1/3 samples in  1A  was  0.05545455 
#% removed by removing sites at least 1/3 samples in  1B  was  0.05886364
#% removed by removing sites at least 1/3 samples in  3A  was  0.05125

#Combine the mean mut freq of 3 genotypes into 1 table.

Summary<-read.csv(paste0("Output_all/Unfiltered/Ts.MinorVariant_Mean_metadata.1A.csv"),stringsAsFactors = F, row.names = 1)
nrow(Summary[!is.na(Summary$mean.1A),])
for (f in 2:3){
        d<-read.csv(paste0("Output_all/Unfiltered/Ts.MinorVariant_Mean_metadata.",geno[f],".csv"),stringsAsFactors = F)
        d<-d[,c(2,ncol(d))]
        print(geno[f])
        print(nrow(d[!is.na(d[,paste0("mean.",geno[f])]),]))
        Summary<-merge(Summary,d,by="merged.pos", all.x = T)
}
write.csv(Summary,"Output_all/Unfiltered/Ts.MinorVariant_Mean_3genotypes.csv")

#################################
merged.meta<-read.csv("Output_all/Unfiltered/merged.metadata.csv",stringsAsFactors = F,row.names = 1)

#2) separate the sites that are Maj!=Ref for comparing the same ancestral nucleotide between the genotypes

for (f in 1:3){
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
                
                dat.same<-dat[dat$MajNt==paste0(dat$ref),]
                dat.same<-dat.same[!is.na(dat.same$merged.pos),]
                dat.s<-merge(pos,dat.same, by="merged.pos", all.x = T)
                
                MutFreq_Ts.same[[i]]<-dat.s[,c("merged.pos","freq.Ts.ref")] 
                names(MutFreq_Ts.same)[i]<-filename
                
                MutFreq_Tvs.same[[i]]<-dat.s[,c("merged.pos","freq.transv.ref")] 
                MutFreq_all.same[[i]]<-dat.s[,c("merged.pos","freq.mutations.ref")] 
                names(MutFreq_Tvs.same)[i]<-filename
                names(MutFreq_all.same)[i]<-filename
                
                dat.diff<-dat[dat$MajNt!=dat$ref,]
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
        
        Ts.same<-MutFreq_Ts.same %>% purrr::reduce(full_join, by='merged.pos')
        Ts.same$mean<-rowMeans(Ts.same[2:ncol(Ts.same)],na.rm=T)
        write.csv(Ts.same,paste0("Output_all/Unfiltered/Ts.same_",geno[f],".csv"))
        
        Ts.diff<-MutFreq_Ts.diff %>% purrr::reduce(full_join, by='merged.pos')
        Ts.diff$mean<-rowMeans(Ts.diff[2:ncol(Ts.diff)],na.rm=T)
        write.csv(Ts.diff,paste0("Output_all/Unfiltered/Ts.diff_",geno[f],".csv"))
        
        Tvs.same<-MutFreq_Tvs.same %>% purrr::reduce(full_join, by='merged.pos')
        Tvs.same$mean<-rowMeans(Tvs.same[2:ncol(Tvs.same)],na.rm=T)
        write.csv(Tvs.same,paste0("Output_all/Unfiltered/Tvs.same_",geno[f],".csv"))
        
        All.same<-MutFreq_all.same %>% purrr::reduce(full_join, by='merged.pos')
        All.same$mean<-rowMeans(All.same[2:ncol(All.same)],na.rm=T)
        write.csv(Ts.same,paste0("Output_all/Unfiltered/All.same_",geno[f],".csv"))
        
        
        s<-length(flist)
        Ts2.same<-Ts.same
        Ts2.same$sum<-apply(Ts2.same[2:(s+1)],1,function(x) sum(!is.na(x)))
        Ts2.same$keep0.3<-(Ts2.same$sum/s)>=1/3
        Ts2.same<-Ts2.same[Ts2.same$keep0.3==T,]
        
        Tvs2.same<-Tvs.same
        Tvs2.same$sum<-apply(Tvs2.same[2:(s+1)],1,function(x) sum(!is.na(x)))
        Tvs2.same$keep0.3<-(Tvs2.same$sum/s)>=1/3
        Tvs2.same<-Tvs2.same[Tvs2.same$keep0.3==T,]
        
        All2.same<-All.same
        All2.same$sum<-apply(All2.same[2:(s+1)],1,function(x) sum(!is.na(x)))
        All2.same$keep0.3<-(All2.same$sum/s)>=1/3
        All2.same<-All2.same[All2.same$keep0.3==T,]
        
        #what % are REMOVED?
        x<-1-sum(Ts2.same$keep0.3==T)/nrow(Ts.same) 
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
        write.csv(Ts2.same,paste0("Output_all/Unfiltered/Ts.same_Mean_metadata.",geno[f],".csv"))
        
        Tvs2.same<-merge(merged.meta,Tvs2.same,by="merged.pos", all.x=T)
        write.csv(Tvs2.same,paste0("Output_all/Unfiltered/Tvs.same_Mean_metadata.",geno[f],".csv"))
        All2.same<-merge(merged.meta,All2.same,by="merged.pos", all.x=T)
        write.csv(All2.same,paste0("Output_all/Unfiltered/All.same_Mean_metadata.",geno[f],".csv"))
        
        Ts2.diff<-Ts.diff[,c("merged.pos","mean")]
        colnames(Ts2.diff)[2]<-paste0("mean.",geno[f]) 
        Ts2.diff<-merge(merged.meta,Ts2.diff,by="merged.pos", all.x=T)
        write.csv(Ts2.diff,paste0("Output_all/Unfiltered/Ts.diff_Mean_metadata.",geno[f],".csv"))
}

#% removed by removing sites at least 1/3 samples in  1A  was  0.04621401
#% removed by removing sites at least 1/3 samples in  1B  was  0.02192203
#% removed by removing sites at least 1/3 samples in  3A  was  0.01967058


#Combine the mean mut freq of 3 genotypes into 1 table.
Summary.same<-read.csv(paste0("Output_all/Unfiltered/Ts.same_Mean_metadata.1A.csv"),stringsAsFactors = F, row.names = 1)
nrow(Summary.same[!is.na(Summary.same$mean.1A),])
for (f in 2:3){
        d<-read.csv(paste0("Output_all/Unfiltered/Ts.same_Mean_metadata.",geno[f],".csv"),stringsAsFactors = F)
        d<-d[,c(2,ncol(d))]
        print(geno[f])
        print(nrow(d[!is.na(d[,paste0("mean.",geno[f])]),]))
        Summary.same<-merge(Summary.same,d,by="merged.pos", all.x = TRUE)
}
write.csv(Summary.same,"Output_all/Unfiltered/Ts.Same_Mean_3genotypes.csv")



Summary.diff<-read.csv(paste0("Output_all/Unfiltered/Ts.diff_Mean_metadata.1A.csv"),stringsAsFactors = F)
Summary.diff<-Summary.diff[,-1]
nrow(Summary.diff[!is.na(Summary.diff$mean.1A),])
for (f in 2:3){
        d<-read.csv(paste0("Output_all/Unfiltered/Ts.diff_Mean_metadata.",geno[f],".csv"),stringsAsFactors = F)
        d<-d[,c(2,ncol(d))]
        print(geno[f])
        print(nrow(d[!is.na(d[,paste0("mean.",geno[f])]),]))
        Summary.diff<-merge(Summary.diff,d,by="merged.pos", all.x = T)
}


write.csv(Summary.diff,"Output_all/Unfiltered/Ts.Diff_Mean_3genotypes.csv")


Summary.sameTv<-read.csv(paste0("Output_all/Unfiltered/Tvs.same_Mean_metadata.1A.csv"),stringsAsFactors = F, row.names = 1)
for (f in 2:3){
        d<-read.csv(paste0("Output_all/Unfiltered/Tvs.same_Mean_metadata.",geno[f],".csv"),stringsAsFactors = F)
        d<-d[,c(2,ncol(d))]
        Summary.sameTv<-merge(Summary.sameTv,d,by="merged.pos", all.x = T)
}
write.csv(Summary.same,"Output_all/Unfiltered/Tvs.Same_Mean_3genotypes.csv")

Summary.All<-read.csv(paste0("Output_all/Unfiltered/All.same_Mean_metadata.1A.csv"),stringsAsFactors = F, row.names = 1)
for (f in 2:3){
        d<-read.csv(paste0("Output_all/Unfiltered/All.same_Mean_metadata.",geno[f],".csv"),stringsAsFactors = F)
        d<-d[,c(2,ncol(d))]
        Summary.All<-merge(Summary.All,d,by="merged.pos", all.x = T)
}
write.csv(Summary.All,"Output_all/Unfiltered/All.Same_Mean_3genotypes.csv")



#############################################################
#############################################################

## **START** ##
Summary<-read.csv("Output_all/Unfiltered/Ts.MinorVariant_Mean_3genotypes.csv",stringsAsFactors = F, row.names = 1)
Summary$gene[Summary$gene=="NS1(P7)"]<-"NS1"

## Plot across genomes all together (Figure 1)
Summary2<-Summary[Summary$merged.pos<=8500 & Summary$merged.pos>=264,]

#colors2<-brewer.pal(n = 8, name = 'Set2')
#colors2<-paste0(cols2,"99")
#brewer.pal(n = 8, name = "Dark2")
#colors2<-qualitative_hcl(6, palette="Dark3")
#hcl_palettes(plot = TRUE)


genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
genenames<-genes$Gene[1:12]

Summary2$roll1001a<-c(rep(NA, times=99), rollmean(Summary2$mean.1A, k=100, na.rm=T))
Summary2$roll1001b<-c(rep(NA, times=99), rollmean(Summary2$mean.1B, k=100, na.rm=T))
Summary2$roll1003a<-c(rep(NA, times=99), rollmean(Summary2$mean.3A, k=100, na.rm=T))


pdf(paste0("Output_all/MVfreq_AcrossGenome_colorspace-dark3-3.pdf"),width=12,height=4.5)

plot(mean.1A~merged.pos, data=Summary2,t="n",log='y',yaxt='n',xlab='Genome position', ylab="Minor variant frequency",
     ylim=c(10^-4,0.4),xlim=c(264,8500))
eaxis(side = 2, at = 10^((-0):(-(4))), cex=2)
for(k in 1:4){abline(h = 1:10 * 10^(-k), col = "gray90",lty=1,lwd=.5)}

points(mean.1B~merged.pos,pch=16, data=Summary2, col=paste0(colors2[3],"B3"), cex=0.2)
points(mean.3A~merged.pos,pch=16, data=Summary2, col=paste0(colors2[5],"B3"), cex=0.2)
points(mean.1A~merged.pos,pch=16, data=Summary2, col=paste0(colors2[1],"B3"), cex=0.2)

lines(roll1001b~merged.pos, data=Summary2, col=colors2[3],lwd=1)
lines(roll1003a~merged.pos, data=Summary2, col=colors2[5],lwd=1)
lines(roll1001a~merged.pos, data=Summary2, col=colors2[1],lwd=1)

legend("topright",legend=c("1A","1B","3A"), col=colors2[c(1,3,5)], pch=16, bty = "n", cex=1)

genes2<-genes[1:12,]
#genes2$Gene[6]<-"NS1"

abline(v=genes2$end, col="gray80", lwd=.5)

ylow<-0.0001;yhigh<-.4
for (j in 1:nrow(genes2)){
        xleft<- genes2$start[j]
        xright<-genes2$start[j+1]
        if (j==1){
                
                rect(-200,ylow,xright,ylow*1.6,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(150,1.27*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else if (j==6){
                rect(xleft,ylow,xright,ylow*1.6,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xleft+100,1.27*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
                
        }
        else if (j==4|j==9){
                rect(xleft,ylow,xright,1.6*ylow,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xleft+50,2.1*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else if (j==12){
                rect(xleft,ylow,genes2$end[j],1.6*ylow,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xleft+600,1.27*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else{rect(xleft,ylow,xright,1.6*ylow,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xright-(xright-xleft)/2,1.27*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)}
}

box()
dev.off()


####


#Tvs
Sum2<-read.csv(paste0("Output_all/Unfiltered/Tvs.MinorVariant_Mean_metadata.1A.csv"),stringsAsFactors = F, row.names = 1)
colnames(Sum2)[ncol(Sum2)]<-paste0("mean.tvs.1A")
for (f in 2:3){
        d<-read.csv(paste0("Output_all/Unfiltered/Tvs.MinorVariant_Mean_metadata.",geno[f],".csv"),stringsAsFactors = F)
        d<-d[,c(2,ncol(d))]
        colnames(d)[2]<-paste0("mean.tvs.",geno[f])
        Sum2<-merge(Sum2,d,by="merged.pos", all.x = T)
}
write.csv(Sum2,"Output_all/Unfiltered/Tvs.MinorVariant_Mean_3genotypes.csv")

#all
Sum3<-read.csv(paste0("Output_all/Unfiltered/All.MinorVariant_Mean_metadata.1A.csv"),stringsAsFactors = F, row.names = 1)
colnames(Sum3)[ncol(Sum3)]<-paste0("mean.all.1A")
for (f in 2:3){
        d<-read.csv(paste0("Output_all/Unfiltered/All.MinorVariant_Mean_metadata.",geno[f],".csv"),stringsAsFactors = F)
        d<-d[,c(2,ncol(d))]
        colnames(d)[2]<-paste0("mean.all.",geno[f])
        Sum3<-merge(Sum3,d,by="merged.pos", all.x = T)
}
write.csv(Sum3,"Output_all/Unfiltered/Tvs.MinorVariant_Mean_3genotypes.csv")

SumT1<-Summary[,c("merged.pos", "mean.1A","mean.1B", "mean.3A")]
SumT<-Sum2[,c("merged.pos","mean.tvs.1A", "mean.tvs.1B", "mean.tvs.3A")]
SumT2<-Sum3[,c("merged.pos","mean.all.1A", "mean.all.1B", "mean.all.3A")]

SumT<-merge(SumT1, SumT, by="merged.pos")
SumT<-merge(SumT, SumT2, by="merged.pos")



tb<-data.frame(type=colnames(SumT)[2:10])
tb$Genotype<-rep(c("1A","1B","3A"),times=3)
tb$Type<-c(rep("Transition",times=3), rep("Transversion",times=3), rep("Total MV",times=3))
for (i in 1:nrow(tb)){
        n<-i+1
        tb$Mean[i]<-mean(SumT[,n], na.rm=T)
        tb$SE[i]<-std.error(SumT[,n], na.rm=T)
}
write.csv(tb, "Output_all/Unfiltered/MVFreq_Summary_3genotypes.csv")


ggplot(tb,aes(x=Type,y=Mean,group=Genotype, color=Genotype))+
        geom_point(position=position_dodge(width=0.3),size =1.5)+scale_color_manual(values=colors2[c(1,3,5)])+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, position=position_dodge(width=0.3))+
        theme(axis.title.x=element_blank())+ylab("Mean minor variant frequency")+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=13),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:2)+0.5,  
                   color = "gray60", size=.5)
        
ggsave(filename="Output_all/Ave.Mutfreq_by.type_by.genotypes.pdf",width = 5, height = 5.8)



#barplot
ggplot(tb,aes(x=Type,y=Mean, ymin=Mean-SE, ymax=Mean+SE, fill=Genotype))+
        geom_bar(position=position_dodge(.9), stat="identity",width=0.8)+
        scale_fill_manual(values=paste0(colors2[c(1,3,5)],"E6"), labels=c("1A","1B","3A"))+
        geom_errorbar(position=position_dodge(.9), width=.2, color="gray30")+
        theme(axis.title.x=element_blank())+ylab("Mean minor variant frequency")+
        theme_linedraw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        geom_vline(xintercept = c(1:2)+0.5,  
                   color = "gray60", size=.4)+
        theme(panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(linetype=2, colour="gray60"),
              panel.grid.minor.y = element_blank())
ggsave(filename="Output_all/Ave.Mutfreq_by.type_by.genotypes_barplot2.pdf",width = 6, height = 4.5)



#boxplot

tb2<-SumT[,-1]
tb2m<-melt(tb2)
tb2m$Genotype<-c(rep("1A", times=8800),rep("1B", times=8800),rep("3A", times=8800),rep("1A", times=8800),rep("1B", times=8800),rep("3A", times=8800),rep("1A", times=8800),rep("1B", times=8800),rep("3A", times=8800))
tb2m$Type<-c(rep("Transitions", times=8800*3), rep("Transversions", times=8800*3), rep("Totals", times=8800*3))
tb2m<-tb2m[!is.na(tb2m$value),]

label_scientific <- function(l) {
        # turn in to character string in scientific notation
        l <- format(l, scientific = TRUE)
        # quote the part before the exponent to keep all the digits
        l <- gsub("^(.*)e", "'\\1'e", l)
        # turn the 'e+' into plotmath format
        l <- gsub("e", "%*%10^", l)
        # return this as an expression
        parse(text=l)
}

col80<- paste0(colors2,"CC")
col80<-col80[c(1,3,5)]


ggplot(tb2m,aes(x=Type, y=value,fill=Genotype))+scale_y_continuous(trans = 'log10', labels=label_scientific)+
        scale_fill_manual(values=paste0(colors2[c(1,3,5)],"CC"), labels=c("1A","1B","3A"))+
        geom_boxplot(aes(middle=mean(value), color=Genotype),outlier.alpha = 0.2)+ 
        labs(x="")+ylab("Mean minor variant frequency")+
        theme_bw()+
        theme(axis.text.x = element_text(size=13),axis.title.y = element_text(size=13))
        
ggsave(filename="Output_all/Ave.Mutfreq_by.type_by.genotypes_boxplot.pdf",width = 7, height = 7)



####
## plot mean and SE for each gene
mf<-data.frame()
#geno<-c("1A","1B","3A")
for(g in 1:3){
        dat<-Summary[,c("merged.pos","gene",paste0("mean.",geno[g]))]
        colnames(dat)[3]<-"mean"
        dat$genotype<-geno[g]
        mf<-rbind(mf,dat)
        }


mf1<-mf[!is.na(mf$mean),]
SumMFGenes<-aggregate(mf1$mean,by=list(mf1$genotype,mf1$gene),FUN=mean)
SumMF.se<-aggregate(mf1$mean,by=list(mf1$genotype,mf1$gene),FUN=std.error)

sumG<-cbind(SumMFGenes, SumMF.se$x)
colnames(sumG)<-c("Genotype","Gene","Mean","SE")
sumG$Gene<-factor(sumG$Gene, levels=c("5' UTR","Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))


ggplot(sumG, aes(x=Gene, y=Mean, group=Genotype, color=Genotype))+
        geom_point(position=position_dodge(width=0.3))+scale_fill_manual(values=colors2[c(1,3,5)])+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, position=position_dodge(width=0.3))+
        theme_bw()+theme(axis.text=element_text(size=11), axis.title=element_text(size=13))+ylab("Minor variant frequency")+
        theme(panel.grid.major.x=element_blank(),axis.title.y = element_text(size=13))+
        geom_vline(xintercept = c(1:11)+0.5,  
                   color = "gray70", size=.5)+
        theme(axis.title.x=element_blank())+
        ylim(0.001,0.022)
#ggsave(filename="Output_all/Unfiltered/Ave.MVfreq_by.gene_by.genotype.pdf",width = 10, height = 8)
ggsave(filename="Output_all/Ave.MVfreq_by.gene_by.genotype.pdf", width = 8.5, height = 5)



ggplot(sumG,aes(x=Gene,y=Mean, ymin=Mean-SE, ymax=Mean+SE, fill=Genotype))+
        geom_bar(position=position_dodge(.9), stat="identity",width=0.8)+
        scale_fill_manual(values=paste0(colors2[c(1,3,5)],"E6"), labels=c("1A","1B","3A"))+
        geom_errorbar(position=position_dodge(.9), width=.2, color="gray40")+
        theme(axis.title.x=element_blank())+ylab("Mean minor variant frequency")+
        theme_linedraw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        geom_vline(xintercept = c(1:11)+0.5,  
                   color = "gray60", size=.4)+
        theme(panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(linetype=2, colour="gray60"),
              panel.grid.minor.y = element_blank())

ggsave(filename="Output_all/Ave.MVfreq_by.gene_by.genotype_Barplot.pdf", width = 8.5, height = 5)



### Plot MVF separately for Nonsyn and syn

mf2<-data.frame()
for(g in 1:3){
        dat<-Summary[,c("merged.pos","gene",paste0("mean.",geno[g]),paste0("Type.",geno[g]))]
        colnames(dat)[3:4]<-c("mean", "Type")
        dat$genotype<-geno[g]
        mf2<-rbind(mf2,dat)
}


mf2<-mf2[!is.na(mf2$mean),]
#remove the stop sites
mfstop<-mf2[mf2$Type=="stop",]
mf3<-mf2
#plot(mfstop$mean[mfstop$genotype=="1A"])
mf2<-mf2[mf2$Type!="stop",]
#remove 5'UTR
mf2<-mf2[mf2$gene!="5' UTR",]
SumMFGenes2<-aggregate(mf2$mean,by=list(mf2$genotype,mf2$gene, mf2$Type),FUN=mean)
SumMF.se2<-aggregate(mf2$mean,by=list(mf2$genotype,mf2$gene, mf2$Type),FUN=std.error)

sumG2<-cbind(SumMFGenes2, SumMF.se2$x)
colnames(sumG2)<-c("Genotype","Gene","Type","Mean","SE")
sumG2$Gene<-factor(sumG2$Gene, levels=c("Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))


ggplot(sumG2, aes(x=Gene, y=Mean, shape=Type,color=Genotype))+
        geom_point(position=position_dodge(width=0.8),size =2)+scale_color_manual(values=colors2[c(1,3,5)])+
        scale_shape_manual(values=c(4,19))+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, size=.2, position=position_dodge(width=0.8))+
        theme(axis.title.x=element_blank())+ylab("Minor variant frequency")+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray60", size=.4)+
        guides(
                shape = guide_legend(order = 2),
                color=guide_legend(order=1)
        )

ggsave(filename="Output_all/Ave.MVfreq_by.gene_by.genotype_N.NS.pdf", width = 8, height = 5)


#### include nonsense (stop) mutations
#remove 5' UTR
mf3.2<-mf3[mf3$gene!="5' UTR",]
SumMFGenes3<-aggregate(mf3.2$mean,by=list(mf3.2$genotype,mf3.2$gene, mf3.2$Type),FUN=mean)
SumMF.se3<-aggregate(mf3.2$mean,by=list(mf3.2$genotype,mf3.2$gene, mf3.2$Type),FUN=std.error)

sumG3<-cbind(SumMFGenes3, SumMF.se3$x)
colnames(sumG3)<-c("Genotype","Gene","Type","Mean","SE")
sumG3$Gene<-factor(sumG3$Gene, levels=c("5' UTR","Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))


ggplot(sumG3, aes(x=Gene, y=Mean, shape=Type,color=Genotype))+
        geom_point(position=position_dodge(width=0.8),size =2)+scale_color_manual(values=colors2[c(1,3,5)])+
        scale_shape_manual(values=c(4,17,19))+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, size=.2, position=position_dodge(width=0.8))+
        theme(axis.title.x=element_blank())+ylab("Minor variant frequency")+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray60", size=.4)+
        guides(
                shape = guide_legend(order = 2),
                color=guide_legend(order=1)
        )

ggsave(filename="Output_all/Ave.MVfreq_by.gene_by.genotype_N.NS.stop.pdf", width = 8.5, height = 5)








  

###
#plot mean https://stackoverflow.com/questions/19876505/boxplot-show-the-value-of-mean   not finished
means <- aggregate(value ~ Type+Genotype, tb2m, mean)

min.mean.sd.max <- function(x) {
        r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
        names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
        r
}


  
ggplot(tb2m,aes(x=Type, y=value,fill=Genotype))+scale_y_continuous(trans = 'log10', labels=label_scientific)+
        scale_fill_manual(values=col80, labels=c("1A","1B","3A"))+
        stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") 

        geom_boxplot(aes(middle=mean(value), color=Genotype),outlier.alpha = 0.2)+ 
        labs(x="")+ylab("Mean minor variant frequency")+
        theme_bw()+
        theme(axis.text.x = element_text(size=14))+
        stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                     width = .75, linetype = "dashed")

geom_boxplot(
        aes(
                x = wind, ymin = min(dust), lower = quantile(dust, .25), 
                middle = mean(dust), upper = quantile(dust, .75),
                ymax = max(dust)
        ), stat = "identity"
)



### 
# stat test
Comb<-t(combn(1:3, 2))
combnames<-c("1A-1B","1A-3A","1B-3A")

results.Ts<-list()
results.Tvs<-list()
results.all<-list()

for (i in 1:3){
        k=1
        n1<-Comb[i,1]+k
        n2<-Comb[i,2]+k
        r1<-wilcox.test(SumT[,n1], SumT[,n2], alternative = "greater", paired = FALSE) 
        results.Ts[[i]]<-r1
        names(results.Ts)[i]<-combnames[i]
        k=4
        n1<-Comb[i,1]+k
        n2<-Comb[i,2]+k
        r2<-wilcox.test(SumT[,n1], SumT[,n2], alternative = "greater", paired = FALSE) 
        results.Tvs[[i]]<-r2
        names(results.Tvs)[i]<-combnames[i]
        k=7
        n1<-Comb[i,1]+k
        n2<-Comb[i,2]+k
        r3<-wilcox.test(SumT[,n1], SumT[,n2], alternative = "greater", paired = FALSE) 
        results.all[[i]]<-r3
        names(results.all)[i]<-combnames[i]
 
}

sink("Output_all/MutFreq_mean_comparison_genotypes.txt")
for (i in 1:3){
        print(names(results.Ts[i]))
        print("Transition")
        print(results.Ts[[i]][3])
        print("Transversion")
        print(results.Tvs[[i]][3])
        print("Total")
        print(results.all[[i]][3])
}
sink(NULL)
#[1] "1A-1B"
#[1] "Transition"
#$p.value
#[1] 4.623662e-18
#
#[1] "Transversion"
#$p.value
#[1] 2.191749e-30
#
#[1] "Total"
#$p.value
#[1] 3.892474e-22
#
#[1] "1A-3A"
#[1] "Transition"
#$p.value
#[1] 5.445805e-23
#
#[1] "Transversion"
#$p.value
#[1] 7.819815e-104
#
#[1] "Total"
#$p.value
#[1] 2.285575e-37
#
#[1] "1B-3A"
#[1] "Transition"
#$p.value
#[1] 0.1606904
#
#[1] "Transversion"
#$p.value
#[1] 1.639254e-30
#
#[1] "Total"
#$p.value
#[1] 0.001043338


sumTs<-Summary[,c("merged.pos", "gene", "mean.1A","mean.1B", "mean.3A")]
genenames<-unique(genes$Gene)
Comb<-t(combn(1:3, 2))
combnames<-c("1A-1B","1A-3A","1B-3A")

wilcox.res<-data.frame(matrix(nrow=12, ncol=3))
rownames(wilcox.res)<-genenames[1:12]
colnames(wilcox.res)<-combnames                       
for (i in 1:3){
        k=2
        n1<-Comb[i,1]+k
        n2<-Comb[i,2]+k
        for (g in 1:12){
                r1<-wilcox.test(sumTs[sumTs$gene==genenames[g],n1], sumTs[sumTs$gene==genenames[g],n2], alternative = "greater", paired = FALSE) 
                wilcox.res[g,i]<-r1[[3]][1]
        }
}

write.csv(wilcox.res,"Output_all/MF_wilcoxTest_byGene.csv",row.names = T)


### 
# stat test for comparing nonsyn vs syn MVF 
genenames<-unique(Summary$gene)

sumTs<-Summary[,c("merged.pos", "gene", "mean.1A","mean.1B", "mean.3A","Type.1A","Type.1B","Type.3A")]

wilcox.resS<-data.frame(matrix(nrow=12, ncol=3))
wilcox.resNS<-data.frame(matrix(nrow=12, ncol=3))
wilcox.res2<-data.frame(matrix(nrow=12, ncol=3))
rownames(wilcox.resS)<-genenames
colnames(wilcox.resS)<- combnames                 
rownames(wilcox.resNS)<-genenames
colnames(wilcox.resNS)<- combnames                   
rownames(wilcox.res2)<-genenames
colnames(wilcox.res2)<- geno                    

for (i in 1:3){
        d<-Summary[,c("merged.pos", "gene", paste0("mean.",geno[i]),paste0("Type.",geno[i]))]
        colnames(d)[3:4]<-c("mean", "Type")
        
        k=2
        n1<-Comb[i,1]+k
        n2<-Comb[i,2]+k
        
        if (i==1) d1<-sumTs[sumTs$Type.1A==sumTs$Type.1B,]
        if (i==2) d1<-sumTs[sumTs$Type.1A==sumTs$Type.3A,]
        if (i==3) d1<-sumTs[sumTs$Type.1B==sumTs$Type.3A,]
        
        
        for (g in 1:12){
                r1s<-wilcox.test(d1[d1$gene==genenames[g]& d1$Type.1A=="syn",n1], d1[d1$gene==genenames[g] & d1$Type.1A=="syn",n2], alternative = "greater", paired = FALSE) 
                wilcox.resS[g,i]<-r1s[[3]][1]
                r1n<-wilcox.test(d1[d1$gene==genenames[g]& d1$Type.1A=="nonsyn",n1], d1[d1$gene==genenames[g] & d1$Type.1A=="nonsyn",n2], alternative = "greater", paired = FALSE) 
                wilcox.resNS[g,i]<-r1n[[3]][1]
                
                r2<-wilcox.test(d$mean[d$gene==genenames[g]& d$Type=="syn"], d$mean[d$gene==genenames[g] & d$Type=="nonsyn"], alternative = "greater", paired = FALSE) 
                wilcox.res2[g,i]<-r2[[3]][1]
                
                
                
        }
}

write.csv(wilcox.resS,"Output_all/MVF_wilcoxTest_byGene.Syn.csv",row.names = T)
write.csv(wilcox.resNS,"Output_all/MVF_wilcoxTest_byGene.NS.csv",row.names = T)
write.csv(wilcox.res2,"Output_all/MVF_wilcoxTest_NSvsS.csv",row.names = T)


