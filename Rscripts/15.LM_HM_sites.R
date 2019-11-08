library(purrr)
library(tidyverse)
library(zoo)
library(colorspace)
source("Rscripts/label_scientific.R")

colors2<-qualitative_hcl(6, palette="Dark3")


source("Rscripts/baseRscript.R")

mfs<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F, row.names = 1)
mfs<-mfs[mfs$pos>351,]

genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector

genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
genenames<-genes$Gene[1:12]

### Plot mutation freq. across the genome based on the mutation types 

n<-data.frame("pos"=c(1:mfs$pos[nrow(mfs)]))
mfs<-merge(n,mfs,by="pos",all.x=T)
mfs2<-mfs[,c(1,197:204)]
hmlm<-read.csv("Data/HMLMsites.csv")

HMLM<-list()
MF<-list()
for (i in 1:nrow(hmlm)){
        df<-mfs2[mfs2$pos>=hmlm$Start[i]&mfs2$pos<=hmlm$End[i],]
        HMLM[[i]]<-df
        names(HMLM)[i]<-i
        MF[[i]]<-mean(df$mean, na.rm = T)
        names(MF)[i]<-paste0(hmlm$Type[i],i)
}

MFave<-as.data.frame(do.call(rbind,MF))
MFave$ID<-rownames(MFave)
MFave$ID[14:19]<-paste0("HM",1:6)
colnames(MFave)[1]<-"MF"
allmean<-mean(mfs$mean,na.rm = T)
ggplot(MFave, aes(x=ID, y=MF))+
        geom_point()+
        theme_bw()+
        ylab("Average mutation frequency")
        theme(axis.text.x=element_text(angle=90, hjust=1))+
        geom_hline(yintercept=allmean, color="gray60", linetype=2)


   
#######     
#Use all mutation freq
#
mfs<-read.csv("Output1A/MutFreq.filtered/Filtered.AllMutFreq.Q35.csv",stringsAsFactors = F, row.names = 1)
mfs<-mfs[mfs$pos>351,]

genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector

genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
genenames<-genes$Gene[1:12]

### Plot mutation freq. across the genome based on the mutation types 

n<-data.frame("pos"=c(1:mfs$pos[nrow(mfs)]))
mfs<-merge(n,mfs,by="pos",all.x=T)
mfs2<-mfs[,c(1,197:203)]
hmlm<-read.csv("Data/HMLMsites.csv")

HMLM<-list()
MF<-list()
for (i in 1:nrow(hmlm)){
        df<-mfs2[mfs2$pos>=hmlm$Start[i]&mfs2$pos<=hmlm$End[i],]
        HMLM[[i]]<-df
        names(HMLM)[i]<-i
        MF[[i]]<-mean(df$mean, na.rm = T)
        names(MF)[i]<-paste0(hmlm$Type[i],i)
}

MFave<-as.data.frame(do.call(rbind,MF))
MFave$ID<-rownames(MFave)
MFave$ID[14:19]<-paste0("HM",1:6)
colnames(MFave)[1]<-"MF"
allmean<-mean(mfs$mean,na.rm = T)
ggplot(MFave, aes(x=ID, y=MF))+
        geom_point()+
        theme_bw()+
        ylab("Average mutation frequency")+
        theme(axis.text.x=element_text(angle=90, hjust=1))+
        geom_hline(yintercept=allmean, color="gray60", linetype=2)+
        theme(axis.title.x=element_blank())
ggsave("Output1A/SummaryFig.Filtered/HM.LM.sites.pdf", width = 5, height = 4.5)
#no patterns compared to in vitro LM HM sites
