library(tidyverse)
library(zoo)
library(ggplot2)
library(reshape)
library(ggpubr)
library(ggthemes)
library(sfsmisc)

########################
filtered<-list()
fnames<-c("Ts", "Ts_NA", "Ts_zero")
for (i in 1:3){
        filename<-fnames[i]
        df<-read.csv(paste0("Output/Mut.freq.filtered/Summary_",filename,".Q35.csv"),stringsAsFactors = F)
        df<-df[,-1]
        filtered[[i]]<-df
        names(filtered)[i]<-filename
}


#Plot summary of  CpG creating vs. non-CpGcreating mutation frequencies   



for (j in 1:length(filtered)){
  TS<-filtered[[j]]
  
  #syn
  k=1
  transMF<-list()
  for (i in c("a","t")) {
          for (cp in c(0,1)){
                  datavector<-TS$mean[TS$Type=="syn" & TS$ref==i &TS$makesCpG==cp]
                  nt<-toupper(i)
                  if (cp==0) cpg<-''
                  if (cp==1) cpg<-"(CpG)"
                  vname<-paste0(nt,cpg)
                  dat<-data.frame(base=rep(vname,times=length(datavector)), mf=datavector)
                  transMF[[k]]<-dat
                  names(transMF)[k]<-vname
                  k=k+1
          }
  }
  
  mfdata1<-do.call(rbind, transMF)
  
  z=c(0.7,0.3,0.7,0.3)
  mfplot1<-ggplot(mfdata1,aes(x=base,y=mf,fill=base))+geom_boxplot(alpha=z)+labs(x="Nucleotide",y="Mutation frequency")+
          scale_fill_manual(values=c("#66CCEEB3","#66CCEE4D","#228833B3","#2288334D")) + theme_classic()+
          guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
          theme(axis.text.y = element_text(size =10))+
          ggtitle("CpG vs. nonCpG syn transition") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
          theme(plot.title = element_text(hjust = 0.5))+
          geom_vline(xintercept = 2.5)
  
  ggsave(filename=paste0("Output/SummaryFig.Filtered/CpGvsNonCpG_syn_",names(filtered)[j], ".Q35.pdf"),width=4, height=4, units='in',device='pdf', plot=mfplot1)
  
  #nonsyn
  k=1
  transMF2<-list()
  for (i in c("a","t")) {
          for (cp in c(0,1)){
                  datavector<-TS$mean[TS$Type=="nonsyn" & TS$ref==i &TS$makesCpG==cp]
                  nt<-toupper(i)
                  if (cp==0) cpg<-''
                  if (cp==1) cpg<-"(CpG)"
                  vname<-paste0(nt,cpg)
                  dat<-data.frame(base=rep(vname,times=length(datavector)), mf=datavector)
                  transMF2[[k]]<-dat
                  names(transMF2)[k]<-vname
                  k=k+1
          }
  }
  
  mfdata2<-do.call(rbind, transMF2)
  
  z=c(0.7,0.3,0.7,0.3)
  mfplot2<-ggplot(mfdata2,aes(x=base,y=mf,fill=base))+geom_boxplot(alpha=z)+labs(x="Nucleotide",y="Mutatin frequency")+
          scale_fill_manual(values=c("#66CCEEB3","#66CCEE4D","#228833B3","#2288334D")) + theme_classic()+
          guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
          theme(axis.text.y = element_text(size =10))+
          ggtitle("CpG vs. nonCpG nonsyn transition") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
          theme(plot.title = element_text(hjust = 0.5))+
          geom_vline(xintercept = 2.5)
  
  ggsave(filename=paste0("Output/SummaryFig.Filtered/CpGvsNonCpG_nonsyn_",names(filtered)[j], ".Q35.pdf"), width=4, height=4, units='in',device='pdf', plot=mfplot2)
}
  
  
  





#####################################





###
# 2.Transversion syn     

k=1
tvMF<-list()
for (i in c("a","t","c","g")) {
        for (cp in c(0,1)){
                datavector1<-Tv1$mean[Tv1$Type.tv1=="syn" & Tv1$ref==i & Tv1$makesCpG.tv1==cp]
                datavector2<-Tv2$mean[Tv2$Type.tv2=="syn" & Tv2$ref==i & Tv2$makesCpG.tv2==cp]
                datavec<-c(datavector1,datavector2)
                
                nt<-toupper(i)
                if (cp==0) cpg<-''
                if (cp==1) cpg<-"(CpG)"
                vname<-paste0(nt,cpg)
                dat<-data.frame(base=rep(vname,times=length(datavec)), mf=datavec)
                tvMF[[k]]<-dat
                names(tvMF)[k]<-vname
                k=k+1
        }
}
mftvdata<-do.call(rbind, tvMF)

z=c(0.7,0.3,0.7,0.3,0.7,0.3,0.7,0.3)
mfplot3<-ggplot(mftvdata,aes(x=base,y=mf,fill=base))+geom_boxplot(alpha=z)+labs(x="Nucleotide",y="Mutation frequency")+
        scale_fill_manual(values=c("#66CCEE","#66CCEE","#228833","#228833","#CCBB44","#CCBB44","#EE6677","#EE6677")) + theme_classic()+
        guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10))+
        ggtitle("CpG vs. nonCpG syn transversion") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
        theme(plot.title = element_text(hjust = 0.5))+
        geom_vline(xintercept = c(2.5,4.5,6.5), color="gray60")

ggsave(filename="Output/SummaryFig.Filtered/MF_CpGvsNonCpG_syn_Transversion.pdf",width=6, height=4, units='in',device='pdf', plot=mfplot3)

#nonsyn transversion
k=1
tvMF2<-list()
for (i in c("a","t","c","g")) {
        for (cp in c(0,1)){
                datavector1<-Tv1$mean[Tv1$Type.tv1=="nonsyn" & Tv1$ref==i &Tv1$makesCpG.tv1==cp]
                datavector2<-Tv2$mean[Tv2$Type.tv2=="nonsyn" & Tv2$ref==i &Tv2$makesCpG.tv2==cp]
                datavec<-c(datavector1,datavector2)
                
                nt<-toupper(i)
                if (cp==0) cpg<-''
                if (cp==1) cpg<-"(CpG)"
                vname<-paste0(nt,cpg)
                dat<-data.frame(base=rep(vname,times=length(datavec)), mf=datavec)
                tvMF2[[k]]<-dat
                names(tvMF2)[k]<-vname
                k=k+1
        }
}
mftvdata2<-do.call(rbind, tvMF2)

z=c(0.7,0.3,0.7,0.3,0.7,0.3,0.7,0.3)
mfplot4<-ggplot(mftvdata2,aes(x=base,y=mf,fill=base))+geom_boxplot(alpha=z)+labs(x="Nucleotide",y="Mutation frequency")+
        scale_fill_manual(values=c("#66CCEE","#66CCEE","#228833","#228833","#CCBB44","#CCBB44","#EE6677","#EE6677")) + theme_classic()+
        guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10))+
        ggtitle("CpG vs. nonCpG nonsyn transversion") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
        theme(plot.title = element_text(hjust = 0.5))+
        geom_vline(xintercept = c(2.5,4.5,6.5), color="gray60")

ggsave(filename="Output/SummaryFig.Filtered/MF_CpGvsNonCpG_nonsyn_Transversion.pdf",width=6, height=4, units='in',device='pdf', plot=mfplot4)



