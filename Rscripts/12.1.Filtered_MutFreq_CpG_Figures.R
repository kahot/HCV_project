library(tidyverse)
library(zoo)
library(ggplot2)
library(reshape)
library(ggpubr)
library(ggthemes)
library(sfsmisc)
source("Rscripts/baseRscript.R")
#Plot summary of  CpG creating vs. non-CpGcreating mutation frequencies   

TS<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F, row.names = 1)
colnames(TS)



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
mfplot1
ggsave(filename="Output1A/SummaryFig.Filtered/CpGvsNonCpG_syn_Ts.Q35.pdf",width=4, height=4, units='in',device='pdf', plot=mfplot1)

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
ggsave(filename="Output/SummaryFig.Filtered/CpGvsNonCpG_nonsyn_Ts.Q35.pdf", width=4, height=4, units='in',device='pdf', plot=mfplot2)


###################  

mf<-data.frame()

dat<-Summary[,c("merged.pos","gene",paste0("mean.",geno[g]))]
        colnames(dat)[3]<-"mean"
        dat$genotype<-geno[g]
        mf<-rbind(mf,dat)


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


####
 
TS$SE<-apply(TS[,2:196],1,function(x) std.error(x,na.rm = T))
TS$SD<-apply(TS[,2:196],1,function(x) sd(x,na.rm = T))

plot(TS$SE, pch=".")
plot(TS$SD, pch=".")

ts2<-TS[,c("pos","mean","SE","SD")]
ts2$CV<-ts2$SD/ts2$mean

ggplot(ts2, aes(x=pos, y=CV))+
        geom_point(size=.2)+ylab("CV")+
        theme_bw()+
        labs(x="")

ggplot(ts2, aes(x=pos, y=CV))+
        geom_point()+
        geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=.2, size=.2)+
        theme(axis.title.x=element_blank())+ylab("Mutation frequency")+
        theme_bw()+
        labs(x="")
        
        
        +
        guides(
                shape = guide_legend(order = 2),
                color=guide_legend(order=1)
        )

