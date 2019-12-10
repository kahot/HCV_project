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

