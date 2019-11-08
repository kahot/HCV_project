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

library(colorspace)
source("Rscripts/label_sceintific.R")

#dir.create("Output_all/Fst_T")
colors2<-qualitative_hcl(6, palette="Dark3")
col2_light<-qualitative_hcl(6, palette="Set3")
div.colors<-c(colors2[1],col2_light[1],colors2[3],col2_light[3],colors2[5],col2_light[5])
color.genes<-qualitative_hcl(11, palette="Dark3")


#Geller's in vitro mut freq
geller<-read.csv("Data/Geller.raw.data.csv")

dt<-geller
dt$TotalReads<-rowSums(dt[,c("Line1_sequence_depth","Line2_sequence_depth","Line3_sequence_depth")], na.rm=T)
dt$a<-rowSums(dt[,c("Line1_substitutions_for_A","Line2_substitutions_for_A","Line3_substitutions_for_A")], na.rm=T)
dt$t<-rowSums(dt[,c("Line1_substitutions_for_U","Line2_substitutions_for_U","Line3_substitutions_for_U")], na.rm=T)
dt$c<-rowSums(dt[,c("Line1_substitutions_for_C","Line2_substitutions_for_C","Line3_substitutions_for_C")], na.rm=T)
dt$g<-rowSums(dt[,c("Line1_substitutions_for_G","Line2_substitutions_for_G","Line3_substitutions_for_G")], na.rm=T)

dt2<-dt[,c(1:5, 43:47)]

#change the sequence letters to lower case
dt2$Sequence<-as.character(sapply(dt2$Sequence, function(v) { return(tolower(v))}))
dt2$Sequence[dt2$Sequence=="u"]<-"t"

for (i in 1:nrow(dt2)){
        MutNum <- dt2[i,paste0(transition(dt2$Sequence[i]))]
        dt2$freq.Ts[i]<-MutNum/dt2$TotalReads[i]
        
        if (dt2$Sequence[i]=="a"|dt2$Sequence[i]=='g'){
                TrvMutNum<-dt2[i,"c"]+dt2[i,"t"]}
        if (dt2$Sequence[i]=="c"|dt2$Sequence[i]=="t"){
                TrvMutNum<-dt2[i,"a"]+dt2[i,"g"]}
        dt2$freq.transv[i]<-TrvMutNum/dt2$TotalReads[i]
} 
        
write.csv(dt2,"Data/Geller.mut.freq.csv")




transmf1<-list()
k<-1
for (i in c("a","t","c","g")) {
        datavector<-dt2$freq.Ts[dt2$Sequence==i]
        nt<-toupper(i)
        dat<-data.frame(NT=rep(nt,times=length(datavector)), mf=datavector)
        transmf1[[k]]<-dat
        names(transmf1)[k]<-nt
        k=k+1
}
mfdata1<-do.call(rbind, transmf1)
mfdata1$Study<-"In vitro"


TS<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F, row.names = 1)

k=1
transmf2<-list()
for (i in  c("a","t","c","g")) {
        datavector<-TS$mean[TS$ref==i]
        nt<-toupper(i)
        dat<-data.frame(NT=rep(nt,times=length(datavector)), mf=datavector)
        transmf2[[k]]<-dat
        names(transmf2)[k]<-nt
        k=k+1
}
mfdata2<-do.call(rbind, transmf2)
mfdata2$Study<-"In vivo"

mfdata<-rbind(mfdata1, mfdata2)
mfdata$Study<-factor(mfdata$Study, level=c("In vivo", "In vitro"))


z=c(0.7,0.3,0.7,0.3)
ggplot(mfdata,aes(x=NT,y=mf,color=Study, fill=Study))+
        scale_y_continuous(trans = 'log10', labels=label_scientific)+
        geom_boxplot(aes(middle=mean(mf), color=Study, fill=Study),outlier.alpha = 0.2)+
        labs(x="",y="Transition mutation frequency")+
        scale_color_manual(values=colors2[c(1,5)]) +
        scale_fill_manual(values=paste0(colors2[c(1,5)],"66" )) +
        theme_bw()+
        theme(axis.text = element_text(size =12, color="black"), axis.title.y = element_text(size=13))+
        theme(legend.title = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
ggsave(filename="Output1A/SummaryFig.Filtered/Invivo.Invitro.comparison.pdf",width=5, height=4, units='in',device='pdf')

nuc<-c("A","T","C","G")
tb<-data.frame(nt=c("A","T","C","G",paste0(nuc,".inVivo")))

for (i in 1:nrow(tb)){
        if (i<=4){
        tb$Mean[i]<-mean(mfdata$mf[mfdata$NT==nuc[i] & mfdata$Study=="In vitro"], na.rm=T)
        tb$SE[i]<-std.error(mfdata$mf[mfdata$NT==nuc[i] & mfdata$Study=="In vitro"], na.rm=T)
        }
        if (i>4){
                
                tb$Mean[i]<-mean(mfdata$mf[mfdata$NT==nuc[i-4] & mfdata$Study=="In vivo"], na.rm=T)
                tb$SE[i]<-std.error(mfdata$mf[mfdata$NT==nuc[i-4] & mfdata$Study=="In vivo"], na.rm=T)
        }
        
}
write.csv(tb, "Output1A/MutFreq.filtered/Invitro.invivo.mf.summary.csv")
tb$NT<-rep(nuc, times=2)
tb$Study<-c(rep("In vitro", times=4), rep("In vivo", times=4))
tb$Study<-factor(tb$Study, level=c("In vivo", "In vitro"))


ggplot()+scale_y_continuous(trans = 'log10', labels=label_scientific)+
        geom_boxplot(data=mfdata, aes(x=NT, y=mf, middle=mean(mf), color=Study, fill=Study),outlier.alpha = 0.2)+
        labs(x="",y="Transition mutation frequency")+
        scale_color_manual(values=colors2[c(1,5)]) +
        scale_fill_manual(values=paste0(colors2[c(1,5)],"66" )) +
        geom_point(data=tb, aes(x=NT, y=Mean, group=Study), position=position_dodge(.75), size=0.6, color="black") +
        theme_bw()+
        theme(axis.text = element_text(size =12, color="black"), axis.title.y = element_text(size=13))+
        theme(legend.title = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())


ggsave(filename="Output1A/SummaryFig.Filtered/Invivo.Invitro.comparison.pdf",width=5, height=4, units='in',device='pdf')






nuc<-c("A","T","C","G")
for (j in 1:4){
        dt<-geller[geller$Sequence==nuc[j],]

        for (i in 1:3){
                dt[,paste0(nuc[j],"->A.freq",i)]<- dt[, paste0("Line",i,"_substitutions_for_A")]/dt[,paste0("Line",i,"_sequence_depth")]
                dt[,paste0(nuc[j],"->G.freq",i)]<- dt[, paste0("Line",i,"_substitutions_for_G")]/dt[,paste0("Line",i,"_sequence_depth")]
                dt[,paste0(nuc[j],"->T.freq",i)]<- dt[, paste0("Line",i,"_substitutions_for_T")]/dt[,paste0("Line",i,"_sequence_depth")]
                dt[,paste0(nuc[j],"->C.freq",i)]<- dt[, paste0("Line",i,"_substitutions_for_C")]/dt[,paste0("Line",i,"_sequence_depth")]
                
        }
        
}
