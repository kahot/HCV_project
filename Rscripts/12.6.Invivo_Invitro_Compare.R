library(tidyverse)
library(zoo)
library(ggplot2)
library(reshape)
library(ggpubr)
library(ggthemes)
library(sfsmisc)
library(colorspace)
source("Rscripts/label_scientific.R")
source("Rscripts/baseRscript.R")

TS<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F, row.names = 1)

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

####
dt2<-read.csv("Data/Geller.mut.freq.csv", stringsAsFactors = F, row.names = 1)


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
        scale_color_manual(values=colors2[c(1,4)]) +
        scale_fill_manual(values=paste0(colors2[c(1,4)],"66" )) +
        theme_bw()+
        theme(axis.text = element_text(size =12, color="black"), axis.title.y = element_text(size=13))+
        theme(legend.title = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
        geom_vline(xintercept = c(1:3)+0.5, color="gray60")
ggsave(filename="Output1A/SummaryFig.Filtered/Invivo.Invitro.comparison.pdf",width=5, height=4, units='in',device='pdf')

####
k=1
mf1<-list()
for (i in  c("a","t","c","g")) {
        for (type in c("syn","nonsyn")){
                datavector<-TS$mean[TS$ref==i & TS$Type==type]
                nt<-toupper(i)
                dat<-data.frame(NT=rep(nt,times=length(datavector)), 
                                type=rep(type, times=length(datavector)),mf=datavector)
                mf1[[k]]<-dat
                names(mf1)[k]<-paste0(nt,".",type)
                k=k+1
        }
}
mfdat1<-do.call(rbind, mf1)

z=rep(c(0.7,0.3),times=4)
x<-1:4
ybreaks<- c(1:10 * 10^c(-4),1:10 * 10^c(-3),1:10 * 10^c(-2),1:10 * 10^c(-1)) 

ggplot(mfdat1,aes(x=NT, y=mf, fill=factor(type)))+geom_boxplot(outlier.alpha =.4, outlier.color = "gray60")+
        scale_y_continuous(trans = 'log10',breaks=c(0.0001,0.001,0.01), minor_breaks=ybreaks, labels=label_scientific2)+
        labs(x="Nucleotide",y="Mutation frequency")+
        scale_fill_manual(values=colors2[c(5,1)]) + theme_bw()+
        theme(legend.title = element_blank()) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10), panel.grid.major.x=element_blank())+
        geom_vline(xintercept = c(1:3)+0.5, color="gray60")
ggsave("Output1A/MutFreq.filtered/MF.byNT.pdf", width = 5,height = 4)





##
#Add among host mf 
among<-read.csv("Output1A/HCV1A_Combined_mutfreq.csv", row.names = 1, stringsAsFactors = F)

k=1
transmf3<-list()
for (i in  c("a","t","c","g")) {
        datavector<-among$freq.Ts[among$ref==i]
        nt<-toupper(i)
        dat<-data.frame(NT=rep(nt,times=length(datavector)), mf=datavector)
        transmf3[[k]]<-dat
        names(transmf3)[k]<-nt
        k=k+1
}
mfdata3<-do.call(rbind, transmf3)
mfdata3$Study<-"Among hosts"

mfdata<-rbind(mfdata1, mfdata2)
mfdata<-rbind(mfdata, mfdata3)
mfdata$Study<-factor(mfdata$Study, level=c( "In vitro","In vivo", "Among hosts"))







###


nuc<-c("A","T","C","G")
tb<-data.frame(nt=c(paste0(nuc, ".inVitro"),paste0(nuc,".inVivo"), paste0(nuc,".AmongHosts")))

for (i in 1:nrow(tb)){
        if (i<=4){
        tb$Mean[i]<-mean(mfdata$mf[mfdata$NT==nuc[i] & mfdata$Study=="In vitro"], na.rm=T)
        tb$SE[i]<-std.error(mfdata$mf[mfdata$NT==nuc[i] & mfdata$Study=="In vitro"], na.rm=T)
        }
        if (i>4& i<=8){
                
                tb$Mean[i]<-mean(mfdata$mf[mfdata$NT==nuc[i-4] & mfdata$Study=="In vivo"], na.rm=T)
                tb$SE[i]<-std.error(mfdata$mf[mfdata$NT==nuc[i-4] & mfdata$Study=="In vivo"], na.rm=T)
        }
        if (i>8){
                
                tb$Mean[i]<-mean(mfdata$mf[mfdata$NT==nuc[i-8] & mfdata$Study=="Among hosts"], na.rm=T)
                tb$SE[i]<-std.error(mfdata$mf[mfdata$NT==nuc[i-8] & mfdata$Study=="Among hosts"], na.rm=T)
        }
}
write.csv(tb, "Output1A/MutFreq.filtered/Invitro.invivo.amonghosts.mf.summary.csv")
tb$NT<-rep(nuc, times=3)
tb$Study<-c(rep("In vitro", times=4), rep("In vivo", times=4),rep("Among hosts", times=4))
tb$Study<-factor(tb$Study, level=c("In vitro","In vivo", "Among hosts"))


ggplot()+scale_y_continuous(trans = 'log10', labels=label_scientific)+
        geom_boxplot(data=mfdata, aes(x=NT, y=mf, middle=mean(mf), color=Study, fill=Study),outlier.alpha = 0.2)+
        labs(x="",y="Transition mutation frequency")+
        scale_color_manual(values=colors2[c(1,3,5)]) +
        scale_fill_manual(values=paste0(colors2[c(1,3,5)],"66" )) +
        geom_point(data=tb, aes(x=NT, y=Mean, group=Study), position=position_dodge(.75), size=0.6, color="black") +
        geom_errorbar(data=tb, aes(x=NT,ymin=Mean-SE, ymax=Mean+SE, group=Study), width=.2, size=.2, position=position_dodge(width=0.75))+
        theme_bw()+
        theme(axis.text = element_text(size =12, color="black"), axis.title.y = element_text(size=13))+
        theme(legend.title = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
        guides(shape = guide_legend(override.aes = list(size = 7)))+
        geom_vline(xintercept = c(1:3)+0.5,  
                   color = "gray60", size=.4)+
ggsave(filename="Output1A/SummaryFig.Filtered/Invivo.Invitro.Amonghosts.comparison.pdf",width=6, height=4, units='in',device='pdf')

## in vivo & in vitro only 
mfdata2<-mfdata[mfdata$Study=="In vitro"|mfdata$Study=="In vivo",]
mfdata2$Study<-factor(mfdata2$Study, levels=c("In vivo", "In vitro"))
tb2<-tb[1:8,]
tb2$Study<-factor(tb2$Study, levels=c("In vivo", "In vitro"))

ggplot()+scale_y_continuous(trans = 'log10', labels=label_scientific)+
        geom_boxplot(data=mfdata2, aes(x=NT, y=mf, middle=mean(mf), color=Study, fill=Study),outlier.alpha = 0.2)+
        labs(x="",y="Transition mutation frequency")+
        scale_color_manual(values=colors2[c(1,4)]) +
        scale_fill_manual(values=paste0(colors2[c(1,4)],"66" )) +
        geom_point(data=tb2, aes(x=NT, y=Mean, group=Study), position=position_dodge(.75), size=0.6, color="black") +
        geom_errorbar(data=tb2, aes(x=NT,ymin=Mean-SE, ymax=Mean+SE, group=Study), width=.2, size=.2, position=position_dodge(width=0.75))+
        theme_bw()+
        theme(axis.text = element_text(size =12, color="black"), axis.title.y = element_text(size=13))+
        theme(legend.title = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
        guides(shape = guide_legend(override.aes = list(size = 7)))+
        geom_vline(xintercept = c(1:3)+0.5,  
                   color = "gray60", size=.4)+
        scale_x_discrete(breaks=c("A","T","C","G"),labels=c(expression(A%->%G),expression("T"%->%C),expression(C%->%"T"),expression(G%->%A)))
ggsave(filename="Output1A/SummaryFig.Filtered/Invivo.Invitro.MFcomparison.pdf",width=5, height=4, units='in',device='pdf')




################
####
#in vivo vs. in vitro by gene
dt2<-read.csv("Data/Geller.mut.freq.csv", stringsAsFactors = F, row.names = 1)
colnames(dt2)[1]<-"pos"

genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
genes$Gene[6]<-"NS1"
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector
end<-TS$pos[nrow(TS)]
genetable<-genetable[genetable$pos>=342&genetable$pos<=end,]

TS<-merge(TS, genetable, by="pos")
TS2<-merge(dt2, genetable, by="pos")

TS<-TS[TS$pos>=TS2$pos[1],]
TS2<-TS2[TS2$pos<=TS$pos[nrow(TS)],]

summary1.mean<-aggregate(TS$mean,by=list(TS$gene),FUN=mean)
summary1.se<-aggregate(TS$mean,by=list(TS$gene),FUN=std.error)
summary2.mean<-aggregate(TS2$freq.Ts,by=list(TS2$gene),FUN=mean, na.rm=T)
summary2.se<-aggregate(TS2$freq.Ts,by=list(TS2$gene),FUN=std.error)

summary<-merge(summary1.mean,summary1.se, by="Group.1")
summary$Study<-"In vivo"

summary2<-merge(summary2.mean,summary2.se, by="Group.1")
summary2$Study<-"In vitro"

summary<-rbind(summary,summary2)
colnames(summary)<-c("Gene", "Mean", "SE", "Study")
summary$Gene<-factor(summary$Gene, levels=c(genes$Gene[2:12]))
summary$Study<-factor(summary$Study, levels=c("In vivo", "In vitro"))

ggplot(summary,aes(x=Gene,y=Mean,color=Study, fill=Study))+
        scale_y_continuous(trans = 'log10', labels=label_scientific)+
        geom_point(data=summary, position=position_dodge(width=0.3), size=2,shape = 21)+scale_color_manual(values=colors2[c(1,4)])+
        scale_fill_manual(values=paste0(colors2[c(1,5)],"66"))+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, position=position_dodge(width=0.3))+
        theme_bw()+theme(axis.text=element_text(size=11), axis.title=element_text(size=13))+ylab("Average mutation frequency")+
        theme(panel.grid.major.x=element_blank(),axis.title.y = element_text(size=13))+
        geom_vline(xintercept = c(1:10)+0.5, color = "gray70", size=.5)+
        theme(axis.title.x=element_blank())
ggsave(filename="Output1A/SummaryFig.Filtered/Invivo.Invitro.byGene.pdf",width=7, height=4)


## Boxplot
transmf1<-list()
k<-1
for (i in 1:11) {
        datavector<-TS2$freq.Ts[TS2$gene==genes$Gene[(i+1)]]
        dat<-data.frame(Gene=rep(genes$Gene[(i+1)],times=length(datavector)), MF=datavector)
        transmf1[[k]]<-dat
        names(transmf1)[k]<-genes$Gene[(i+1)]
        k=k+1
}
mfdata1<-do.call(rbind, transmf1)
mfdata1$Study<-"In vitro"


k=1
transmf2<-list()
for (i in  1:11) {
        datavector<-TS$mean[TS$gene==genes$Gene[(i+1)]]
        dat<-data.frame(Gene=rep(genes$Gene[(i+1)],times=length(datavector)), MF=datavector)
        transmf2[[k]]<-dat
        names(transmf2)[k]<-genes$Gene[(i+1)]
        k=k+1
}
mfdata2<-do.call(rbind, transmf2)
mfdata2$Study<-"In vivo"

mfdata<-rbind(mfdata1, mfdata2)
mfdata$Study<-factor(mfdata$Study, level=c("In vivo", "In vitro"))
mfdata$Gene<-factor(summary$Gene, levels=c(genes$Gene[2:12]))

z=c(0.7,0.3,0.7,0.3)
ggplot(mfdata,aes(x=Gene,y=MF,color=Study, fill=Study))+
        scale_y_continuous(trans = 'log10', labels=label_scientific)+
        geom_boxplot(aes(middle=mean(MF), color=Study, fill=Study),outlier.alpha = 0.2)+
        labs(x="",y="Transition mutation frequency")+
        scale_color_manual(values=colors2[c(1,5)]) +
        scale_fill_manual(values=paste0(colors2[c(1,5)],"66" )) +
        theme_bw()+
        theme(axis.text = element_text(size =12, color="black"), axis.title.y = element_text(size=13))+
        theme(legend.title = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())


ggsave(filename="Output1A/SummaryFig.Filtered/Invivo.Invitro.comparison.pdf",width=5, height=4, units='in',device='pdf')


#### Average number of samples per site in the filtered dataset used to calculate the mean mut freqs
TS<-TS[TS$pos>=342,]
TS$Total<-apply(TS[,2:196], 1, function(x) sum(!is.na(x)) )
mean(TS$Total, na.rm=T) #141.6

6200*mean(TS$Total, na.rm=T) #877903
