library(stringr)
#library(tidyverse)
#library(zoo)
library(reshape2)
library(ggplot2)
library(dplyr)
#library(DescTools)
library(colorspace)
source("Rscripts/baseRscript.R")
source("Rscripts/label_scientific.R")


colors2<-qualitative_hcl(6, palette="Dark3")
col2_light<-qualitative_hcl(6, palette="Set3")
div.colors<-c(colors2[1],col2_light[1],colors2[3],col2_light[3],colors2[5],col2_light[5])
color.genes<-qualitative_hcl(11, palette="Dark3")

genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
geno<-c("1A","1B","3A")

######
topFst<-list()
fsts<-list.files("Output_all/Fst_T/", pattern=glob2rx("*^HighFst*.csv*"))
for (i in 1:3){
        topFst[[i]]<-read.csv(paste0("Output_all/Fst_T/",fsts[i]),stringsAsFactors = F, row.names = 1)
        names(topFst)[i]<-substr(fsts[i], start=14, stop = 18)
}
summaryFst<-list()
fsts2<-list.files("Output_all/Fst_T/", pattern=glob2rx("*^Summary.Fst*.csv*"))
for (i in 1:3){
        summaryFst[[i]]<-read.csv(paste0("Output_all/Fst_T/",fsts2[i]),stringsAsFactors = F)
        names(summaryFst)[i]<-substr(fsts[i], start=14, stop = 18)
}

Summary<-read.csv("Output_all/Unfiltered/Ts.MinorVariant_Mean_3genotypes.csv",stringsAsFactors = F, row.names = 1)
Summary$gene[Summary$gene=="NS1(P7)"]<-"NS1"
Summary<-Summary[Summary$merged.pos>263,]


#extract mutation frequenices at the highest Fst sites and plot them
MF1A<-read.csv("Output_all/Unfiltered/Ts.same_1A.csv", stringsAsFactors = F, row.names = 1)
MF1B<-read.csv("Output_all/Unfiltered/Ts.same_1B.csv", stringsAsFactors = F, row.names = 1)
MF3A<-read.csv("Output_all/Unfiltered/Ts.same_3A.csv", stringsAsFactors = F, row.names = 1)
colnames(MF1)[2:196]<-1:195

Mvf1A<-read.csv("Output_all/Unfiltered/Ts.MinorVarient_1A.csv", stringsAsFactors = F, row.names = 1)
Mvf1B<-read.csv("Output_all/Unfiltered/Ts.MinorVarient_1B.csv", stringsAsFactors = F, row.names = 1)
Mvf3A<-read.csv("Output_all/Unfiltered/Ts.MinorVarient_3A.csv", stringsAsFactors = F, row.names = 1)
#meta<-read.csv("Output_all/Unfiltered/Ts.MinorVariant_Mean_metadata.1A.csv", stringsAsFactors = F, row.names = 1)




for (g in 1:3){  
g=2
        dt<-summaryFst[[g]]
        dt<-dt[,c("merged.pos", "Fst")]
        dt.top<-topFst[[g]]
        fname<-names(summaryFst)[g]
        g1<-substr(fname, start=1, stop = 2)
        g2<-substr(fname, start=4, stop = 5)
        
        dt<-merge(dt, Summary, by="merged.pos", all.y = T)
        dt$diff.site<-NA
        dt$diff.site[dt[paste0("ref.", g1)]!=dt[,paste0("ref.", g2)]]<-0
        dt$same.site<-NA
        dt$same.site[dt[paste0("ref.", g1)]==dt[,paste0("ref.", g2)]]<-0
        
        if (g==1) {col1<-colors2[1]; col2<-colors2[3]; col3<-"blue"}
        if (g==2) {col1<-colors2[1]; col2<-colors2[5]; col3<-"green"}
        if (g==3) {col1<-colors2[3]; col2<-colors2[5]; col3<-"red"}
        
        
        dt.top1<-dt.top[dt.top$diff<0,] #costly in g1
        dt.top1<-dt.top1[!is.na(dt.top1$diff),]
       
        dt.top2<-dt.top[dt.top$diff>0,]  #costly in g2
        dt.top2<-dt.top2[!is.na(dt.top2$diff),]
 
        Pos<-dt.top$merged.pos[dt.top$Fst>=0.01]  #29 positions
        
        #Pos<-dt.top$merged.pos[dt.top$Fst<0.02&dt.top$Fst>=0.01]
        
        MF1<-get(paste0("MF", g1))
        MF2<-get(paste0("MF", g2))
        Mvf1<-get(paste0("Mvf", g1))
        Mvf2<-get(paste0("Mvf", g2))
        
        s1<- ncol(MF1)-2
        s2<- ncol(MF2)-2
        HiFstPos<-c()
        for (j in 1:length(Pos)){
                dt1<-data.frame(t(MF1[MF1$merged.pos==Pos[j],2:(s1+1)]))
                dt1$Genotype<-g1
                dt1<-dt1[!is.na(dt1[,1]),]
                
                dt2<-data.frame(t(MF2[MF2$merged.pos==Pos[j],2:(s2+1)]))
                dt2$Genotype<-g2
                dt2<-dt2[!is.na(dt2[,1]),]
                
                if (nrow(dt1)<s1/3 | nrow(dt2)<s2/3){
                        next
                }
                else {
                        mfs<-rbind(dt1, dt2)
                        colnames(mfs)[1]<-"MF"
                        
                        Sum<-data.frame(mean=c(mean(dt1[,1], na.rm=T),mean(dt2[,1], na.rm=T)),se=c(std.error(dt1[,1], na.rm=T),std.error(dt2[,1], na.rm=T)))
                        Sum$Genotype<-c(g1,g2)
                        
                        ggplot()+scale_y_continuous(trans='log10', label=label_scientific)+
                                ggtitle(paste0("Position ",Pos[j]," (", dt.top$Type.1A[dt.top$merged.pos==Pos[j]],")"))+
                                geom_point(stat="identity",data=mfs, aes(x=Genotype, y=MF, color=Genotype),position_jitter(0.05), size=0.5)+
                                ylab("Mutation freqeuncy")+
                                geom_point(data=Sum, aes(x=Genotype, y=mean))+
                                geom_errorbar(data=Sum,aes(x=Genotype, ymin=mean-se, ymax=mean+se), width=.1, size=.2)+
                                scale_color_manual(values=c(col1, col2))+
                                theme_bw()+
                                theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(color="black", size=11))+
                                geom_vline(xintercept = 1.5, color = "gray60", size=.4)+
                                theme(legend.position = "none")
                        ggsave(filename=paste0("Output_all/Fst_T/Plots/",fname, ".Position-",Pos[j], ".plot.pdf"),width=4, height=5.5)
                
                        k<-dt.top$codon[dt.top$merged.pos==Pos[j]]
                        
                        
                        Sdt1<-data.frame(t(Mvf1[Mvf1$merged.pos>(Pos[j]-k-3) & Mvf1$merged.pos<(Pos[j]+7-k),2:196]))
                        colnames(Sdt1)<-c((Pos[j]-k-2):(Pos[j]+6-k))
                        Sdt1m<-melt(Sdt1)
                        colnames(Sdt1m)<-c("Position", "MF")
                        Sdt1m$Genotype<-g1
                        
                        Sdt2<-data.frame(t(Mvf2[Mvf2$merged.pos>(Pos[j]-k-3) & Mvf2$merged.pos<(Pos[j]+7-k),2:40]))
                        colnames(Sdt2)<-c((Pos[j]-k-2):(Pos[j]+6-k))
                        Sdt2m<-melt(Sdt2)
                        colnames(Sdt2m)<-c("Position", "MF")
                        Sdt2m$Genotype<-g2
                        
                        DF<-rbind(Sdt1m,Sdt2m)
                        
                        
                        summary<-data.frame(Position=colnames(Sdt1))
                        summary$mean<-colMeans(Sdt1, na.rm = T)
                        for (i in 1:ncol(Sdt1)){
                                summary$se[i]<-std.error(Sdt1[,i], na.rm=T)
                                
                        }
                        summary$Genotype<-g1
                        
                        summary2<-data.frame(Position=colnames(Sdt2))
                        summary2$mean<-colMeans(Sdt2, na.rm = T)
                        for (i in 1:ncol(Sdt2)){
                                summary2$se[i]<-std.error(Sdt2[,i], na.rm=T)
                                
                        }
                        summary2$Genotype<-g2
                        summary<-rbind(summary,summary2)
                        
                        c1<-which(colnames(meta)==paste0("ref.",g1))
                        c2<-which(colnames(meta)==paste0("ref.",g2))
                        
                        
                        nt1<-meta[meta$merged.pos>=(Pos[j]-k-2)& meta$merged.pos<=(Pos[j]+6-k),c1]
                        nt2<-meta[meta$merged.pos>=(Pos[j]-k-2)& meta$merged.pos<=(Pos[j]+6-k),c2]
                        nt1<-toupper(nt1)
                        nt2<-toupper(nt2)
                        rec1<-data.frame(xmin=(k+2.5), xmax=(k+3.5), ymin=-Inf, ymax=Inf)
                        ggplot()+
                                geom_point(data=DF, aes(x=Position, y=MF, color=Genotype),stat="identity",position=position_dodge(width = .6), size=0.5)+
                                scale_y_continuous(trans='log10', label=label_scientific, limits=c(0.00001,0.5))+
                                ylab("Mutation freqeuncy")+
                                scale_color_manual(values=c(col1, col2))+
                                geom_point(data=summary, aes(x=Position, y=mean, group=Genotype),position=position_dodge(width = .6))+
                                geom_errorbar(data=summary,aes(x=Position, ymin=mean-se, ymax=mean+se, group=Genotype), width=.3, size=.2, position=position_dodge(width = .6))+
                                theme_bw()+ggtitle(paste0("Position ",Pos[j]," (", dt.top$Type.1A[dt.top$merged.pos==Pos[j]],")"))+
                                theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(color="black", size=10))+
                                geom_vline(xintercept = c(1:8)+0.5, color = "gray80", size=.2, linetype=2)+
                                geom_vline(xintercept = c(3,6)+0.5, color = "gray40", size=.4)+
                                theme(legend.position = "none")+
                                annotate("Text", x=1:9, y=0.000025,label=nt1, color=col1)+
                                annotate("Text", x=1:9,Fst y=0.00001,label=nt2, color=col2)+
                                geom_rect(data=rec1, aes(xmin=xmin, xmax=xmax, ymin=0, ymax=ymax), fill="gray", alpha=0.2, stat="identity") 
                                
                        ggsave(filename=paste0("Output_all/Fst_T/Plots/",fname,"/", fname,".9Bases-",Pos[j], ".plot.pdf"),width=6, height=5.5)
                               
                       HiFstPos<-c(HiFstPos,Pos[j])
                }
                vname<-paste0("HiFstPos.", fname)
                assign(vname, HiFstPos)        
        }
        
        
## which ones are CpG Creating?
HiFstPos1<-`HiFstPos.1A-3A` #51 sites
top.1A3A<-dt.top[dt.top$merged.pos %in% HiFstPos1,]
top.1A3A<-top.1A3A[top.1A3A$Type.1A==top.1A3A$Type.3A,]

write.csv(top.1A3A, "Output_all/Fst_T/Plots/1A-3A/HighFstSites.summary.csv")

table(top.1A3A$ref.1A)
#a  c  g  t 
#8 21 13  8

table(top.1A3A$codon)
# 1  2  3 
#11  7 32 
table(top.1A3A$codon, top.1A3A$Type.1A)
