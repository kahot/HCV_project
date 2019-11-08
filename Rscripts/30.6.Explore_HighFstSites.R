library(stringr)
#library(tidyverse)
#library(zoo)
library(reshape2)
library(ggplot2)
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


#Start with 1A vs. 3A
#for (g in 1:3){  
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
        
        #if (g==1) {col1<-colors2[1]; col2<-colors2[3]; col3<-"blue"}
        #if (g==2) {col1<-colors2[1]; col2<-colors2[5]; col3<-"green"}
        #if (g==3) {col1<-colors2[3]; col2<-colors2[5]; col3<-"red"}
        
        
        dt.top1A<-dt.top[dt.top$diff<0,]
        dt.top1A<-dt.top1A[!is.na(dt.top1A$diff),]
       
        dt.top3A<-dt.top[dt.top$diff>0,]
        dt.top3A<-dt.top3A[!is.na(dt.top3A$diff),]
        
        
     
#extract mutation frequenices at the highest Fst sites and plot them
MF1a<-read.csv("Output_all/Unfiltered/Ts.same_1A.csv", stringsAsFactors = F, row.names = 1)
MF3a<-read.csv("Output_all/Unfiltered/Ts.same_3A.csv", stringsAsFactors = F, row.names = 1)
colnames(MF1a)[2:196]<-1:195

Mvf1a<-read.csv("Output_all/Unfiltered/Ts.MinorVarient_1A.csv", stringsAsFactors = F, row.names = 1)
Mvf3a<-read.csv("Output_all/Unfiltered/Ts.MinorVarient_3A.csv", stringsAsFactors = F, row.names = 1)
meta<-read.csv("Output_all/Unfiltered/Ts.MinorVariant_Mean_metadata.1A.csv", stringsAsFactors = F, row.names = 1)

Pos<-dt.top$merged.pos[dt.top$Fst>=0.02]  #22 positions

for (j in 1:length(Pos)){
        dt1<-data.frame(t(MF1a[MF1a$merged.pos==Pos[j],2:196]))
        dt1$Genotype<-"1A"
        dt1<-dt1[!is.na(dt1[,1]),]
        
        dt2<-data.frame(t(MF3a[MF3a$merged.pos==Pos[j],2:40]))
        dt2$Genotype<-"3A"
        dt2<-dt2[!is.na(dt2[,1]),]
        
        if (nrow(dt1)<195/3 | nrow(dt2)<39/3){
                next
        }
        else {
                mfs<-rbind(dt1, dt2)
                colnames(mfs)[1]<-"MF"
                
                Sum<-data.frame(mean=c(mean(dt1[,1], na.rm=T),mean(dt2[,1], na.rm=T)),se=c(std.error(dt1[,1], na.rm=T),std.error(dt2[,1], na.rm=T)))
                Sum$Genotype<-c("1A","3A")
                
                ggplot()+scale_y_continuous(trans='log10', label=label_scientific)+
                        ggtitle(paste0("Position ",Pos[j]," (", dt.top$Type.1A[dt.top$merged.pos==Pos[j]],")"))+
                        geom_point(stat="identity",data=mfs, aes(x=Genotype, y=MF, color=Genotype),position_jitter(0.05), size=0.5)+
                        ylab("Mutation freqeuncy")+
                        geom_point(data=Sum, aes(x=Genotype, y=mean))+
                        geom_errorbar(data=Sum,aes(x=Genotype, ymin=mean-se, ymax=mean+se), width=.1, size=.2)+
                        scale_color_manual(values=colors2[c(1,5)])+
                        theme_bw()+
                        theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(color="black", size=11))+
                        geom_vline(xintercept = 1.5, color = "gray60", size=.4)+
                        theme(legend.position = "none")
                ggsave(filename=paste0("Output_all/Fst_T/Plots/Position-",Pos[j], ".",fname,".plot.pdf"),width=4, height=5.5)
        
                k<-dt.top$codon[dt.top$merged.pos==Pos[j]]
                if (k==3){ low<-6; high=4}
                if (k==2){ low<-5; high=5}
                if (k==1){ low<-4; high=6}
                
                Sdt1<-data.frame(t(Mvf1a[Mvf1a$merged.pos>(Pos[j]-k-3) & Mvf1a$merged.pos<(Pos[j]+7-k),2:196]))
                colnames(Sdt1)<-c((Pos[j]-k-2):(Pos[j]+6-k))
                Sdt1m<-melt(Sdt1)
                colnames(Sdt1m)<-c("Position", "MF")
                Sdt1m$Genotype<-"1A"
                
                Sdt2<-data.frame(t(Mvf3a[Mvf1a$merged.pos>(Pos[j]-k-3) & Mvf3a$merged.pos<(Pos[j]+7-k),2:40]))
                colnames(Sdt2)<-c((Pos[j]-k-2):(Pos[j]+6-k))
                Sdt2m<-melt(Sdt2)
                colnames(Sdt2m)<-c("Position", "MF")
                Sdt2m$Genotype<-"3A"
                
                DF<-rbind(Sdt1m,Sdt2m)
                
                
                summary<-data.frame(Position=colnames(Sdt1))
                summary$mean<-colMeans(Sdt1, na.rm = T)
                for (i in 1:ncol(Sdt1)){
                        summary$se[i]<-std.error(Sdt1[,i], na.rm=T)
                        
                }
                summary$Genotype<-"1A"
                
                summary2<-data.frame(Position=colnames(Sdt2))
                summary2$mean<-colMeans(Sdt2, na.rm = T)
                for (i in 1:ncol(Sdt2)){
                        summary2$se[i]<-std.error(Sdt2[,i], na.rm=T)
                        
                }
                summary2$Genotype<-"3A"
                summary<-rbind(summary,summary2)
                
                nt1<-meta$ref.1A[meta$merged.pos>=(Pos[j]-k-2)& meta$merged.pos<=(Pos[j]+6-k)]
                nt2<-meta$ref.3A[meta$merged.pos>=(Pos[j]-k-2)& meta$merged.pos<=(Pos[j]+6-k)]
                nt1<-toupper(nt1)
                nt2<-toupper(nt2)
                
                ggplot()+
                        geom_point(data=DF, aes(x=Position, y=MF, color=Genotype),stat="identity",position=position_dodge(width = .6), size=0.5)+
                        scale_y_continuous(trans='log10', label=label_scientific, limits=c(0.00001,0.5))+
                        ylab("Mutation freqeuncy")+
                        scale_color_manual(values=colors2[c(1,5)])+
                        geom_point(data=summary, aes(x=Position, y=mean, group=Genotype),position=position_dodge(width = .6))+
                        geom_errorbar(data=summary,aes(x=Position, ymin=mean-se, ymax=mean+se, group=Genotype), width=.3, size=.2, position=position_dodge(width = .6))+
                        theme_bw()+ggtitle(paste0("Position ",Pos[j]," (", dt.top$Type.1A[dt.top$merged.pos==Pos[j]],")"))+
                        theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(color="black", size=10))+
                        geom_vline(xintercept = c(1:8)+0.5, color = "gray80", size=.2, linetype=2)+
                        geom_vline(xintercept = c(3,6)+0.5, color = "gray40", size=.4)+
                        theme(legend.position = "none")+
                        annotate("Text", x=1:9, y=0.000025,label=nt1, color=colors2[1])+
                        annotate("Text", x=1:9, y=0.00001,label=nt2, color=colors2[5])
                ggsave(filename=paste0("Output_all/Fst_T/Plots/9Bases-",Pos[j], ".",fname,".plot.pdf"),width=6, height=5.5)
                       
               
        }
                
}
        
        

##
Sdt1<-data.frame(t(Mvf1a[Mvf1a$merged.pos>1202 & Mvf1a$merged.pos<1212,2:196]))
colnames(Sdt1)<-1203:1211
Sdt1m<-melt(Sdt1)
colnames(Sdt1m)<-c("Position", "MF")
Sdt1m$Genotype<-"1A"

Sdt2<-data.frame(t(Mvf3a[Mvf3a$merged.pos>1202 & Mvf3a$merged.pos<1212,2:40]))
colnames(Sdt2)<-1203:1211
Sdt2m<-melt(Sdt2)
colnames(Sdt2m)<-c("Position", "MF")
Sdt2m$Genotype<-"3A"

DF<-rbind(Sdt1m,Sdt2m)


summary<-data.frame(Position=colnames(Sdt1))
summary$mean<-colMeans(Sdt1, na.rm = T)
for (i in 1:ncol(Sdt1)){
        summary$se[i]<-std.error(Sdt1[,i], na.rm=T)
        
}
summary$Genotype<-"1A"
        
summary2<-data.frame(Position=colnames(Sdt2))
summary2$mean<-colMeans(Sdt2, na.rm = T)
for (i in 1:ncol(Sdt2)){
        summary2$se[i]<-std.error(Sdt2[,i], na.rm=T)
        
}
summary2$Genotype<-"3A"
summary<-rbind(summary,summary2)


ggplot()+
        geom_point(data=DF, aes(x=Position, y=MF, color=Genotype),stat="identity",position=position_dodge(width = .6), size=0.5)+
        scale_y_continuous(trans='log10', label=label_scientific, limits=c(0.00004,0.5))+
        ylab("Mutation freqeuncy")+
        scale_color_manual(values=colors2[c(1,5)])+
        geom_point(data=summary, aes(x=Position, y=mean, group=Genotype),position=position_dodge(width = .6))+
        geom_errorbar(data=summary,aes(x=Position, ymin=mean-se, ymax=mean+se, group=Genotype), width=.3, size=.2, position=position_dodge(width = .6))+
        theme_bw()+
        theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(color="black", size=10))+
        geom_vline(xintercept = c(1:8)+0.5, color = "gray80", size=.2, linetype=2)+
        geom_vline(xintercept = c(3,6)+0.5, color = "gray40", size=.4)+
        theme(legend.position = "none")+
        annotate("Text", x=1:9, y=0.00004,label=c("G","G","T","C","A","A","C","T","G"))+
        annotate("Text", x=1:9, y=0.00007,label=c("G","G","A","C","A","A","G","C","C"))
        
        annotate("rect",xmin=4.5, xmax=5.5,ymin=0.000036, ymax=0.6, alpha=.2 )



###
# 1545
pos<-1545
dt1<-data.frame(t(MF1a[MF1a$merged.pos==pos,2:196]))
dt1$Genotype<-"1A"
row.names(dt1)<-1:195

dt2<-data.frame(t(MF3a[MF3a$merged.pos==pos,2:40]))
dt2$Genotype<-"3A"

nrow(dt1[!is.na(dt1[,1]),]) #69
nrow(dt2[!is.na(dt2[,1]),]) #18
range(dt1[,1], na.rm=T) # 0.0007745933 0.3876041
range(dt2[,1], na.rm=T) # 0.0000000  0.03113525
mean(dt1 [,1], na.rm=T) # 0.05966333
mean(dt2 [,1], na.rm=T) # 0.003184129

C1545<-rbind(dt1, dt2)
colnames(C1545)[1]<-"mf"

C1545Sum<-data.frame(Average=c(mean(dt1[,1], na.rm=T),mean(dt2[,1], na.rm=T)),se=c(std.error(dt1[,1], na.rm=T),std.error(dt2[,1], na.rm=T)))
C1545Sum$Genotype<-c("1A","3A")


ggplot()+scale_y_continuous(trans='log10', label=label_scientific, expand = expand_scale(mult = c(0.2, .1)))+
        geom_point(stat="identity",data=C1545, aes(x=Genotype, y=mf, color=Genotype),position_jitter(0.05), size=0.5)+
        ylab("Mutation freqeuncy")+
        geom_point   (data=C1545Sum, aes(x=Genotype, y=Average))+
        geom_errorbar(data=C1545Sum,aes(x=Genotype, ymin=Average-se, ymax=Average+se), width=.1, size=.2)+
        scale_color_manual(values=colors2[c(1,5)])+
        theme_bw()+
        theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(color="black", size=11))+
        geom_vline(xintercept = 1.5, color = "gray60", size=.4)+
        theme(legend.position = "none")
ggsave(filename="Output_all/Fst_T/C1545.plot.pdf",width=4, height=5.5)
