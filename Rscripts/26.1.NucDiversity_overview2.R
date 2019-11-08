#Use overveiw2 data to avoid losing data in the filtered overview (overview3) for Fst calculation

library(ggplot2)
library(reshape)
library(ggpubr)
library(ggthemes)
library(plotrix)
library(grid)
library(purrr)
library(zoo)
library(tidyverse)
library(gplots)
library(colorspace)

source("Rscripts/baseRscript.R")
geno<-c("1A","1B","3A")
colors2<-qualitative_hcl(6, palette="Dark3")
col2_light<-qualitative_hcl(6, palette="Set3")
div.colors<-c(colors2[1],col2_light[1],colors2[3],col2_light[3],colors2[5],col2_light[5])



#If run previously, skip to ## **START** ##

#dir.create("Output3A/Overview_D2/")
#dir.create("Output1B/Overview_D2/")
#dir.create("Output1A/Overview_D2/")

for (f in 1:3){
        flist1<-list.files(paste0("Output",geno[f],"/Overview2/"),pattern="overview2.csv")
        flist2<-list.files(paste0("Output",geno[f],"/SeqData/"),pattern="SeqData")
        
        Overview1<-list()
        for (i in 31:length(flist1)){ 
       
                overview<-read.csv(paste0("Output",geno[f],"/Overview2/",flist1[i]),stringsAsFactors=FALSE, row.names=1)
                #transition info only
                overview<-overview[,c(1:9,20,23,27,31,32,33,36,39)]
                low_reads<-which(overview$TotalReads<1000) 
                overview[low_reads,c(8:ncol(overview))]<-NA
                filename<-substr(paste(flist1[i]),start=1,stop=6)
                
                seq<-read.csv(paste0("Output",geno[f],"/SeqData/",flist2[i]),stringsAsFactors=FALSE, row.names=1)
                seq<-seq[,c(7,1:4)]
                dat<-merge(seq, overview, by="pos", all.y=T )
                dat$D<-""
                for (k in 1:nrow(dat)){
                        if (is.na(dat$freq.Ts[k]))  dat$D[k]<-NA
                        else{
                                ac<-dat[k, "a"]*dat[k,"c"]
                                ag<-dat[k, "a"]*dat[k,"g"]
                                at<-dat[k, "a"]*dat[k,"t"]
                                cg<-dat[k, "c"]*dat[k,"g"]
                                ct<-dat[k, "c"]*dat[k,"t"]
                                gt<-dat[k, "g"]*dat[k,"t"]
                                m<-(dat[k,"TotalReads"]^2-dat[k,"TotalReads"])/2
                                dat$D[k]<-(ac+ag+at+cg+ct+gt)/m
                        }
                }
                
                Overview1[[i]]<-dat
                names(Overview1)[i]<-filename
                print(filename)
                write.csv(dat,paste0("Output",geno[f],"/Overview_D2/",filename,"_overviewD.csv"))
        }
        
        listname<-paste0("Overview1.", geno[f])
        assign(listname, Overview1)
        
        pi<-data.frame(SampleID=names(Overview1))
        for (i in 1:length(Overview1)){
                        df<-Overview1[[i]]
                        n<-nrow(df[!is.na(df$D),])
                        df$D<-as.numeric(df$D)
                        pi$Pi[i]<-(sum(df$D, na.rm=T))/n
                        
        }
        piname<-paste0("Pi.",geno[f])
        assign(piname, pi)
        write.csv(pi, paste0("Output_all/Diversity/NucDiversity.per.sample2.",geno[f],".csv"))
        
        pdf(paste0("Output_all/Diversity/Pi2__",geno[f],".pdf"), height = 5, width=length(flist1)*0.0192+6.1973)
        plot(pi$Pi, xaxt="n", xlab='', ylab='Nucleotide diversity' , pch=16, ylim=c(0,0.03), cex=.8)
        xlabel<-as.character(pi$SampleID)
        xloc<-seq(1,length(flist1), by=1)
        mtext(xlabel, side=1, line=0.5, at=xloc, las=2, cex=0.5)
        dev.off()
}

#Recreate the plots again using saved data

len<-c(195,21,39)
for (f in 1:3) {
        pi<-read.csv(paste0("Output_all/Diversity/NucDiversity.per.sample2.",geno[f],".csv"),stringsAsFactors = F,row.names = 1) 
        
        pdf(paste0("Output_all/Diversity/Pi_overview2_",geno[f],".pdf"), height = 4, width=len[f]*0.0192+6.1973)
        plot(pi$Pi, xaxt="n", xlab='', ylab='Nucleotide diversity' , pch=16, ylim=c(0,0.04), cex=.7)
        xlabel<-as.character(pi$SampleID)
        xloc<-seq(1,len[f], by=1)
        mtext(xlabel, side=1, line=0.5, at=xloc, las=2, cex=0.5)
        dev.off()
}



## **START** ##

# plot nucloetide diversity ave, and se
pi.sum<-data.frame(Genotype=geno)
for (f in 1:3) {
        pi<-read.csv(paste0("Output_all/Diversity/NucDiversity.per.sample2.",geno[f],".csv"),stringsAsFactors = F,row.names = 1) 
        pi.sum$mean[f]<-mean(pi$Pi, na.rm=T)
        pi.sum$se[f]<-std.error(pi$Pi, na.rm=T)
}

ggplot(pi.sum,aes(x=Genotype,y=mean,color=Genotype))+
        geom_point(position=position_dodge(width=0.3),size =2.5)+scale_color_manual(values=colors2[c(1,3,5)])+
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(width=0.3))+
        theme(axis.title.x=element_blank())+ylab("Average nucleotide diversity")+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=13),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:2)+0.5,  
                   color = "gray60", size=.5)

ggsave(filename="Output_all/Diversity/Ave.Pie_by.genotypes.pdf",width = 5, height = 4)


## Calculate pi (nuc diversity) per gene by genotype and plot them
### addd gene annotation info for all genes
Genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
Genes$Gene[6]<-"NS1"
gene.vector<-c()
for (i in 1:12){
        gene.vector<-c(gene.vector, rep(Genes$Gene[i],times=Genes$start[i+1]-Genes$start[i]))
}

n<-data.frame(pos=1:length(gene.vector))
gene.vec<-cbind(n,gene.vector)
colnames(gene.vec)[2]<-"gene"


for(g in 1:3){
        flist<-list.files(paste0("Output",geno[g],"/Overview_D2/"), pattern="overviewD.csv")
        
        Dlist<-list()
        for (i in 1:length(flist)){
                dt<-read.csv(paste0("Output",geno[g],"/Overview_D2/",flist[i]), stringsAsFactors = F, row.names = 1)
                Dlist[[i]]<-dt[,c("pos","D")] 
                fname<-substr(flist[i], start = 1, stop = 6)
                names(Dlist)[i]<-fname
        }
        
        for (i in 1:length(flist)) {
                colnames(Dlist[[i]])<-c("pos",paste0(names(Dlist[i])))
        }
        
        Ds<-Dlist%>% purrr::reduce(full_join, by='pos')
        Ds<-merge(gene.vec, Ds, "pos")
        lname<-paste0("Ds.",geno[g])
        assign(lname, Ds)
        write.csv(Ds, paste0("Output_all/Diversity/Ds_",geno[g],".csv"))

}


SummaryPi<-data.frame()
for (g in 1:3){
        df<-read.csv(paste0("Output_all/Diversity/Ds_",geno[g],".csv"), stringsAsFactors = F, row.names = 1)
        colnames(df)[2]<-"gene"
        lname<-paste0("Ds.",geno[g])
        assign(lname, df)
        
        pi.genes<-data.frame(Gene=Genes$Gene[1:12])
        for (i in 1:12){
               df2<-df[df$gene==Genes$Gene[i],]
               pi.genes[i, "Mean"]<-sum(df2[,3:ncol(df2)], na.rm=T)/ sum(!is.na(df2[,3:ncol(df2)]))
               pi.genes[i,"SE"]<-mean(std.error(df2[,3:ncol(df2)], na.rm=T), na.rm=T) 
        }
        pi.genes$Genotype<-geno[g]
        SummaryPi<-rbind(SummaryPi, pi.genes)
        
}

#write.csv(SummaryPi, "Output_all/Diversity/Pi_byGene_byGenotype.csv")

SummaryPi<-read.csv("Output_all/Diversity/Pi_byGene_byGenotype.csv", row.names = 1)
SummaryPi$Gene<-factor(SummaryPi$Gene, levels=c("5' UTR","Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

#plot the pit calcualted mannually (not by SNPGenie)
ggplot(SummaryPi,aes(x=Gene,y=Mean,group=Genotype, color=Genotype))+
        geom_point(position=position_dodge(width=0.3),size =1.5)+scale_color_manual(values=colors2[c(1,3,5)])+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, position=position_dodge(width=0.3))+
        theme(axis.title.x=element_blank())+ylab(expression(paste("Nucleotide diversity (",pi,")")))+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:11)+0.5,  
                   color = "gray60", size=.4)

ggsave(filename="Output_all/Diversity/MeanPi.byGene.byGenotype.pdf",width = 8.5, height = 5)


####
merged.meta<-read.csv("Output_all/merged.metadata.csv", stringsAsFactors = F, row.names = 1)
meta1A<-merged.meta[,c("org.pos.1A","Type.1A")]
colnames(meta1A)[1]<-"pos"
Pi1A<-merge(meta1A, Ds.1A, by="pos")
colnames(Pi1A)[2:3]<-c("Type","gene")

meta1B<-merged.meta[,c("org.pos.1B","Type.1B")]
colnames(meta1B)[1]<-"pos"
Pi1B<-merge(meta1B, Ds.1B, by="pos")
colnames(Pi1B)[2:3]<-c("Type","gene")

meta3A<-merged.meta[,c("org.pos.3A","Type.3A")]
colnames(meta3A)[1]<-"pos"
Pi3A<-merge(meta3A, Ds.3A, by="pos")
colnames(Pi3A)[2:3]<-c("Type","gene")


#remove 5'UTR
PiNS<-data.frame()
for (g in 1:3){

        df<-get(paste0("Pi",geno[g]))
        dfN<-df[df$gene!="5' UTR" & df$Type=="nonsyn",]
        dfS<-df[df$gene!="5' UTR" & df$Type=="syn",]
        
        piN<-data.frame(Gene=Genes$Gene[2:12])     
        piS<-data.frame(Gene=Genes$Gene[2:12])     

        for (i in 1:11){
                dfN2<-dfN[dfN$gene==Genes$Gene[(i+1)],]
                piN[i, "Mean"]<-sum(dfN2[,4:ncol(dfN2)], na.rm=T)/ sum(!is.na(dfN2[,3:ncol(dfN2)]))
                piN[i,"SE"]<-mean(std.error(dfN2[,3:ncol(dfN2)], na.rm=T), na.rm=T)
                dfS2<-dfS[dfS$gene==Genes$Gene[(i+1)],]
                piS[i, "Mean"]<-sum(dfS2[,4:ncol(dfS2)], na.rm=T)/ sum(!is.na(dfS2[,3:ncol(dfS2)]))
                piS[i,"SE"]<-mean(std.error(dfS2[,3:ncol(dfS2)], na.rm=T), na.rm=T)
        }
       
        piN$Group<-paste0(geno[g],".Nonsyn")
        piN$Type<-"Nonsyn"
        piS$Group<-paste0(geno[g],".Syn")
        piS$Type<-"Syn"
        pi<-rbind(piN,piS)
        pi$Genotype<-geno[g]
        PiNS<-rbind(PiNS, pi)
        
}

PiNS$Gene<-factor(PiNS$Gene, levels=c("Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))
write.csv(PiNS, "Output_all/Diversity/PiNS_Summary.csv")

ggplot(PiNS, aes(x=Gene, y=Mean, group=Group,fill=Group))+
        geom_bar(position=position_dodge(.9), stat="identity",width=0.8,colour="gray60",size=.5)+
        scale_fill_manual(values=div.colors)+
        ylab("nucleotide diversity")+
        geom_errorbar(position=position_dodge(.9),aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, col="gray60")+
        theme_bw()+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13),axis.title.x = element_text(size=13))+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray70", size=.5)+
        theme(panel.grid.major.x = element_blank())

ggsave("Output_all/Diversity/Pi.NS.bygene_bygenotype2.pdf", width = 9, height = 4.5)


## point plot

ggplot(PiNS, aes(x=Gene, y=Mean, shape=Type,color=Genotype))+
        geom_point(position=position_dodge(width=0.8),size =2)+scale_color_manual(values=colors2[c(1,3,5)])+
        scale_shape_manual(values=c(4,19))+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, size=.2, position=position_dodge(width=0.8))+
        theme(axis.title.x=element_blank())+ylab(expression(paste("Nucleotide diversity (",pi,")")))+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray60", size=.4)

ggsave(filename="Output_all/Diversity/PiNS.byGene.byGenotype.pdf",width = 8, height = 5)

#plot PiN/PiS
Pitb<-data.frame()
for (g in 1:3){
        pi.tb<-data.frame(Gene=Genes$Gene[2:12])
        for ( i in 1:11) { 
            pi.tb$pNpS[i]<-PiNS$Mean[PiNS$Genotype==geno[g] & PiNS$Gene==Genes$Gene[i+1] & PiNS$Type=="Nonsyn"]/PiNS$Mean[PiNS$Genotype==geno[g] & PiNS$Gene==Genes$Gene[i+1] & PiNS$Type=="Syn"]
        }
        pi.tb$Genotype<-geno[g]
        Pitb<-rbind(Pitb, pi.tb)
}

Pitb$Gene<-factor(Pitb$Gene, levels=c("Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

ggplot(Pitb, aes(x=Gene, y=pNpS, group=Genotype,fill=Genotype))+
        geom_bar(position=position_dodge(.9), stat="identity",width=0.8,size=.5)+
        scale_fill_manual(values=colors2[c(1,3,5)])+
        ylab(expression(paste(pi[N],"/",pi[S])))+xlab("")+
         theme_bw()+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray70", size=.5)+
        theme(panel.grid.major.x = element_blank())

ggsave("Output_all/Diversity/PiN.over.piS.bygene_bygenotype_bar2.pdf", width = 9, height = 4.5)

#point plot
ggplot(Pitb, aes(x=Gene, y=pNpS, color=Genotype))+
        geom_point(position=position_dodge(width=0.8),size =3)+scale_color_manual(values=colors2[c(1,3,5)])+
        scale_shape_manual(values=c(4,19))+
        theme(axis.title.x=element_blank())+ylab(expression(paste(pi[N],"/",pi[S])))+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  color = "gray60", size=.4)+
        geom_hline(yintercept = 1,  color = "gray30", size=.4)

ggsave(filename="Output_all/Diversity/PiN.over.piS.bygene_bygenotype.point2.pdf",width = 8, height = 5)





# Stats --wilcoxin test 
Comb<-t(combn(1:3, 2))
combnames<-c("1A-1B","1A-3A","1B-3A")
genenames<-as.character(unique(Pitb$Gene))

#data formatting
Pi1A$mean.1A<-rowMeans(Pi1A[4:ncol(Pi1A)], na.rm = T)
meta1a<-merged.meta[,c("merged.pos","org.pos.1A","Type.1A")]
colnames(meta1a)[2]<-"pos"
nuc.div<-merge(meta1a, Pi1A[,c(1:3,199)], by="pos", all.x=T)
nuc.div<-nuc.div[,c("merged.pos", "gene", "Type.1A","mean.1A")]

Pi1B$mean.1B<-rowMeans(Pi1B[4:ncol(Pi1B)], na.rm = T)
meta1b<-merged.meta[,c("merged.pos","org.pos.1B","Type.1B")]
colnames(meta1b)[2]<-"pos"
nuc.div2<-merge(meta1b, Pi1B[,c(1:3,ncol(Pi1B))], by="pos", all.x=T)
nuc.div2<-nuc.div2[,-c(1,4)]

nuc.div<-merge(nuc.div, nuc.div2[,c(1,2,4)], by="merged.pos")

Pi3A$mean.3A<-rowMeans(Pi3A[4:ncol(Pi3A)], na.rm = T)
meta3a<-merged.meta[,c("merged.pos","org.pos.3A","Type.3A")]
colnames(meta3a)[2]<-"pos"
nuc.div2<-merge(meta3a, Pi3A[,c(1:3,ncol(Pi3A))], by="pos", all.x=T)
nuc.div2<-nuc.div2[,-c(1,4)]

nuc.div<-merge(nuc.div, nuc.div2[,c(1,2,4)], by="merged.pos")
nuc.div<-nuc.div[,c("merged.pos","gene","mean.1A","mean.1B","mean.3A","Type.1A","Type.1B","Type.3A")]

#prepare empty data frames
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
        d<-nuc.div[,c("merged.pos", "gene", paste0("mean.",geno[i]),paste0("Type.",geno[i]))]
        colnames(d)[3:4]<-c("mean", "Type")
        
        k=2
        n1<-Comb[i,1]+k
        n2<-Comb[i,2]+k
        
        if (i==1) d1<-nuc.div[nuc.div$Type.1A==nuc.div$Type.1B,]
        if (i==2) d1<-nuc.div[nuc.div$Type.1A==nuc.div$Type.3A,]
        if (i==3) d1<-nuc.div[nuc.div$Type.1B==nuc.div$Type.3A,]
        
        
        for (g in 1:12){
                r1s<-wilcox.test(d1[d1$gene==genenames[g]& d1$Type.1A=="syn",n1], d1[d1$gene==genenames[g] & d1$Type.1A=="syn",n2], alternative = "greater", paired = FALSE) 
                wilcox.resS[g,i]<-r1s[[3]][1]
                r1n<-wilcox.test(d1[d1$gene==genenames[g]& d1$Type.1A=="nonsyn",n1], d1[d1$gene==genenames[g] & d1$Type.1A=="nonsyn",n2], alternative = "greater", paired = FALSE) 
                wilcox.resNS[g,i]<-r1n[[3]][1]
                
                r2<-wilcox.test(d$mean[d$gene==genenames[g]& d$Type=="syn"], d$mean[d$gene==genenames[g] & d$Type=="nonsyn"], alternative = "greater", paired = FALSE) 
                wilcox.res2[g,i]<-r2[[3]][1]
                
                
                
        }
}

write.csv(wilcox.resS,"Output_all/Diversity/WilcoxTest.Pi_byGene.Syn.csv",row.names = T)
write.csv(wilcox.resNS,"Output_all/Diversity/WilcoxTest.Pi_byGene.NS.csv",row.names = T)
write.csv(wilcox.res2,"Output_all/Diversity/WilcoxTest.Pi_NSvsS.csv",row.names = T)




####
# plot sliding window D values across the genome

#merged.meta<-read.csv("Output_all/merged.metadata.csv", stringsAsFactors = F, row.names = 1)

#select sliding window size:  9 codons
ws<-9*3-1

for (g in 1:2){
        df<-get(paste0("Pi",geno[g]))
        df$gene<-as.character(df$gene)
        df$Mean<-rowMeans(df[4:ncol(df)], na.rm=T)
        
        if (g==1|g==2) df<-df[df$pos>=342,];n<-data.frame(pos=342:8641)
        if (g==3) df<-df[df$pos>=340,];n<-data.frame(pos=340:8639)
        df<-df[1:8300,]
        df$SE<-apply(df[4:ncol(df)], 1, function(x) std.error(x,na.rm=T))
        df2<-df[,c("pos","Type","gene","Mean","SE")]
        
        n<-data.frame(pos=342:8641)
        df2<-merge(n,df2, by="pos",all.x=T)
        
        #df2<-df2[1:8250,]
        #rown<-as.integer(nrow(df2)/3)
        codonN<-seq(1,nrow(df2),by=3)
        windDF<-data.frame(window=1:length(codonN))
        for (k in 1:length(codonN)){
                i<-codonN[k]
                df3<-df2[i:(i+ws),]
                windDF$no.N[k]<-table(df3$Type)[1]
                windDF$no.S[k]<-table(df3$Type)[2]
                windDF$pi.N[k]<-mean(df3$Mean[df3$Type=="nonsyn"], na.rm=T)
                windDF$pi.S[k]<-mean(df3$Mean[df3$Type=="syn"], na.rm=T)
        }

        windDF$N_over_S<-windDF$pi.N/windDF$pi.S
        windDF<-windDF[1:2760,]
        
        write.csv(windDF,paste0("Output_all/Diversity/sliding.window.9.", geno[g],".csv"))
        
        ggplot(windDF, aes(x=window, y=N_over_S))+geom_line(size=.2)+
                ylab(expression(paste(pi[N],"/",pi[S])))+
                theme_bw()+
                labs(x="Codon position")+
                theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=12))+
                theme(panel.grid.major.x = element_blank())
        ggsave(paste0("Output_all/Diversity/PiN.PiS.ratio.across.genome.",geno[g],".pdf"), width=12, height=5)
        
        ggplot(windDF, aes(x=window, y=N_over_S))+geom_bar(stat="identity")+
                ylab(expression(paste(pi[N],"/",pi[S])))+
                theme_bw()+
                labs(x="Codon position")+
                theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=12))+
                theme(panel.grid.major.x = element_blank())
        
        
        
}



#######
#Average pi per sample/patient/pop
Pi<-read.csv(paste0("Output_all/Diversity/NucDiversity.per.sample2.1A.csv"), stringsAsFactors = F, row.names = 1)
Pi$Genotype<-"1A"
for(g in 2:3){
        dt<-read.csv(paste0("Output_all/Diversity/NucDiversity.per.sample2.",geno[g],".csv"), stringsAsFactors = F, row.names = 1)
        dt$Genotype<-geno[g]
        Pi<-rbind(Pi, dt)
}









###### Calculate Fst between the samples within genotypes

for (g in 1:1){
        Overview1<-get(paste0("Overview1.",geno[g]))
        ids<-names(Overview1)
        Comb<-t(combn(ids,2))
        FstDF<-data.frame(matrix(ncol=0,nrow=nrow(Comb)))
        pi<-read.csv(paste0("Output_all/Diversity/NucDiversity.per.sample2.",geno[g],".csv"), stringsAsFactors = F, row.names = 1)
        
        for (i in 1:nrow(Comb)){
                samp1<-Comb[i,1]
                samp2<-Comb[i,2]
                df1<-Overview1[[samp1]]
                df2<-Overview1[[samp2]]
                fname<-paste(samp1,samp2)
                FstDF$comb[i]<-fname
                
                df1$DT<-""
                for (k in 1:nrow(df1)){
                        if (is.na(df1$freq.Ts[k])|is.na(df2$freq.Ts[k])){
                                df1$DT[k]<-NA
                        }
                        else{
                                ac<-(df1[k, "a"]+df2[k, "a"])*(df1[k,"c"]+df2[k,"c"])
                                ag<-(df1[k, "a"]+df2[k, "a"])*(df1[k,"g"]+df2[k,"g"])
                                at<-(df1[k, "a"]+df2[k, "a"])*(df1[k,"t"]+df2[k,"t"])
                                cg<-(df1[k, "c"]+df2[k, "c"])*(df1[k,"g"]+df2[k,"g"])
                                ct<-(df1[k, "c"]+df2[k, "c"])*(df1[k,"t"]+df2[k,"t"])
                                gt<-(df1[k, "g"]+df2[k, "g"])*(df1[k,"t"]+df2[k,"t"])
                                m<-((df1[k,"TotalReads"]+df2[k,"TotalReads"])^2-df1[k,"TotalReads"]-df2[k,"TotalReads"])/2
                                df1$DT[k]<-(ac+ag+at+cg+ct+gt)/m
                        }
                }
                df1$DT<-as.numeric(df1$DT)
                n<-nrow(df1[!is.na(df1$DT),])
                FstDF$piT[i]<-sum(df1$DT, na.rm=T)/n
                FstDF$piS.ave[i]<-mean(c(pi$Pi[pi$SampleID==samp1],pi$Pi[pi$SampleID==samp2]),na.rm=T )
                FstDF$Fst[i]<-(FstDF$piT[i]- FstDF$piS.ave[i])/ FstDF$piT[i]
                print(i)
        }
        tname<-paste0("Fst_",geno[g])
        assign(tname,FstDF)
        write.csv(FstDF,paste0("Output_all/Diversity/Fst2_",geno[g],".csv"))
}


#######

save(Overview1.1A, file="Overview1.1A.RData")
save(Overview1.1B, file="Overview1.1B.RData")
save(Overview1.3A, file="Overview1.3A.RData")
load("Overview1.1A.RData")




#############
#
library(ggplot2)
ident<-read.csv("Data/identity.csv",stringsAsFactors = F)
cols3<-c("#66CCEE","#66CCEE","#EE66774D","#EE6677")
colors2<-qualitative_hcl(6, palette="Dark3")


#GenePropm<-melt(GeneProp)
ggplot(ident,aes(x=Gene, y=Pairwise.identity, fill=Genotype))+geom_bar(stat='identity', position = "dodge", width = 0.7)+
        theme_classic()+labs(x="", y="Pairwise identity (%)")+
        #ggtitle(paste0("Over/Under-represented AAs of highly different MF sites (5%) between ",g1," & ", g2))+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = paste0(colors2[c(1,3,5)],"E6"))+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/pairwise-identity.pdf"), width = 7, height = 4.5)


ggplot(ident,aes(x=Gene, y=Identical.Sites, fill=Genotype))+geom_bar(stat='identity', position = "dodge", width = 0.7)+
        theme_classic()+labs(x="", y="% Identical sites")+
        scale_fill_manual("", values = paste0(colors2[c(1,3,5)],"E6"))+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/identical.sites.pdf"), width = 7, height = 4.5)

ggplot(ident,aes(x=Gene, y=GC, fill=Genotype))+geom_bar(stat='identity', position = "dodge", width = 0.7)+
        theme_classic()+labs(x="", y="GC content (%)")+
        scale_fill_manual("", values = paste0(colors2[c(1,3,5)],"E6"))+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/GCcontnet.pdf"), width = 7, height = 4.5)


##
dat<-aggregate(mf[,c(paste0("mean.",geno[g]))], list(mf[,paste0("Type.",geno[g])]), mean, na.rm=T)

#ident.mean<-aggregate(ident[,3:5], list(ident$Genotype) , mean )
ident.sum<-read.csv("Data/identity_summary.csv",stringsAsFactors = F)
ident.sum2<-ident.sum
ident.sum2$Value[ident.sum2$Stat=="Pi"]<-ident.sum$Value[ident.sum$Stat=="Pi"]/10
ident.sum2$Value[ident.sum2$Stat=="ThetaS"]<-ident.sum$Value[ident.sum$Stat=="ThetaS"]/10

ggplot(ident.sum2, aes(x=Genotype, y=Value, group=Stat))+ 
        geom_point(aes(color=Stat))+geom_line(aes(color=Stat))+
        scale_y_continuous(name = "%", sec.axis = sec_axis(~.*10, name = "Theta")) +
        theme_bw()+
        theme(legend.title=element_blank())
ggsave("Output_all/Diversity/Stats_compare.pdf")       



############# calculate pairwise % nucleodiversity by genes: not finished   ###

library(Biostrings)

genes<-read.csv("Data/HCV_annotations2.csv", stringsAsFactors = F)
genenames<-genes$Gene
gene<-c()
for (i in 1:12){
        gene<-c(gene, rep(genenames[i], times=genes$end[i]-genes$start[i]+1))
}

n<-data.frame(pos=1:(length(gene)))
g<-cbind(n,gene)


hcv1a<-read.fasta("Data/HCV1A_Con_Align.fasta", as.string=TRUE)


ids<-names(hcv1a)
Comb<-t(combn(ids,2))
df<-data.frame(matrix(ncol=0,nrow=nrow(Comb)))

PairIdent<-list()
for (i in 1:nrow(Comb)){
        samp1<-Comb[i,1]
        samp2<-Comb[i,2]
        s1<-DNAString(paste0(hcv1a[samp1]))
        s2<-DNAString(paste0(hcv1a[samp2]))
        pair<-pairwiseAlignment(s1,s2)
        PairIdent[[i]]<-pid(pair)
        
        seq1<-paste0(hcv1a[samp1])
        seq2<-paste0(hcv1a[samp2])
        seq11<-unlist(strsplit(seq1, ""))
        seq22<-unlist(strsplit(seq2, ""))
        
        dat<-g
        dat$seq1<-seq11[1:nrow(g)]
        dat$seq2<-seq22[1:nrow(g)]
        for (j in 1:genes)
                seg1<-DNAString(paste0(dat$seq1[dat$gene==genenames[j]]))
                seg2<-DNAString(paste0(dat$seq2[dat$gene==genenames[j]]))
                pair<-pairwiseAlignment(seg1,seg2)
                PairIdent[[j]]<-pid(pair)
                
        
}
s1<-DNAString(paste0(hcv1a[[1]]))
names(hcv1a)
