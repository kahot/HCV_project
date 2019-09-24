library(stringr)
library(tidyverse)
library(zoo)
library(reshape2)
library(ggplot2)
library(DescTools)
library(reshape2)
library(colorspace)
source("Rscripts/baseRscript.R")

#dir.create("Output_all/Fst/")
colors2<-qualitative_hcl(6, palette="Dark3")
col2_light<-qualitative_hcl(6, palette="Set3")
div.colors<-c(colors2[1],col2_light[1],colors2[3],col2_light[3],colors2[5],col2_light[5])
color.genes<-qualitative_hcl(11, palette="Dark3")

#cols3<-c("#009988CC" ,"#66CCEECC", "#EE6677CC", "#4477AACC")

genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
geno<-c("1A","1B","3A")

FstList<-list()
fstfiles<-list.files("Output_all/Fst_T/", pattern=glob2rx("*^Fst*.csv*"))
for (i in 1:3){
        FstList[[i]]<-read.csv(paste0("Output_all/Fst_T/",fstfiles[i]),stringsAsFactors = F, row.names = 1)
        names(FstList)[i]<-substr(fstfiles[i], start=1, stop = 9)
}

Summary<-read.csv("Output_all/Unfiltered/Ts.MinorVariant_Mean_3genotypes.csv",stringsAsFactors = F, row.names = 1)
Summary$gene[Summary$gene=="NS1(P7)"]<-"NS1"

#for (g in 1:3){  
g=1
        dt<-FstList[[g]]
        fname<-substr(names(FstList)[g], start=5, stop = 9)
        g1<-substr(names(FstList)[g], start=5, stop = 6)
        g2<-substr(names(FstList)[g], start=8, stop = 9)
        
        #attached the metadata
        dt<-merge(dt, Summary, by="merged.pos")
        dt<-dt[dt$Type.1A!="stop",]
        dt<-dt[dt$Type.1B!="stop",]
        dt<-dt[dt$Type.3A!="stop",]
        
        #top5%
        s<-as.integer(nrow(dt)*0.05)
        dt$diff<-dt[,paste0("mean.",g1)]-dt[,paste0("mean.",g2)]
        
        
        #select highest 5% Fst
        dt.top<-dt[order(dt$Fst,decreasing = T,na.last = T),]
        dt.top<-dt.top[c(1:s),]
        
        #select sites with MF g1> g2 (costly sites in g2)
        dt1<-dt[dt$diff>0,]
        dt1.top<-dt1[order(dt1$Fst,decreasing = T,na.last = T),]
        dt1.top<-dt1.top[c(1:s),]
        
        #select sites with MF g1 < g2 (costly sites in g1)
        dt2<-dt[dt$diff<0,]
        dt2.top<-dt2[order(dt2$Fst,decreasing = T,na.last = T),]
        dt2.top<-dt2.top[c(1:s),]

        #volcano plot like
        plot(dt$diff, dt$Fst, pch=16, cex=.5, col='grey', xlim=c(-0.08,0.08), xlab=paste0("Difference in mutation frequency (", fname,")"),
             ylab="Fst")
        points(dt.top$diff[dt.top$diff>0], dt.top$Fst[dt.top$diff>0], pch=16, cex=.5,col=colors2[3])
        points(dt.top$diff[dt.top$diff<0], dt.top$Fst[dt.top$diff<0], pch=16, cex=.5,col=colors2[5])
        abline(v=0, col='grey50', lty=3,lwd=1.3)
        
        
        plot(dt$diff, dt$Fst, pch=16, cex=.5, col='grey', xlim=c(-0.08,0.08), xlab=paste0("Difference in mutation frequency (", fname,")"),
             ylab="Fst")
        abline(v=0, col='grey50', lty=3,lwd=1.3)
        abline(h=min(dt.top$Fst), col='grey50', lty=3,lwd=1.3)

        for (i in 1:11){
                points(dt.top$diff[dt.top$gene==genes$Gene[i+1]],dt.top$Fst[dt.top$gene==genes$Gene[i+1]],pch=16, cex=.5,col=color.genes[i])
                
        }
        legend("topleft", col=color.genes, legend=genes$Gene[2:12], pch=16, cex=.5, bty="n")
        
        
        pdf(paste0("Output_all/Fst2/Volvano.plot.",fname,  "_red.pdf"), width=8, height = 18)
        par(mfrow = c(6,2))
        for (i in 1:11){
                plot(dt$diff, dt$Fst, pch=16, cex=.5, col='grey', xlim=c(-0.08,0.08), xlab=paste0("Difference in mutation frequency (", fname,")"),
                     ylab="Fst", main =genes$Gene[i+1])
                abline(v=0, col='grey50', lty=3,lwd=1.3)
                abline(h=min(dt.top$Fst), col='grey50', lty=3,lwd=1.3)
                points(dt.top$diff[dt.top$gene==genes$Gene[i+1]&dt.top$diff>0],dt.top$Fst[dt.top$gene==genes$Gene[i+1]&dt.top$diff>0],pch=16, cex=.5,col="red")
                points(dt.top$diff[dt.top$gene==genes$Gene[i+1]&dt.top$diff<0],dt.top$Fst[dt.top$gene==genes$Gene[i+1]&dt.top$diff<0],pch=16, cex=.5,col="blue")
                n1<-nrow(dt.top[dt.top$gene==genes$Gene[i+1]&dt.top$diff>0,])
                n2<-nrow(dt.top[dt.top$gene==genes$Gene[i+1]&dt.top$diff<0,])
                
                ymax<-max(dt$Fst, na.rm = T)
                text(0.076,(ymax-0.03), n1, cex=0.7, col="red")
                text(-0.076,(ymax-0.03), n2, cex=0.7,col="blue")
                
                if (g==1){
                mtext(text=expression(""%->% "costly in 1B"), side=1, adj=1, cex=.6, line=0.5, col="red")
                mtext(text=expression("costly in 1A" %<-%""), side=1, adj=0, cex=.6,line=0.5, col="blue")
                }
                if (g==2){
                        mtext(text=expression(""%->% "costly in 3A"), side=1, adj=1, cex=.6, line=0.5, col="red")
                        mtext(text=expression("costly in 1A" %<-%""), side=1, adj=0, cex=.6, line=0.5, col="blue")
                }
                if (g==3){
                        mtext(text=expression(""%->% "costly in 3A"), side=1, adj=1, cex=.6, line=0.5, col="red")
                        mtext(text=expression("costly in 1B" %<-%""), side=1, adj=0, cex=.6, line=0.5, col="blue")
                }
        }
        dev.off()

        
        
}
        