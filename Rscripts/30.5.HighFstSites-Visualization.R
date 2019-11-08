library(stringr)
#library(tidyverse)
#library(zoo)
#library(reshape2)
library(ggplot2)
#library(DescTools)
library(colorspace)
source("Rscripts/baseRscript.R")

#dir.create("Output_all/Fst/")
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

#Summary<-read.csv("Output_all/Unfiltered/Ts.Same_Mean_3genotypes.csv",stringsAsFactors = F, row.names = 1)
#Summary$gene[Summary$gene=="NS1(P7)"]<-"NS1"

#n<-data.frame(merge.pos=264:8650)
# Plot the locations of highFst sites
for (g in 1:3){  
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
        
        
        pdf(paste0("Output_all/Fst_T/Diff.site.Fst.byGenes.",fname,  ".pdf"), width=10, height = 30)
        par(mfrow = c(11,1))
        for (i in 1:11){
                tb1<-dt.top[dt.top$gene==genes$Gene[i+1]&dt.top$diff<0,]
                tb1$plot<-0
                tb2<-dt.top[dt.top$gene==genes$Gene[i+1]&dt.top$diff>0,]
                #tb2<-tb2[!is.na(tb2$diff),]
                tb2$plot<-0
                
                dt1<-dt[dt$gene==genes$Gene[i+1],]
                
                #convert nucleotide to number to plot
                #dt1$ref1<-as.numeric(factor(dt1[, paste0("ref.",g1)], levels=c("a","t","c","g")))
            
                plot(dt1$Fst~dt1$merged.pos, pch=16, col="gray70", ylim=c(-0.01,0.08),  main =genes$Gene[i+1],cex=0.7, xlab="Genome position", ylab="Fst")
                points(dt1$diff.site~dt1$merged.pos, pch="|", col=col3)
                points(tb1$Fst~tb1$merged.pos,ylim=c(0,0.08), pch=16, col=col1, cex=0.8)
                points(tb2$Fst~tb2$merged.pos,ylim=c(0,0.08), pch=16, col=col2,cex=0.8)
                points(tb1$plot~tb1$merged.pos,,ylim=c(0,0.08), pch="|", col=col1, cex=0.8 )
                points(tb2$plot~tb2$merged.pos,,ylim=c(0,0.08), pch="|", col=col2, cex=0.8 )
                
                legend("topright", col=c(col1,col2), legend=c(g1,g2), pch=16, cex=.8, bty="n")
                        }
        dev.off()
        
        
        pdf(paste0("Output_all/Fst_T/Diff.site.MF.byGenes.",fname,  ".pdf"), width=10, height = 45)
        par(mfrow = c(22,1))
        for (i in 1:11){
                tb1<-dt.top[dt.top$gene==genes$Gene[i+1]&dt.top$diff<0,]
                tb1$plot<-0
                tb2<-dt.top[dt.top$gene==genes$Gene[i+1]&dt.top$diff>0,]
                tb2$plot<-0
                
                dt1<-dt[dt$gene==genes$Gene[i+1],]
                
                #dt1<-Summary[Summary$gene==genes$Gene[i+1],]
                #dt2<-dt[dt$gene==genes$Gene[i+1],]
                
                
                plot(dt1[,paste0("mean.",g1)]~dt1$merged.pos, pch=16, col="gray70",  main =paste(genes$Gene[i+1], g1),
                     ylim=c(0,0.08), cex=0.6,xlab="Genome position", ylab="Mutation frequency")
                points(tb1[,paste0("mean.",g1)]~tb1$merged.pos,ylim=c(0,0.08), pch=16, col=col1, cex=0.9)
                points(dt1$diff.site~dt1$merged.pos, pch="|", col=col3)
                points(tb1$plot~tb1$merged.pos,ylim=c(0,0.08), pch="|", col=col1, cex=0.8 )
                
                plot(dt1[,paste0("mean.",g2)]~dt1$merged.pos, pch=16, col="gray70",  main =paste(genes$Gene[i+1], g2),
                     ylim=c(0,0.08), cex=0.6,xlab="Genome position", ylab="Mutation frequency")
                points(dt1[,paste0("mean.",g2)]~dt1$merged.pos, pch=16, col="gray70", cex=0.6)
                points(tb2[,paste0("mean.",g2)]~tb2$merged.pos,ylim=c(0,0.08), pch=16, col=col2, cex=0.9)
                points(dt1$diff.site~dt1$merged.pos, pch="|", col=col3)
                points(tb2$plot~tb2$merged.pos,ylim=c(0,0.08), pch="|", col=col2, cex=0.8 )
                
                
        }
        dev.off()
        
        
                
              
        pdf(paste0("Output_all/Fst_T/Diff.site.MF.byGenes2.",fname,  ".pdf"), width=10, height = 45)
        par(mfrow = c(22,1))
        for (i in 1:11){
                tb1<-dt.top[dt.top$gene==genes$Gene[i+1]&dt.top$diff<0,]
                tb1$plot<-0
                tb2<-dt.top[dt.top$gene==genes$Gene[i+1]&dt.top$diff>0,]
                tb2$plot<-0
                
                dt1<-dt[dt$gene==genes$Gene[i+1],]
                
                #dt1<-Summary[Summary$gene==genes$Gene[i+1],]
                #dt2<-dt[dt$gene==genes$Gene[i+1],]
                
                
                plot(dt1[,paste0("mean.",g1)]~dt1$merged.pos, pch=16, col="gray70",  main =paste(genes$Gene[i+1], g1),
                     ylim=c(0,0.08), cex=0.6,xlab="Genome position", ylab="Mutation frequency")
                points(tb1[,paste0("mean.",g1)]~tb1$merged.pos,ylim=c(0,0.08), pch=16, col=col1, cex=0.9)
                points(dt1$same.site~dt1$merged.pos, pch="|", col="gray70",cex=0.7)
                points(dt1$diff.site~dt1$merged.pos, pch="|", col=col3, cex=.7)
                points(tb1$plot~tb1$merged.pos,ylim=c(0,0.08), pch="|", col=col1, cex=0.8 )
                
                plot(dt1[,paste0("mean.",g2)]~dt1$merged.pos, pch=16, col="gray70",  main =paste(genes$Gene[i+1], g2),
                     ylim=c(0,0.08), cex=0.6,xlab="Genome position", ylab="Mutation frequency")
                points(dt1[,paste0("mean.",g2)]~dt1$merged.pos, pch=16, col="gray70", cex=0.6)
                points(tb2[,paste0("mean.",g2)]~tb2$merged.pos,ylim=c(0,0.08), pch=16, col=col2, cex=0.9)
                points(dt1$same.site~dt1$merged.pos, pch="|", col="gray70",cex=0.7)
                points(dt1$diff.site~dt1$merged.pos, pch="|", col=col3)
                points(tb2$plot~tb2$merged.pos,ylim=c(0,0.08), pch="|", col=col2, cex=0.8 )
                
                
        }
        dev.off()
       
                
        pdf(paste0("Output_all/Fst_T/HighFst_allMF.",fname,  "_genes.pdf"), width=8, height = 18)
        par(mfrow = c(6,2))
        for (i in 1:11){
                
                tb1<-dt.top[dt.top$gene==genes$Gene[i+1]&dt.top$diff<0,]
                tb1<-tb1[!is.na(tb1$diff),]
                tb2<-dt.top[dt.top$gene==genes$Gene[i+1]&dt.top$diff>0,]
                tb2<-tb2[!is.na(tb2$diff),]
                
                dt1<-Summary[Summary$gene==genes$Gene[i+1],]
                dt2<-dt[dt$gene==genes$Gene[i+1],]
                nrow(dt1)
                nrow(dt2)
                
                plot(dt1[,paste0("mean.",g1)]~dt1$merged.pos, pch=16, col="gray70",  main =paste(genes$Gene[i+1], g1),
                     ylim=c(0,0.08), cex=0.6,xlab="Genome position", ylab="Mutation frequency")
                points(dt2[,paste0("mean.",g1)]~dt2$merged.pos, pch=16, col="lightskyblue1", cex=0.6)
                points(tb1[,paste0("mean.",g1)]~tb1$merged.pos,ylim=c(0,0.08), pch=16, col=col1, cex=0.9)
                
                plot(dt1[,paste0("mean.",g2)]~dt1$merged.pos, pch=16, col="gray70",  main =paste(genes$Gene[i+1], g2),
                     ylim=c(0,0.08), cex=0.6,xlab="Genome position", ylab="Mutation frequency")
                points(dt2[,paste0("mean.",g2)]~dt2$merged.pos, pch=16, col="lightskyblue1", cex=0.6)
                points(tb2[,paste0("mean.",g2)]~tb2$merged.pos,ylim=c(0,0.08), pch=16, col=col1, cex=0.9)
                
        }
        dev.off() 
        
        
        
        
        
}


