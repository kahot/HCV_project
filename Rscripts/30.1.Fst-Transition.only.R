#Calculate between-Patients and between-genotypes Fst using Dvalues
source("Rscripts/baseRscript.R")
library(grid)
library(purrr)
library(zoo)
library(tidyverse)
library(gplots)


geno<-c("1A","1B","3A")

#load the seqdata saved as Overview1
load("Overview1.1A.RData")
load("Overview1.1B.RData")
load("Overview1.3A.RData")

merged.meta<-read.csv("Output_all/merged.metadata.csv", stringsAsFactors = F, row.names = 1)
mergepo<-data.frame(merged.pos=merged.meta[264:8699,"merged.pos"])

for (f in 2:3){
        Overview1<-get(paste0("Overview1.",geno[f]))
        overview1<-Overview1[[1]]
        
        #remove teh sites maj!=ref
        overview1<-overview1[overview1$MajNt==paste0(overview1$ref),]
        overview1<-overview1[!is.na(overview1$pos),]
        #remove the sites with total reads <1000 & 
        for (j in 1:nrow(overview1)){
                if (is.na(overview1$freq.Ts[j]))  {
                        overview1$a[j]<-NA
                        overview1$c[j]<-NA
                        overview1$g[j]<-NA
                        overview1$t[j]<-NA  
                } 
        }

        a<-overview1[,c("pos", "a")]
        c<-overview1[,c("pos", "c")]
        g<-overview1[,c("pos", "g")]
        t<-overview1[,c("pos", "t")]
        
        
        for (i in 2:length(Overview1)){ 
                overview<-Overview1[[i]]
                overview<-overview[overview$MajNt==paste0(overview$ref),]
                overview<-overview[!is.na(overview$pos),]
                #remove the sites with total reads <1000 & 
                for (j in 1:nrow(overview)){
                        if (is.na(overview$freq.Ts[j]))  {
                                overview$a[j]<-NA
                                overview$c[j]<-NA
                                overview$g[j]<-NA
                                overview$t[j]<-NA  } 
                }
                a<-merge(a, overview[,c("pos","a")], by="pos", all = T)
                c<-merge(c, overview[,c("pos","c")], by="pos", all = T)
                g<-merge(g, overview[,c("pos","g")], by="pos", all = T)
                t<-merge(t, overview[,c("pos","t")], by="pos", all = T)
        }
        s<-length(Overview1)
        a$aSum<-rowSums(a[2:(s+1)], na.rm=T)
        c$cSum<-rowSums(c[2:(s+1)], na.rm=T)
        g$gSum<-rowSums(g[2:(s+1)], na.rm=T)
        t$tSum<-rowSums(t[2:(s+1)], na.rm=T)
        
        dat<-merge(a[,c("pos","aSum")],c[,c("pos","cSum")], by="pos", all=T)
        dat<-merge(dat,g[,c("pos","gSum")], by="pos", all=T)
        dat<-merge(dat,t[,c("pos","tSum")], by="pos", all=T)
        colnames(dat)<-c("pos","a","c","g","t")
        dat$Dt<-""
        
        for (k in 1:nrow(dat)){
                if (rowSums(dat[k,2:5])<=3000)  dat$Dt[k]<-NA
                else{
                        ag<-dat[k, "a"]*dat[k,"g"]
                        ct<-dat[k, "c"]*dat[k,"t"]
                        m<-(rowSums(dat[k,2:5])^2-rowSums(dat[k,2:5]))/2
                        dat$Dt[k]<-(ag+ct)/m
                }
        }
        
        dt<-merge(dat, merged.meta, by.x="pos", by.y=paste0("org.pos.",geno[f]), all=T )
        dt<-dt[,-1]
        write.csv(dt, paste0("Output_all/Fst_T/D_",geno[f],".csv"))
}        



 
#data trimming (remove extra rows -run only once)

for (i in 1:3){
        dt<-read.csv(paste0("Output_all/Fst_T/D_", geno[i],".csv"), stringsAsFactors = F, row.names = 1)
        dt<-dt[,c(6,1:5,7:31)]
        dt<-merge(dt, mergepo, by="merged.pos", all.y = T)
        #dt$D[rowSums(dt[,2:5])<3000]<-NA
          
        if (i==1) dt<-dt[dt$merged.pos>263&dt$merged.pos<8645,]
        if (i==2) dt<-dt[dt$merged.pos>263&dt$merged.pos<8602,]
        if (i==3) dt<-dt[dt$merged.pos>263&dt$merged.pos<8636,]
        
        write.csv(dt,paste0("Output_all/Fst_T/D2_",geno[i],".csv") )
}

       
###### Calculate Fst between genotypes
geno<-c("1A","1B","3A")
CombG<-t(combn(geno,2))
FstList<-list()

#normalize the uneven reads between genotypes
for (f in 1:3){
        g1<-CombG[f,1]
        g2<-CombG[f,2]
        fname<-paste0(g1,"-",g2)
        print(fname)
        
        df1<-read.csv(paste0("Output_all/Fst_T/D2_", g1,".csv"), stringsAsFactors = F, row.names = 1)
        df2<-read.csv(paste0("Output_all/Fst_T/D2_", g2,".csv"), stringsAsFactors = F, row.names = 1)
        n<-min(nrow(df1),nrow(df2))
        df1<-df1[1:n,]
        df2<-df2[1:n,]
        
        samesites<-which(df1[,paste0("ref.",g1)]==df1[,paste0("ref.",g2)])
        
        df1<-df1[samesites,]
        df2<-df2[samesites,]
        for (i in 2:5){
                df1[,i]<-as.numeric(df1[,i])
                df2[,i]<-as.numeric(df2[,i])
        }
        
        mergepo<-data.frame(merged.pos=merged.meta[264:n,"merged.pos"])
        df1<-merge(mergepo, df1, by="merged.pos", all.x=T)
        df2<-merge(mergepo, df2, by="merged.pos", all.x=T)
        
        FstDF<-data.frame(merged.pos=df1$merged.pos)
        FstDF$DtT<-""
        
        #divider (m) =(m^2-m)/2   For m=2*m1, divider(m)= 2m1^2-m1
        for (k in 1:nrow(FstDF)){
                if (is.na(df1$Dt[k])|is.na(df2$Dt[k])) FstDF$DtT[k]<-NA
                else{
                        #normalize the reads between the genotypes
                        adj<-sum(df1[k,2:5])/sum(df2[k,2:5])

                        ag<-(df1[k, "a"]+df2[k, "a"]*adj)*(df1[k,"g"]+df2[k,"g"]*adj)
                        ct<-(df1[k, "c"]+df2[k, "c"]*adj)*(df1[k,"t"]+df2[k,"t"]*adj)
                        m<-2*(sum(df1[k,2:5]))^2-sum(df1[k,2:5])
                        FstDF$DtT[k]<-(ag+ct)/m
                }
        }
        
        
        FstDF$DtT<-as.numeric(FstDF$DtT)
        
        for (k in 1:nrow(FstDF)){
                FstDF$Dt.ave[k]<-mean(c(df1$Dt[k],df2$Dt[k]))
                FstDF$Fst[k]<-(FstDF$DtT[k]- FstDF$Dt.ave[k])/ FstDF$DtT[k]
        }
        
        tname<-paste0("Fst.",fname)
        assign(tname,FstDF)
        
        FstList[[f]]<-FstDF
        names(FstList)[f]<-tname
        
        write.csv(FstDF,paste0("Output_all/Fst_T/",tname,".csv"))
}




### Plot Fst

colors2<-qualitative_hcl(6, palette="Dark3")
col2_2<-paste0(colors2,"99")

genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
genenames<-genes$Gene[1:12]

nrow(FstList[[1]]) #6525
nrow(FstList[[2]]) #5649
nrow(FstList[[3]]) #5716




labels<-substr(names(FstList), start = 5,stop=9) 

pdf("Output_all/Fst_T/Fst_across.genome.pdf", width=10, height=4)

plot(Fst~merged.pos, data=FstList[[1]], t="n", xlab='Genome position', ylab="Fst",
     ylim=c(-1,1), xlim=c(264,8400),col=colors2[3], cex=0.2)
points(Fst~merged.pos, data=FstList[[3]],pch=16,col=col2_2[5], ylim=c(-1,1), xlim=c(264,8400),cex=0.2)
points(Fst~merged.pos, data=FstList[[2]],pch=16,col=col2_2[3], ylim=c(-1,1), xlim=c(264,8400),cex=0.2)
points(Fst~merged.pos, data=FstList[[1]],pch=16,col=col2_2[1], ylim=c(-1,1), xlim=c(264,8400),cex=0.2)

#roll100.1<-rollmean(Fst1$Fst, k=100, na.rm=T)
#roll100.2<-rollmean(Fst2$Fst, k=100, na.rm=T)
#roll100.3<-rollmean(Fst3$Fst, k=100, na.rm=T)
#
#roll100.1a<-c(rep(NA, times=(264+100)),roll100.1)
#roll100.2b<-c(rep(NA, times=(264+100)),roll100.2)
#roll100.3a<-c(rep(NA, times=(264+100)),roll100.3)
#
#lines(roll100.2b,  col=colors2[3],lwd=1)
#lines(roll100.3a,  col=colors2[5],lwd=1)
#lines(roll100.1a,  col=colors2[1],lwd=1)


legend("topright", legend=labels, col=colors2[c(1,3,5)], pch=16, bty="n", cex=0.8)


genes2<-genes[1:12,]

ylow<--1.14
for (j in 1:nrow(genes2)){
        xleft<- genes2$start[j]
        xright<-genes2$start[j+1]
        if (j==1){
                
                rect(-200,ylow,xright,ylow*1.6,density = NULL, angle = 45,col="white",border ="gray60",lwd=1.5)
                text(150,1.2*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else if (j==6){
                rect(xleft,ylow,xright,ylow*1.6,density = NULL, angle = 45,col="white",border ="gray60",lwd=1.5)
                text(xleft+80,1.2*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
                
        }
        else if (j==4|j==9){
                rect(xleft,ylow,xright,1.6*ylow,density = NULL, angle = 45,col="white",border ="gray60",lwd=1.5)
                text(xleft+50,0.8*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else if (j==12){
                rect(xleft,ylow,genes2$end[j],1.6*ylow,density = NULL, angle = 45,col="white",border ="gray60",lwd=1.5)
                text(xleft+600,1.2*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else{rect(xleft,ylow,xright,1.6*ylow,density = NULL, angle = 45,col="white",border ="gray60",lwd=1.5)
                text(xright-(xright-xleft)/2,1.2*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)}
}


box()

######
types<-merged.meta[,c("merged.pos", "Type.1A", "Type.1B","Type.3A")]
Comb<-t(combn(1:3, 2))
combnames<-c("1A-1B","1A-3A","1B-3A")

fst.geno<-data.frame(pop=combnames)
for (i in 1:3){
        dat<-FstList[[i]]
        dat<-merge(types, dat, by="merged.pos")
        
        n1<-Comb[i,1]+1
        n2<-Comb[i,2]+1
        
        fst.geno$f[i]<-mean(dat$Fst, na.rm = T)
        fst.geno$fs[i]<-mean(dat$Fst[dat[,n1]=="syn" & dat[,n2]=="syn" ], na.rm = T)
        fst.geno$fn[i]<-mean(dat$Fst[dat[,n1]=="nonsyn" & dat[,n2]=="nonsyn"  ], na.rm = T)
        
}


colorRampBlue <- colorRampPalette(c("white", "steelblue1", "blue3"))

pairpi<-read.csv("Output_all/Diversity/Pairwise.pi.csv", stringsAsFactors = F)

fst<-pairpi[,2:3]
fst$Fst[c(1,5,6)]<-0
fst$Fst[c(2,3,4)]<-fst.geno$fs[1:3]
fst[7,1:2]<-c("1B","1A")
fst[8,1:2]<-c("3A","1A")
fst[9,1:2]<-c("3A","1B")
fst$Fst[c(7:9)]<-fst.geno$fn[1:3]
fst$type<-c(NA,"syn","syn","syn",NA,NA,"nonsyn","nonsyn","nonsyn")

ggplot(data=fst, aes(Sample1,Sample2, fill=Fst)) +geom_tile(color="gray50", size=.3)+ 
        scale_fill_gradientn(colors=colorRampBlue(64), limit=c(0,0.07))+
        theme_bw()+
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
              panel.grid.major = element_blank(),panel.border = element_blank())+
        geom_rect(aes(xmin = 1 - 0.5, xmax = 3 + 0.5, ymin = 1 - 0.5, ymax = 3 + 0.5),
                  fill = "transparent", color = "gray40", size = .3)+
        geom_text(aes(label = round(Fst, 4)))+
        geom_text(aes(label=type), nudge_y=0.2)
ggsave("Output_all/Fst_T/Fst.genotypes.heatmap.pdf",width=5, height=4.3, units='in',device='pdf')





