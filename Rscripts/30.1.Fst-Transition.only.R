#Calculate between-Patients and between-genotypes Fst using Dvalues
source("Rscripts/baseRscript.R")
library(grid)
library(purrr)
library(zoo)
library(tidyverse)
library(gplots)
library(colorspace)

geno<-c("1A","1B","3A")

#create Overview2.2 list for each genotype
for (g in 2:3){
        HCVFiles<-list.files(paste0("Output", geno[g],"/SeqData/"),pattern="SeqData")
        SeqData<-list()
        for (i in 1:length(HCVFiles)){ 
                df<-read.csv(paste0("Output", geno[g],"/SeqData/",HCVFiles[i]),stringsAsFactors=FALSE, row.names = 1)
                SeqData[[i]]<-df
                names(SeqData)[i]<-substr(paste(HCVFiles[i]),start=1,stop=7)
                listname<-paste0("SeqData",geno[g])
                assign(listname, SeqData)
        }
        
}

save(SeqData1A,file="SeqData1A.RData")
save(SeqData1B,file="SeqData1B.RData")
save(SeqData3A,file="SeqData3A.RData")

#load the seqdata saved as Overview1
load("SeqData1A.RData")
load("SeqData1B.RData")
load("SeqData3A.RData")

Overview2.2.3A<-Overview3A

merged.meta<-read.csv("Output_all/merged.metadata.csv", stringsAsFactors = F, row.names = 1)
mergepo<-data.frame(merged.pos=merged.meta[264:8699,"merged.pos"])

for (f in 1:3){
        SeqData1<-get(paste0("SeqData",geno[f]))
        seqdata1<-SeqData1[[1]]
        
        #remove the sites maj!=ref
        seqdata1<-seqdata1[seqdata1$MajNt==paste0(seqdata1$ref),]
        #remove the sites with total reads <1000 
        low_reads<-which(seqdata1$TotalReads<1000) 
        seqdata1[low_reads,c("a","c","g","t")]<-NA
        
        seqdata1<-seqdata1[!is.na(seqdata1$a),]
        a<-seqdata1[,c("pos", "a")]
        c<-seqdata1[,c("pos", "c")]
        g<-seqdata1[,c("pos", "g")]
        t<-seqdata1[,c("pos", "t")]
        
        
        for (i in 2:length(SeqData1)){ 
                seqdata<-SeqData1[[i]]
                seqdata<-seqdata[seqdata$MajNt==paste0(seqdata$ref),]
                low_reads<-which(seqdata$TotalReads<1000) 
                seqdata[low_reads,c("a","c","g","t")]<-NA
                seqdata<-seqdata[!is.na(seqdata$a),]
                
                a<-merge(a, seqdata[,c("pos","a")], by="pos", all = T)
                c<-merge(c, seqdata[,c("pos","c")], by="pos", all = T)
                g<-merge(g, seqdata[,c("pos","g")], by="pos", all = T)
                t<-merge(t, seqdata[,c("pos","t")], by="pos", all = T)
        }
        
        s<-length(SeqData1)
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
        colnames(dat)[1]<-paste0("org.pos.",geno[f])
        dt<-merge(dat, merged.meta, by=paste0("org.pos.",geno[f]), all.y=T )
        dt<-dt[order(dt$merged.pos),]
        write.csv(dt, paste0("Output_all/Fst_T/D_",geno[f],".csv"))
}        

#data trimming (remove extra rows -run only once)
for (i in 1:3){
        dt<-read.csv(paste0("Output_all/Fst_T/D_", geno[i],".csv"), stringsAsFactors = F, row.names = 1)
        dt<-dt[,c(7,2:6,1,8:31)]
        #dt$D[rowSums(dt[,2:5])<3000]<-NA
          
        if (i==1) dt<-dt[dt$merged.pos>263&dt$merged.pos<8679,]
        if (i==2) dt<-dt[dt$merged.pos>263&dt$merged.pos<8639,]
        if (i==3) dt<-dt[dt$merged.pos>263&dt$merged.pos<8642,]
        
        write.csv(dt,paste0("Output_all/Fst_T/D2_",geno[i],".csv") )
}

       
###### Calculate Fst between genotypes
geno<-c("1A","1B","3A")
CombG<-t(combn(geno,2))
FstList<-list()

#normalize the uneven reads between genotypes
for (f in 2:3){
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
        
        #mergepo<-data.frame(merged.pos=merged.meta[264:n,"merged.pos"])
        #df1<-merge(mergepo, df1, by="merged.pos", all.x=T)
        #df2<-merge(mergepo, df2, by="merged.pos", all.x=T)
        
        FstDF<-data.frame(merged.pos=df1$merged.pos)
        FstDF$DtT<-""
        
        #divider (m) =(m^2-m)/2   For m=2*m1, divider(m)= 2m1^2-m1
        for (k in 1:nrow(FstDF)){
                if (is.na(df1$Dt[k])|is.na(df2$Dt[k])) FstDF$DtT[k]<-NA
                else{
                        #normalize the reads between the genotypes
                        adj<-sum(df1[k,c("a","g","c","t")])/sum(df2[k,c("a","g","c","t")])

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



#####
FstList<-list()
fstfiles<-list.files("Output_all/Fst_T/", pattern=glob2rx("*^Fst*.csv*"))
for (i in 1:3){
        FstList[[i]]<-read.csv(paste0("Output_all/Fst_T/",fstfiles[i]),stringsAsFactors = F, row.names = 1)
        names(FstList)[i]<-substr(fstfiles[i], start=1, stop = 9)
}

### Plot Fst

colors2<-qualitative_hcl(6, palette="Dark3")
col2_2<-paste0(colors2,"99")

genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
genenames<-genes$Gene[1:12]

nrow(FstList[[1]]) #6556
nrow(FstList[[2]]) #5652
nrow(FstList[[3]]) #5734

d<-FstList[[2]]
min(d$Fst, na.rm=T)

labels<-substr(names(FstList), start = 5,stop=9) 

pdf("Output_all/Fst_T/Fst_across.genome.pdf", width=10, height=4)

plot(Fst~merged.pos, data=FstList[[1]], t="n", xlab='Genome position', ylab="Fst",
     ylim=c(-.05,.1), xlim=c(264,8400),col=colors2[3], cex=0.2)
points(Fst~merged.pos, data=FstList[[3]],pch=16,col=col2_2[5], ylim=c(-1,1), xlim=c(264,8400),cex=0.2)
points(Fst~merged.pos, data=FstList[[2]],pch=16,col=col2_2[3], ylim=c(-1,1), xlim=c(264,8400),cex=0.2)
points(Fst~merged.pos, data=FstList[[1]],pch=16,col=col2_2[1], ylim=c(-1,1), xlim=c(264,8400),cex=0.2)

legend("topright", legend=labels, col=colors2[c(1,3,5)], pch=16, bty="n", cex=0.8)
genes2<-genes[1:12,]

ylow<--0.046
for (j in 1:nrow(genes2)){
        xleft<- genes2$start[j]
        xright<-genes2$start[j+1]
        if (j==1){
                
                rect(-200,ylow,xright,ylow*1.3,density = NULL, angle = 45,col="white",border ="gray60",lwd=1.5)
                text(150,1.1*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else if (j==6){
                rect(xleft,ylow,xright,ylow*1.3,density = NULL, angle = 45,col="white",border ="gray60",lwd=1.5)
                text(xleft+80,1.1*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
                
        }
        else if (j==4|j==9){
                rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="gray60",lwd=1.5)
                text(xleft+50,0.8*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else if (j==12){
                rect(xleft,ylow,genes2$end[j],1.3*ylow,density = NULL, angle = 45,col="white",border ="gray60",lwd=1.5)
                text(xleft+600,1.1*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else{rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="gray60",lwd=1.5)
                text(xright-(xright-xleft)/2,1.1*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)}
}


box()
dev.off()


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
        scale_fill_gradientn(colors=colorRampBlue(64), limit=c(0,0.01))+
        theme_bw()+
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
              panel.grid.major = element_blank(),panel.border = element_blank())+
        geom_rect(aes(xmin = 1 - 0.5, xmax = 3 + 0.5, ymin = 1 - 0.5, ymax = 3 + 0.5),
                  fill = "transparent", color = "gray40", size = .3)+
        geom_text(aes(label = round(Fst, 4)))+
        geom_text(aes(label=type), nudge_y=0.2)
ggsave("Output_all/Fst_T/Fst.genotypes.heatmap.pdf",width=5, height=4.3, units='in',device='pdf')





