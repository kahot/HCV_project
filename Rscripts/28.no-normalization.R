#Calculate between-Patients and between-genotypes Fst using Dvalues
source("Rscripts/baseRscript.R")
library(grid)
library(purrr)
library(zoo)
library(tidyverse)
library(gplots)


geno<-c("1A","1B","3A")

flist1<-list.files(paste0("Output1A/Overview_D2/"),pattern="overviewD.csv")
Overview1<-list()
for (i in 1:length(flist1)){ 
        overviews<-read.csv(paste0("Output1A/Overview_D2/",flist1[i]),stringsAsFactors=FALSE, row.names = 1)
        Overview1[[i]]<-overviews
        names(Overview1)[i]<-substr(paste(flist1[i]),start=1,stop=7)
}
Overview1.1A<-Overview1

load("Overview1.1A.RData")
load("Overview1.1B.RData")
load("Overview1.3A.RData")

geno<-c("1A","1B","3A")

for (f in 1:3){
        
        Overview1<-get(paste0("Overview1.",geno[f]))
        overview1<-Overview1[[1]]
        #remove the sites with total reads <1000
        for (j in 1:nrow(overview1)){
                if (is.na(overview1$freq.Ts[j]))  {overview1$a[j]<-NA
                overview1$c[j]<-NA
                overview1$g[j]<-NA
                overview1$t[j]<-NA  } 
        }
        a<-overview1[,c("pos", "a")]
        c<-overview1[,c("pos", "c")]
        g<-overview1[,c("pos", "g")]
        t<-overview1[,c("pos", "t")]
        
        
        for (i in 2:length(Overview1)){ 
                overview<-Overview1[[i]]
                for (j in 1:nrow(overview)){
                        if (is.na(overview$freq.Ts[j]))  {overview$a[j]<-NA
                        overview$c[j]<-NA
                        overview$g[j]<-NA
                        overview$t[j]<-NA  } 
                }
                a<-merge(a, overview[,c(1,2)], by="pos", all.x = T)
                c<-merge(c, overview[,c(1,3)], by="pos", all.x = T)
                g<-merge(g, overview[,c(1,4)], by="pos", all.x = T)
                t<-merge(t, overview[,c(1,5)], by="pos", all.x = T)
        }
        s<-length(Overview1)
        a$aSum<-rowSums(a[2:(s+1)], na.rm=T)
        c$cSum<-rowSums(c[2:(s+1)], na.rm=T)
        g$gSum<-rowSums(g[2:(s+1)], na.rm=T)
        t$tSum<-rowSums(t[2:(s+1)], na.rm=T)
        
        dat<-cbind(a[,c("pos","aSum")],c[,"cSum"],g[,"gSum"],t[,"tSum"])
        colnames(dat)<-c("pos","a","c","g","t")
        dat$D<-""
        for (k in 1:nrow(dat)){
                if (rowSums(dat[k,2:5])<=3000)  dat$D[k]<-NA
                else{
                        ac<-dat[k, "a"]*dat[k,"c"]
                        ag<-dat[k, "a"]*dat[k,"g"]
                        at<-dat[k, "a"]*dat[k,"t"]
                        cg<-dat[k, "c"]*dat[k,"g"]
                        ct<-dat[k, "c"]*dat[k,"t"]
                        gt<-dat[k, "g"]*dat[k,"t"]
                        m<-(rowSums(dat[k,2:5])^2-rowSums(dat[k,2:5]))/2
                        dat$D[k]<-(ac+ag+at+cg+ct+gt)/m
                }
        }
        
        merged.meta<-read.csv("Output_all/merged.metadata.csv", stringsAsFactors = F, row.names = 1)
        dt<-merge(dat, merged.meta, by.x="pos", by.y=paste0("org.pos.",geno[f]), all=T )
        write.csv(dt, paste0("Output_all/Fst/D_",geno[f],".csv"))
}        

#data trimming (remove extra rows -run only once)
merged.meta<-read.csv("Output_all/merged.metadata.csv", stringsAsFactors = F, row.names = 1)
mergepo<-data.frame(merged.pos=merged.meta[1:8699,"merged.pos"])

for (i in 1:3){
        dt<-read.csv(paste0("Output_all/Fst/D_", geno[i],".csv"), stringsAsFactors = F, row.names = 1)
        dt<-dt[,c(7,2:5,1,6,8:32)]
        dt<-merge(dt, mergepo, by="merged.pos", all.y = T)
        #dt$D[rowSums(dt[,2:5])<3000]<-NA
        
        if (i==1) dt<-dt[dt$merged.pos>263&dt$merged.pos<8679,]
        if (i==2) dt<-dt[dt$merged.pos>263&dt$merged.pos<8639,]
        if (i==3) dt<-dt[dt$merged.pos>263&dt$merged.pos<8642,]
        
        write.csv(dt,paste0("Output_all/Fst/D2_",geno[i],".csv") )
}


###### Calculate Fst between the samples within genotypes
geno<-c("1A","1B","3A")
CombG<-t(combn(geno,2))
FstList<-list()

for (f in 2:3){
        g1<-CombG[f,1]
        g2<-CombG[f,2]
        fname<-paste0(g1,"-",g2)
        print(fname)
        
        df1<-read.csv(paste0("Output_all/Fst/D2_", g1,".csv"), stringsAsFactors = F, row.names = 1)
        df2<-read.csv(paste0("Output_all/Fst/D2_", g2,".csv"), stringsAsFactors = F, row.names = 1)
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
        
        
        FstDF<-data.frame(merged.pos=df1$merged.pos)
        FstDF$DT<-""
        for (k in 1:nrow(FstDF)){
                if (is.na(df1$D[k])|is.na(df2$D[k])) FstDF$DT[k]<-NA
                else{
                        ac<-(df1[k, "a"]+df2[k, "a"])*(df1[k,"c"]+df2[k,"c"])
                        ag<-(df1[k, "a"]+df2[k, "a"])*(df1[k,"g"]+df2[k,"g"])
                        at<-(df1[k, "a"]+df2[k, "a"])*(df1[k,"t"]+df2[k,"t"])
                        cg<-(df1[k, "c"]+df2[k, "c"])*(df1[k,"g"]+df2[k,"g"])
                        ct<-(df1[k, "c"]+df2[k, "c"])*(df1[k,"t"]+df2[k,"t"])
                        gt<-(df1[k, "g"]+df2[k, "g"])*(df1[k,"t"]+df2[k,"t"])
                        m<-as.numeric(((rowSums(df1[k,2:5])+rowSums(df2[k,2:5]))^2-(rowSums(df1[k,2:5])+rowSums(df2[k,2:5])))/2)
                        FstDF$DT[k]<-(ac+ag+at+cg+ct+gt)/m
                }
        }
        
        
        FstDF$DT<-as.numeric(FstDF$DT)
        
        for (k in 1:nrow(FstDF)){
                FstDF$D.ave[k]<-mean(c(df1$D[k],df2$D[k]))
                FstDF$Fst[k]<-(FstDF$DT[k]- FstDF$D.ave[k])/ FstDF$DT[k]
        }
        
        tname<-paste0("Fst.",fname)
        assign(tname,FstDF)
        
        FstList[[f]]<-FstDF
        names(FstList)[f]<-tname
        
        write.csv(FstDF,paste0("Output_all/Fst/",tname,".csv"))
}

### Plot Fst


colors2<-paste0(cols2,"99")

genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
genenames<-genes$Gene[1:12]

DF<-FstList[[1]]
hist(DF$Fst)



labels<-substr(names(FstList), start = 5,stop=9) 

pdf("Output_all/Fst/Fst_across.genome.pdf", width=10, height=4)

plot(Fst~merged.pos, data=FstList[[3]],pch=16,xlab='Genome position', ylab="Fst per site",
     ylim=c(-1.3,1), xlim=c(264,8400),col=colors2[3], cex=0.3)
points(Fst~merged.pos, data=FstList[[1]],pch=16,col=colors2[1], ylim=c(-1,1), xlim=c(264,8400),cex=0.3)
points(Fst~merged.pos, data=FstList[[2]],pch=16,col=colors2[2], ylim=c(-1,1), xlim=c(264,8400),cex=0.3)
legend("topright", legend=labels, col=cols2, pch=16, bty="n", cex=0.8)


genes2<-genes[1:12,]
genes2$Gene[6]<-"NS1"
ylow<-1
for (j in 1:nrow(genes2)){
        xleft<- genes2$start[j]
        xright<-genes2$start[j+1]
        if (j==1){
                
                rect(-200,ylow,xright,ylow*1.6,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(150,1.2*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else if (j==6){
                rect(xleft,ylow,xright,ylow*1.6,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(xleft+80,1.2*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
                
        }
        else if (j==4|j==9){
                rect(xleft,ylow,xright,1.6*ylow,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(xleft+50,2.0*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else if (j==12){
                rect(xleft,ylow,genes2$end[j],1.6*ylow,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(xleft+600,1.2*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else{rect(xleft,ylow,xright,1.6*ylow,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(xright-(xright-xleft)/2,1.2*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)}
}


box()



plot(Fst~merged.pos, data=FstList[[2]],pch=16,xlab='Genome position', ylab="Fst per site",
     ylim=c(0,1), xlim=c(264,8400),col=colors2[1], cex=0.3)

plot(Fst~merged.pos, data=FstList[[3]],pch=16,xlab='Genome position', ylab="Fst per site",
     ylim=c(0,1), xlim=c(264,8400),col=colors2[1], cex=0.3)




plot(DT~merged.pos, data=FstList[[3]],pch=16,xlab='Genome position', ylab="Fst per site",
     ylim=c(-0.1,1), xlim=c(264,8400),col=colors2[3], cex=0.3)
