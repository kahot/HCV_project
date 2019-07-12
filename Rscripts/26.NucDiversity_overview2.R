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

source("Rscripts/baseRscript.R")

dir.create("Output3A/Overview_D2/")
dir.create("Output1B/Overview_D2/")
dir.create("Output1A/Overview_D2/")

geno<-c("1A","1B","3A")
f=2
f=1

for (f in 1:1){
        flist1<-list.files(paste0("Output",geno[f],"/Overview2/"),pattern="overview2.csv")
        flist2<-list.files(paste0("Output",geno[f],"/SeqData/"),pattern="SeqData")
        
        Overview1<-list()
        #for (i in 1:length(flist1)){ 
        for (i in 1:30){ 
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


