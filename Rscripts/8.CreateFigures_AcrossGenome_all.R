library(purrr)
library(tidyverse)
library(zoo)
library(plotrix)
library(sfsmisc)
#Script to analyse the frequency data and associate with features. 


#Prep data
source("Rscripts/baseRscript.R")

HCVFiles_overview<-list.files("Output/Overview2/",pattern="overview2.csv")

Overview_summary<-list()
for (i in 1:length(HCVFiles_overview)){ 
        overviews<-read.csv(paste0("Output/Overview2/",HCVFiles_overview[i]),stringsAsFactors=FALSE)
        overviews<-overviews[,-1]
        Overview_summary[[i]]<-overviews
        names(Overview_summary)[i]<-substr(paste(HCVFiles_overview[i]),start=1,stop=7)
}

source("Rscripts/MutationFreqSum.R")


MutFreq_Ts<-list()
MutFreq_tv1<-list()
MutFreq_tv2<-list()
MutFreq_tvs<-list()
MutFreq_all<-list()

for (i in 1:length(Overview_summary)){
        dat<-Overview_summary[[i]]
        filename<-names(Overview_summary)[i]
        
        MutFreq_Ts[[i]]<-dat[,c("pos","freq.Ts")] 
        MutFreq_tv1[[i]]<-dat[,c("pos","freq.transv1")] 
        MutFreq_tv2[[i]]<-dat[,c("pos","freq.transv2")] 
        MutFreq_tvs[[i]]<-dat[,c("pos","freq.transv")] 
        MutFreq_all[[i]]<-dat[,c("pos","freq.mutations")] 
        
        names(MutFreq_Ts)[i]<-filename
        names(MutFreq_tv1)[i]<-filename
        names(MutFreq_tv2)[i]<-filename
        names(MutFreq_tvs)[i]<-filename
        names(MutFreq_all)[i]<-filename
        
        
}
#assign column names for the list
for (i in 1:length(MutFreq_Ts)) {
        colnames(MutFreq_Ts[[i]])<-c("pos",paste0(names(MutFreq_Ts[i])))
        colnames(MutFreq_tv1[[i]])<-c("pos",paste0(names(MutFreq_tv1[i])))
        colnames(MutFreq_tv2[[i]])<-c("pos",paste0(names(MutFreq_tv2[i])))
        colnames(MutFreq_tvs[[i]])<-c("pos",paste0(names(MutFreq_tvs[i])))
        colnames(MutFreq_all[[i]])<-c("pos",paste0(names(MutFreq_all[i])))
}

TsMutFreq<-MutFreq_Ts%>% purrr::reduce(full_join, by='pos')
Tv1.MutFreq<-MutFreq_tv1 %>% purrr::reduce(full_join, by='pos')
Tv2.MutFreq<-MutFreq_tv2 %>% purrr::reduce(full_join, by='pos')
Tvs.MutFreq<-MutFreq_tvs %>% purrr::reduce(full_join, by='pos')
AllMutFreq<-MutFreq_all %>% purrr::reduce(full_join, by='pos')


s<-length(Overview_summary)
mean(rowMeans(TsMutFreq[2:(s+1)],na.rm=T),na.rm=T) #[1] 0.008879231
mean(rowMeans(Tvs.MutFreq[2:(s+1)],na.rm=T),na.rm=T) # 0.001690682
mean(rowMeans(AllMutFreq[2:(s+1)],na.rm=T),na.rm=T) # 0.01056991

#####
# Write csv files
write.csv(TsMutFreq, "Output/MutFreq/TsMutFreq_maj_summary.csv")
write.csv(Tv1.MutFreq, "Output/MutFreq/Tv1MutFreq_maj_summary.csv")
write.csv(Tv2.MutFreq, "Output/MutFreq/Tv2MutFreq_maj_summary.csv")
write.csv(Tvs.MutFreq, "Output/MutFreq/TvsMutFreq_maj_summary.csv")
write.csv(AllMutFreq, "Output/MutFreq/AllMutFreq_maj_summary.csv")


#############################
##### read the saved csv files
TsMutFreq<-read.csv("Output/MutFreq/TsMutFreq_maj_summary.csv",stringsAsFactors = F)
Tv1.MutFreq<-read.csv("Output/MutFreq/Tv1MutFreq_maj_summary.csv",stringsAsFactors = F)
Tv2.MutFreq<-read.csv("Output/MutFreq/Tv2MutFreq_maj_summary.csv",stringsAsFactors = F)
Tvs.MutFreq<-read.csv("Output/MutFreq/TvsMutFreq_maj_summary.csv",stringsAsFactors = F)
AllMutFreq<-read.csv("Output/MutFreq/AllMutFreq_maj_summary.csv",stringsAsFactors = F)

AllMutFreq<-AllMutFreq[,-1]
TsMutFreq<-TsMutFreq[,-1]
Tvs.MutFreq<-Tvs.MutFreq[,-1]
###################################



##Create SNV (mutation) frequency summary 
#total number of samples 
HCVFiles_overview<-list.files("Output/Overview2/",pattern="overview2.csv")
s<-length(HCVFiles_overview)

genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)

MFsummary<-data.frame("Mean"=matrix(nrow=3))
rownames(MFsummary)<-c("AllMF","Transition","Transversion")
for (i in 1:3){
        if (i==1) {SNVFreq<-AllMutFreq; title<-'Average Total Mutation Freq (1A)';name<-"Total_Mutations_1A";yax<-"Average mutation frequency"; ylow<-0.0001;yhigh<-1}
        if (i==2) {SNVFreq<-TsMutFreq; title<-'Average Transition Mutation Freq (1A)';name<-"Transition_1A";yax<- "Average transition mutation frequency";ylow<-0.0001;yhigh<-1}
        if (i==3) {SNVFreq<-Tvs.MutFreq; title<-'Average Transversion Mutation Freq (1A)';name<-"Transversion_1A";yax<- "Average transversion mutation frequency";ylow<-0.00001;yhigh<-0.1}        
         
        #Count the number of NAs per row
        SNVFreq$NofNA<-apply(SNVFreq[2:(s+1)], 1, function(x)sum(is.na(x))) #8371 sites
        SNVFreq1<-SNVFreq[which(SNVFreq$NofNA<0.5*s),]
        SNVFreq1$mean<-apply(SNVFreq1[2:(s+1)], 1, mean, na.rm=T)
        print(nrow(SNVFreq1))
        
        #2. Filter out sites with only one data points rows
        SNVFreq2<-SNVFreq[-which(SNVFreq$NofNA==(s-1)|SNVFreq$NofNA==s),] 
        SNVFreq2$mean<-apply(SNVFreq2[2:(s+1)], 1, mean, na.rm=T)
        
        DFmean<-SNVFreq1[,c("pos","mean")]
        write.csv(DFmean, paste0("Output/MutFreq/Ave.MF.",name,'.csv'))
        
        
        startnuc<-SNVFreq1$pos[1]
        endnuc<-SNVFreq1$pos[nrow(SNVFreq1)] #8635
        
        MFsummary$Mean[i]<-mean(SNVFreq1$mean)
        MFsummary$range.low[i]<-range(SNVFreq1$mean)[1]
        MFsummary$range.high[i]<-range(SNVFreq1$mean)[2]
        MFsummary$SE[i]<-std.error(SNVFreq1$mean)
        MFsummary$SD[i]<-sd(SNVFreq1$mean)
        
        
        MFsummary$Mean2[i]<-mean(SNVFreq2$mean)
        MFsummary$range.low2[i]<-range(SNVFreq2$mean)[1]
        MFsummary$range.high2[i]<-range(SNVFreq2$mean)[2]
        MFsummary$SE2[i]<-std.error(SNVFreq2$mean)
        MFsummary$SD2[i]<-sd(SNVFreq2$mean)
       
        n<-data.frame("pos"=c(startnuc:(endnuc-20)))
        SNVFreqs<-merge(n,SNVFreq1,by="pos",all.x=T)
        
        cols <- c("#66CCEE","#228833","#CCBB44","#EE6677","#AA3377","#4477AA","#BBBBBB")
        
        #plot the average SNV frequency across the genome (based on H77)
        pdf(paste0("Output/SummaryFigures/",title,".pdf"),width=15,height=7.5)
        plot(mean~pos, data=SNVFreqs,t="n",log='y',yaxt='n',xlab='Genome position (H77)',ylab=paste0(yax),
             main=paste0(title),ylim=c(ylow,yhigh),xlim=c(340,8500))
        eaxis(side = 2, at = 10^((-0):(-(5))), cex=2)
        for(k in 1:5){abline(h = 1:10 * 10^(-k), col = "gray80")}
        points(mean~pos, data=SNVFreqs,pch=20,col="#66CCEE",cex=0.5)
        
        #add rolling average
        roll100<-rollmean(SNVFreqs$mean, k=100, na.rm=T)
        SNVFreqs$roll100<-c(rep(NA, times=50),roll100,c(rep(NA, times=49)))
        lines(roll100~pos,data=SNVFreqs, col="#4477AA",lwd=1.5)
        
        
        for (j in 1:(nrow(genes)-1)){
                xleft<-genes$start[j]
                xright<-genes$start[j+1]
                if (j==1){
                        rect(0,ylow,xright,ylow*1.4,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                        text(200,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                }
                else if (j==6){
                        rect(xleft,ylow,xright,ylow*1.4,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                        text(xleft+80,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                        
                }
                else if (j==4|j==9){
                        rect(xleft,ylow,xright,1.4*ylow,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                        text(xleft+50,.9*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                }
                else{rect(xleft,ylow,xright,1.4*ylow,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                        text(xleft+200,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)}
        }
        
        
        box()
        dev.off()
        print(title)
        print(mean(SNVFreqs$mean, na.rm=T))
        
        
}

write.csv(MFsummary,"Output/SummaryStats/MutFreqSummary_Table.csv")

## Create a condensed figure for a poster

SNVFreq<-AllMutFreq
yax<-"Average mutation frequency"; ylow<-0.0001;yhigh<-1
SNVFreq$mean<-apply(SNVFreq[2:(s+1)], 1, mean, na.rm=T)
SNVFreq$NofNA<-apply(SNVFreq[2:(s+1)], 1, function(x)sum(is.na(x))) #8371 sites
SNVFreq1<-SNVFreq[which(SNVFreq$NofNA<0.5*s),]
print(nrow(SNVFreq1))
startnuc<-SNVFreq1$pos[1]
endnuc<-8609

MFsummary$Mean[i]<-mean(SNVFreq1$mean)
MFsummary$range.low[i]<-range(SNVFreq1$mean)[1]
MFsummary$range.high[i]<-range(SNVFreq1$mean)[2]
MFsummary$SE[i]<-std.error(SNVFreq1$mean)

n<-data.frame("pos"=c(startnuc:endnuc))
SNVFreqs<-merge(n,SNVFreq1,by="pos",all.x=T)

cols <- c("#66CCEE","#228833","#CCBB44","#EE6677","#AA3377","#4477AA","#BBBBBB")

par(mar=c(4,4,1,1))

pdf(paste0("Output/SummaryFigures/MF_acrossGenome1.pdf"),width=12,height=4.6)
plot(mean~pos, data=SNVFreqs,t="n",log='y',yaxt='n', xlab='',ylab='',
     ylim=c(ylow,yhigh),xlim=c(340,8500), cex=3)
eaxis(side = 2, at = 10^((-0):(-(5))), cex=3)
for(k in 1:5){abline(h = 1:10 * 10^(-k), col = "gray80")}
points(mean~pos, data=SNVFreqs,pch=20,col="#66CCEE",cex=0.5)

mtext("Genome position", side=1, line=2.5, cex=1.4)
mtext("Average mutation frequency", side=2, line=2.5, cex=1.4)
#add rolling average
roll100<-rollmean(SNVFreqs$mean, k=100, na.rm=T)
SNVFreqs$roll100<-c(rep(NA, times=50),roll100,c(rep(NA, times=49)))
lines(roll100~pos,data=SNVFreqs, col="#4477AA",lwd=1.5)
        
        
for (j in 1:(nrow(genes)-1)){
        xleft<-genes$start[j]
        xright<-genes$start[j+1]
        if (j==1){
                rect(0,ylow,xright,ylow*1.6,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(200,1.24*ylow,paste0(genes$Gene[j]),col="black", cex=.8)
        }
       # else if (j==6){
        #        rect(xleft,ylow,xright,ylow*1.6,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
         #       text(xleft+80,1.24*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                
        #}
        else if (j==4|j==6|j==9){
                rect(xleft,ylow,xright,1.6*ylow,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(xleft+50,.82*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        }
        else{rect(xleft,ylow,xright,1.6*ylow,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(xleft+200,1.24*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)}
}

box()
dev.off()




########################################
### Plot mutation freq. across the genome based on the mutation types 
dat<-Overview_summary[[3]]
mutationtypes<-dat[,c("pos","MajNt","ref","Type","Type.tv1","Type.tv2")]
muttypes2<-dat[,c("pos","MajNt","Type","Type.tv1","Type.tv2","makesCpG","makesCpG.tv1","makesCpG.tv2")]

t1<-TsMutFreq
t1$mean<-apply(t1[2:(s+1)], 1, mean, na.rm=T)
t1$NofNA<-apply(t1[2:(s+1)], 1, function(x)sum(is.na(x)))
t2<-t1[which(t1$NofNA<0.5*s),]


mfs<-merge(mutationtypes,t2,by='pos')

maxpos<-mfs$pos[nrow(mfs)]
n<-data.frame("pos"=c(1:maxpos))
MF3<-merge(n,mfs,by="pos",all.x=T)

pdf(paste0("Output/SummaryFigures/Transition_Mut_Freq_based_on_Mutation_Types.pdf"),width=15,height=7.5)
maxnuc=8600
par(mar = c(3,5,1,2))
#selcoeffcolumn <-SC3$mean 

plot(MF3$pos[1:maxnuc],MF3$mean[1:maxnuc],
     log="y", ylab="Average mutation frequency",cex.lab=1.4,
     yaxt="n", xlab="",xaxt='n',
     col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-4,0.1),xlim=c(340,maxnuc))
axis(1,at=c(seq(500,8500,by=1000)), labels=c(seq(500,8500,by=1000)))
eaxis(side = 2, at = 10^((-1):(-(4))), cex=2)

for(i in 2:4){abline(h = 1:10 * 10^(-i), col = "gray41")}

for (i in 1:maxnuc){
        c=0; co = 1
        if (is.na(MF3$Type[i])==T) next
        if (MF3$Type[i]=="stop") {c=1;p=21}
        if (MF3$Type[i]=="syn") {c=cols[3];p=21}
        if (MF3$Type[i]=="nonsyn"&MF3$ref[i]%in%c("c","g")) {c=cols[4];p=21}
        if (MF3$Type[i]=="nonsyn"&MF3$ref[i]%in%c("a","t")) {c=cols[5];p=21}
        if (c!=0) points(MF3$pos[i],MF3$mean[i],pch=p,col='gray30',lwd=0.5,
                         bg=rgb(red=col2rgb(c)[1]/255,
                                green=col2rgb(c)[2]/255,
                                blue=col2rgb(c)[3]/255,
                                maxColorValue = 1,alpha=0.8),cex=.6)
}
#Add legend
legpos=7000; legposV=0.13
rect(legpos, 0.5*legposV, (legpos+1800), 1.05*legposV, density = NULL, angle = 45,col=alpha("white",1))
points((legpos+100),legposV*0.8,pch=21,bg=1,col=1,cex=1)
text((legpos+150),legposV*0.8,"Nonsense",adj=0, cex=1)
points((legpos+900),legposV*0.8,pch=21,bg=cols[3],col=1,cex=1)
text((legpos+1050),legposV*0.8,"Syn",adj=0, cex=1)

points((legpos+100),legposV*0.6,pch=21,bg=cols[5],col=1,cex=1)
text((legpos+150),legposV*0.6,"Non-syn, A/T",adj=0, cex=1)
points((legpos+900),legposV*0.6,pch=21,bg=cols[4],col=1,cex=1)
text((legpos+1050),legposV*0.6,"Non-syn, C/G",adj=0, cex=1)



ylow=0.000275
for (j in 1:(nrow(genes)-1)){
        xleft<-genes$start[j]
        xright<-genes$start[j+1]
        if (j==1){
                rect(0,ylow,xright,ylow*1.2,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(200,1.1*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        }
        else if (j==6){
                rect(xleft,ylow,xright,ylow*1.2,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(xleft+80,1.1*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                
        }
        else if (j==4|j==9){
                rect(xleft,ylow,xright,1.2*ylow,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(xleft+50,1.3*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        }
        else{rect(xleft,ylow,xright,1.2*ylow,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(xleft+200,1.1*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)}
}

box()

dev.off()



#2) ############
#Transversion mutationsfigure
#Format the transversion mutation types to syn, nonsyn, mixed, and nonsesne.
tv1<-Tv1.MutFreq
tv2<-Tv2.MutFreq
s<-length(Overview_summary)
tv1$mean<-apply(tv1[2:(s+1)], 1, mean, na.rm=T)
tv2$mean<-apply(tv2[2:(s+1)], 1, mean, na.rm=T)
tv1$NofNA<-apply(tv1[2:(s+1)], 1, function(x)sum(is.na(x)))
tv2$NofNA<-apply(tv2[2:(s+1)], 1, function(x)sum(is.na(x)))

tv1.1<-tv1[which(tv1$NofNA<0.5*s),]
tv2.1<-tv2[which(tv2$NofNA<0.5*s),]
tv1.2<-tv1.1[,c("pos","mean")]
tv2.2<-tv2.1[,c("pos","mean")]
colnames(tv2.2)[2]<-"mean2"

tv1.3<-merge(mutationtypes,tv1.2,by='pos')
tv2.3<-merge(mutationtypes,tv2.2,by='pos')
MFtv<-merge(tv1.3,tv2.2,by='pos')

MFtv$Type.tvs<-""
for (i in 1:nrow(MFtv)){ 
        if (is.na(MFtv$Type.tv1[i])|is.na(MFtv$Type.tv2[i])) MFtv$Type.tvs[i]<-NA
        else 
               {if (MFtv$Type.tv1[i] == MFtv$Type.tv2[i]) MFtv$Type.tvs[i]<-MFtv$Type.tv1[i]
                if (MFtv$Type.tv1[i] != MFtv$Type.tv2[i]) MFtv$Type.tvs[i]<-'mixed'}
}        

######
#2.1 plot the summary trnasversion mutations 

maxpos<-MFtv$pos[nrow(MFtv)]
n<-data.frame("pos"=c(1:maxpos))
MFtv.1<-merge(n,MFtv,by="pos",all.x=T)


col80<-c("#FFFFFFCC","#228833CC","#CCBB44CC","#EE6677CC","#AA3377CC")

pdf(paste0("Output/SummaryFigures/Transv_Mut_Freq_based_on_Mutation_Types.pdf"),width=15,height=7.5)
maxnuc=8600
par(mar = c(3,5,1,2))

plot(MFtv.1$pos[1:maxnuc],MFtv.1$mean[1:maxnuc],
     log="y", ylab="Average mutation frequency",cex.lab=1.4,
     yaxt="n", xlab="",xaxt='n',
     col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-6,0.1),xlim=c(340,maxnuc))
axis(1,at=c(seq(500,8500,by=1000)), labels=c(seq(500,8500,by=1000)))
eaxis(side = 2, at = 10^((-1):(-(6))), cex=2)

for(i in 1:6){abline(h = 1:10 * 10^(-i), col = "gray60")}

for (i in 1:maxnuc){
        if (is.na(MFtv.1$Type.tvs[i])) next
        if (MFtv.1$Type.tvs[i]=="stop") {c=1;p=21}
        if (MFtv.1$Type.tvs[i]=="syn") {c=col80[3];p=21}
        if (MFtv.1$Type.tvs[i]=="nonsyn"&MFtv.1$ref[i]%in%c("c","g")) {c=col80[4];p=21}
        if (MFtv.1$Type.tvs[i]=="nonsyn"&MFtv.1$ref[i]%in%c("a","t")) {c=col80[5];p=21}
        if (MFtv.1$Type.tvs[i]=="mixed"){c=col80[2];p=21}
        points(MFtv.1$pos[i],mean(MFtv.1$mean[i],MFtv.1$mean2[i]),pch=p,
               col='gray30',lwd = .5, bg=c,cex=.7)
        
}

#Add legend
legpos=7700; legposV=0.05
rect(legpos, .29*legposV, (legpos+1000), 1.7*legposV, density = NULL, angle = 45,col=alpha("white",1))
points((legpos+100),legposV/0.7,pch=21,bg=1,col=col80[1],cex=1)
text((legpos+150),legposV/0.7,"Nonsense",adj=0, cex=1)
points((legpos+100),legposV,pch=21,bg=col80[4],col=1,cex=1)
text((legpos+150),legposV,"Non-syn, C/G",adj=0, cex=1)
points((legpos+100),legposV*0.7,pch=21,bg=col80[5],col=1,cex=1)
text((legpos+150),legposV*0.7,"Non-syn, A/T",adj=0, cex=1)
points((legpos+100),legposV*0.49,pch=21,bg=col80[2],col=1,cex=1)
text((legpos+150),legposV*0.49,"Mixed",adj=0, cex=1)
points((legpos+100),legposV*0.35,pch=21,bg=cols[3],col=1,cex=1)
text((legpos+150),legposV*0.35,"Syn",adj=0, cex=1)

dev.off()

#2.2. plot trnasversion mutations separately (Type 1 and Type 2)


maxpos<-tv1.3$pos[nrow(tv1.3)]
n<-data.frame("pos"=c(1:maxpos))
tv1.4<-merge(n,tv1.3,by="pos",all.x=T)
tv2.4<-merge(n,tv2.3,by="pos",all.x=T)

col80<-c("#FFFFFFCC","#228833CC","#CCBB44CC","#EE6677CC","#AA3377CC")

pdf(paste0("Output/SummaryFigures/Transv_MutFreq_Mutation_Types_separate.pdf"),width=15,height=7.5)
maxnuc=8600
par(mar = c(3,5,1,2))

plot(tv1.4$pos[1:maxnuc],tv1.4$mean[1:maxnuc],
     log="y", ylab="Average mutation frequency",cex.lab=1.4,
     yaxt="n", xlab="",xaxt='n',
     col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-6,.1),xlim=c(340,maxnuc))
axis(1,at=c(seq(500,8500,by=1000)), labels=c(seq(500,8500,by=1000)))
eaxis(side = 2, at = 10^((-0):(-(6))), cex=2)

for(i in 1:6){abline(h = 1:10 * 10^(-i), col = "gray60")}

for (i in 1:maxnuc){
        if (is.na(tv1.4$Type.tv1[i])==T) next
        if (tv1.4$Type.tv1[i]=="stop") {c=1;p=21}
        if (tv1.4$Type.tv1[i]=="syn") {c=col80[3];p=21}
        if (tv1.4$Type.tv1[i]=="nonsyn"&tv1.4$ref[i]%in%c("c","g")) {c=col80[4];p=21}
        if (tv1.4$Type.tv1[i]=="nonsyn"&tv1.4$ref[i]%in%c("a","t")) {c=col80[5];p=21}
        if (is.na(tv2.4$Type.tv2[i])==T) next
        if (tv2.4$Type.tv2[i]=="stop") {c=1;p=21}
        if (tv2.4$Type.tv2[i]=="syn") {c=col80[3];p=21}
        if (tv2.4$Type.tv2[i]=="nonsyn"&tv2.4$ref[i]%in%c("c","g")) {c=col80[4];p=21}
        if (tv2.4$Type.tv2[i]=="nonsyn"&tv2.4$ref[i]%in%c("a","t")) {c=col80[5];p=21}
        
        points(tv1.4$pos[i],tv1.4$mean[i],pch=p,
               col='gray30',lwd = .5, bg=c,cex=.7)
        points(tv2.4$pos[i],tv2.4$mean[i],pch=p,
               col='gray30',lwd = .5, bg=c,cex=.7)
}

#Add legend
legpos=7700; legposV=0.05
rect(legpos, .29*legposV, (legpos+1000), 1.2*legposV, density = NULL, angle = 45,col=alpha("white",1))
points((legpos+100),legposV,pch=21,bg=1,col=col80[1],cex=1)
text((legpos+150),legposV,"Nonsense",adj=0, cex=1)
points((legpos+100),legposV*0.7,pch=21,bg=col80[4],col=1,cex=1)
text((legpos+150),legposV*0.7,"Non-syn, C/G",adj=0, cex=1)
points((legpos+100),legposV*0.49,pch=21,bg=col80[5],col=1,cex=1)
text((legpos+150),legposV*0.49,"Non-syn, A/T",adj=0, cex=1)
points((legpos+100),legposV*0.35,pch=21,bg=cols[3],col=1,cex=1)
text((legpos+150),legposV*0.35,"Syn",adj=0, cex=1)

ylow=0.0000024
for (j in 1:(nrow(genes)-1)){
        xleft<-genes$start[j]
        xright<-genes$start[j+1]
        if (j==1){
                rect(0,ylow,xright,ylow*1.3,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(200,1.15*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)}
        else if (j==6){
                rect(xleft,ylow,xright,ylow*1.3,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(xleft+80,1.15*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)}
        else if (j==4|j==9){
                rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(xleft+50,1.4*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)}
        else{rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="#4477AA",lwd=1.5)
                text(xleft+200,1.15*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)}
}

box()


dev.off()





##### Test if mean frequenceis differ between mutation types
#testing the mean of means
#1.Mutation frequencies (all together)  

for (i in c("t2","tv1.2","tv2.2")) {
        dat<-get(i) 
        FreqData<-merge(muttypes2,dat,by='pos')
        if (i=="t2") {type<-3;cpg<-6; filename<-"transitions"}
        if (i=="tv1.2") {type<-4;cpg<-7; filename<-"transv1"}
        if (i=="tv2.2") {type<-5;cpg<-8; filename<-"transv2"}
        
        types<-FreqData[type]
        typevector1<-types=="syn"
        typevector2<-types=="nonsyn"
        typevector3<-types=="stop"
        mfSyn<-FreqData$mean[typevector1==T]
        mfNonSyn<-FreqData$mean[typevector2==T]
        mfStop<-FreqData$mean[typevector3==T]
        cpgmake<-FreqData[cpg]
        cpgvector1<-cpgmake==0
        cpgvector2<-cpgmake==1
        
        result1<-wilcox.test(mfSyn, mfNonSyn,alternative = "greater", paired = FALSE) #p-value < 2.2e-16
        result2<-wilcox.test(mfNonSyn,mfStop,alternative = "greater", paired = FALSE) # p-value = 0.952
        
        #Test if CpG making mutation frequencies are lower than non-CpGmaking 
        mfSyn.CpG<-FreqData$mean[typevector1==T&cpgvector2==T]
        mfSyn.nonCpG<-FreqData$mean[typevector1==T&cpgvector1==T]
        mfNonSyn.CpG<-FreqData$mean[typevector2==T&cpgvector2==T]
        mfNonSyn.nonCpG<-FreqData$mean[typevector2==T&cpgvector1==T]
        
        result3<-wilcox.test(mfSyn.CpG, mfSyn.nonCpG, alternative = "greater", paired = FALSE) #p-value < 2.2e-16
        result4<-wilcox.test(mfNonSyn.CpG, mfNonSyn.nonCpG, alternative = "greater", paired = FALSE)  #p-value < 2.2e-16 
        
        #Test whether synonymous C to T  and G to A mutations are more costly than A to G 
        CostSynNoCpGA2G<-FreqData$mean[typevector1==T&cpgvector1==T&FreqData$MajNt=="a"]
        CostSynNoCpGC2T<-FreqData$mean[typevector1==T&cpgvector1==T&FreqData$MajNt=="c"]
        CostSynNoCpGG2A<-FreqData$mean[typevector1==T&cpgvector1==T&FreqData$MajNt=="g"]
        result5<-wilcox.test(CostSynNoCpGC2T,CostSynNoCpGA2G,alternative = "greater", paired = FALSE)  #p-value =1
        result6<-wilcox.test(CostSynNoCpGG2A,CostSynNoCpGA2G,alternative = "greater", paired = FALSE)   #p-value =1
        
        #Test whether non-synonymous C to T  muts and G to A , are more costly than A to G 
        CostnonSynNoCpGA2G<-FreqData$mean[typevector2==T&cpgvector1==T&FreqData$MajNt=="a"]
        CostnonSynNoCpGC2T<-FreqData$mean[typevector2==T&cpgvector1==T&FreqData$MajNt=="c"]
        CostnonSynNoCpGG2A<-FreqData$mean[typevector2==T&cpgvector1==T&FreqData$MajNt=="g"]
        result7<-wilcox.test(CostnonSynNoCpGC2T,CostnonSynNoCpGA2G,alternative = "greater", paired = FALSE) #p-value =1
        result8<-wilcox.test(CostnonSynNoCpGG2A,CostnonSynNoCpGA2G,alternative = "greater", paired = FALSE) #p-value =1
        
        WilcoxTest.results<-data.frame(matrix(ncol=2,nrow=6))
        colnames(WilcoxTest.results)<-c("test","P.value")
        for (r in 1:6){
                result<-get(paste0('result',r))
                WilcoxTest.results$test[r]<-result[[7]]
                WilcoxTest.results$P.value[r]<-result[[3]]
                
        }
        write.csv(WilcoxTest.results,paste0("Output/SummaryStats/WilcoxTestResults_MajMutFreq_",filename,".csv"))
        
}



#2.Mutation frequencies - base by base 
source("Rscripts/MutationFreqSum.R")
dat<-Overview_summary[[3]]
muttypes2<-dat[,c("pos","MajNt","Type","Type.tv1","Type.tv2","makesCpG","makesCpG.tv1","makesCpG.tv2")]

WilcoxTest.results<-data.frame(Nt="",test="",P.value="")
WilcoxTest.results<-data.frame(matrix(ncol=3,nrow=2))
colnames(WilcoxTest.results)<-c("nt","test","P.value")

for (i in c("A","T")) {
        scpg<-get(paste0(i,"_syn_CpG")) 
        sncpg<-get(paste0(i,"_syn_NonCpG")) 
        nscpg<-get(paste0(i,"_nonsyn_CpG")) 
        nsncpg<-get(paste0(i,"_nonsyn_NonCpG")) 
        fname<-i
        result1<-wilcox.test(scpg, sncpg, alternative = "greater", paired = FALSE) 
        result2<-wilcox.test(nscpg,nsncpg,alternative = "greater", paired = FALSE) # p-value = 0.952
        
        for (r in 1:2){
                result<-get(paste0('result',r))
                WilcoxTest.results$nt[r]<-i
                WilcoxTest.results$test[r]<-result[[7]]
                WilcoxTest.results$P.value[r]<-result[[3]]
                }
        write.csv(WilcoxTest.results,paste0("Output/SummaryStats/Wilcox_Maj_TsMutFreq_",fname,".csv"))
}       
 


for (i in c("A","G")) {
        scpg<-get(paste0(i,"_tv1_syn_CpG")) 
        sncpg<-get(paste0(i,"_tv1_syn_NonCpG")) 
        nscpg<-get(paste0(i,"_tv1_nonsyn_CpG")) 
        nsncpg<-get(paste0(i,"_tv1_nonsyn_NonCpG")) 
        fname<-paste0(i,"_Tv")
        result1<-wilcox.test(scpg, sncpg, alternative = "greater", paired = FALSE) 
        result2<-wilcox.test(nscpg,nsncpg,alternative = "greater", paired = FALSE) # p-value = 0.952
        
        for (r in 1:2){
                result<-get(paste0('result',r))
                WilcoxTest.results$nt[r]<-i
                WilcoxTest.results$test[r]<-result[[7]]
                WilcoxTest.results$P.value[r]<-result[[3]]
        }
        write.csv(WilcoxTest.results,paste0("Output/SummaryStats/Wilcox_Maj_TsMutFreq_",fname,".csv"))
}



for (i in c("C","T")) {
        scpg<-get(paste0(i,"_tv2_syn_CpG")) 
        sncpg<-get(paste0(i,"_tv2_syn_NonCpG")) 
        nscpg<-get(paste0(i,"_tv2_nonsyn_CpG")) 
        nsncpg<-get(paste0(i,"_tv2_nonsyn_NonCpG")) 
        fname<-paste0(i,"_Tv")
        result1<-wilcox.test(scpg, sncpg, alternative = "greater", paired = FALSE) 
        result2<-wilcox.test(nscpg,nsncpg,alternative = "greater", paired = FALSE) # p-value = 0.952
        
        for (r in 1:2){
                result<-get(paste0('result',r))
                WilcoxTest.results$nt[r]<-i
                WilcoxTest.results$test[r]<-result[[7]]
                WilcoxTest.results$P.value[r]<-result[[3]]
        }
        write.csv(WilcoxTest.results,paste0("Output/SummaryStats/Wilcox_Maj_TsMutFreq_",fname,".csv"))
}

hist(G_tv1_syn_CpG)
hist(G_tv1_syn_NonCpG)
boxplot(log(G_tv1_syn_CpG),log(G_tv1_syn_NonCpG))
       
