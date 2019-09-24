library(reshape)
library(tidyverse)
library(zoo)
library(purrr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(sfsmisc)
source("Rscripts/baseRscript.R")


#####

HCVFiles_overview3<-list.files("Output1A/Overview3/",pattern="overview3.csv")
s<-length(HCVFiles_overview3)


##### read the saved summary mutation frequency files
TsMutFreq<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F,row.names = 1)
mean(TsMutFreq$mean) #0.004804249

TsMutFreq<-TsMutFreq[TsMutFreq$pos>=342,]
mean(TsMutFreq$mean) #0.004828831

TsMutFreq2<-TsMutFreq
TsMutFreq3<-TsMutFreq2

#remove the low viral load samples
info<-read.csv("Data/CoInfectionstudy_HCV_Genomes.csv",stringsAsFactors = F)

low_ids<-info$Enumber[info$Sequence.Genotype=="1A"&info$HCV_load<100000]
low_ids<-low_ids[!is.na(low_ids)]
low_ids2<- paste0(low_ids, ".")
low_ids2<-strtrim(low_ids2,7) #18

Trans<-TsMutFreq2[,- which(names(TsMutFreq2) %in% low_ids2)]  # 183 samples (12 removed)
LowViralCount<-TsMutFreq2[, which(names(TsMutFreq2) %in% low_ids2)]   #12 samples


# compare the low viral count samples vs. those that are not:
mean(rowMeans(Trans[2:184],na.rm=T)) #0.004845944
LowViralCount$means<-rowMeans(LowViralCount,na.rm=T)
mean(LowViralCount$means, na.rm=T)  #0.004535194

#Difference between all 195 samples vs. 183 samples  -> very small
mean(TsMutFreq$mean)- mean(rowMeans(Trans[2:184],na.rm=T)) ## -1.452781e-05

############################################################################################

# 1.Filter out Mut Freq < 0.001  -> replaced with NA
mf<-0.001
TsMutFreq2[, c(2:(s+1))][TsMutFreq2[, 2:(s+1)] <mf]<-NA
#Tv1.MutFreq2[Tv1.MutFreq2<mf]<-NA
#Tv2.MutFreq2[Tv2.MutFreq2<mf]<-NA
#Tvs.MutFreq2[Tvs.MutFreq2<mf]<-NA
#AllMutFreq2[AllMutFreq2<mf]<-NA

#Filter out Mut Freq < 0.001  -> replaced with 0
mf<-0.001
TsMutFreq3[, c(2:(s+1))][TsMutFreq3[, c(2:(s+1))]<mf]<-0
#Tv1.MutFreq3[Tv1.MutFreq3<mf]<-0
#Tv2.MutFreq3[Tv2.MutFreq3<mf]<-0
#Tvs.MutFreq3[Tvs.MutFreq3<mf]<-0
#AllMutFreq3[AllMutFreq3<mf]<-0


### Look at Transition mutations
Ts_NA<-TsMutFreq2
Ts_zero<-TsMutFreq3

#mean
Ts_NA$mean<-rowMeans(Ts_NA[2:(s+1)], na.rm=T)
mean(Ts_NA$mean) #0.005236761

Ts_zero$mean<-rowMeans(Ts_zero[2:(s+1)],na.rm=T)
mean(Ts_zero$mean) #0.004751198

#save coding regions only
write.csv(TsMutFreq, "Output1A/Q35Compare/Summary_Ts.Q35.csv")
write.csv(Ts_NA, "Output1A/Q35Compare/Summary_Ts_NA.Q35.csv")
write.csv(Ts_zero, "Output1A/Q35Compare/Summary_Ts_zero.Q35.csv")


#########################################################################################
##Plot across the genome

endnuc<-TsMutFreq2$pos[nrow(TsMutFreq2)]
SNPFreq<-Ts_NA

n<-data.frame("pos"=c(342:endnuc))
SNPFreqs<-merge(n,SNPFreq,by="pos",all.x=T)
pdf(paste0("Output1A/Q35Compare/Ts.mutfreq.NAreplaced.pdf"),width=15,height=7.5)
plot(mean~pos, data=SNPFreqs,t="n",log='y',yaxt='n',xlab='Genome position (H77)',ylab="Average transition mutation frequency",
     main="Transition muttaion frequency ",ylim=c(0.0001,0.1),xlim=c(340,8500))
eaxis(side = 2, at = 10^((-1):(-(4))), cex=2)
for(i in 2:4){abline(h = 1:10 * 10^(-i), col = "gray80")}
points(mean~pos, data=SNPFreqs,pch=20,col="#EE667766",cex=0.5)

#add rolling average
roll100<-rollmean(SNPFreqs$mean, k=100, na.rm=T)
SNPFreqs$roll100<-c(rep(NA, times=50),roll100,c(rep(NA, times=49)))
lines(roll100~pos,data=SNPFreqs, col="#AA3377")

dev.off()

###
SNPFreq<-Ts_zero

n<-data.frame("pos"=c(1:endnuc))
SNPFreqs<-merge(n,SNPFreq,by="pos",all.x=T)
pdf(paste0("Output1A/Q35Compare/Ts.mutfreq.Zero.replaced.pdf"),width=15,height=7.5)
plot(mean~pos, data=SNPFreqs,t="n",log='y',yaxt='n',xlab='Genome position (H77)',ylab="Average transition mutation frequency",
     main="Transition muttaion frequency ",ylim=c(0.0001,0.1),xlim=c(340,8500))
eaxis(side = 2, at = 10^((-1):(-(4))), cex=2)
for(i in 2:4){abline(h = 1:10 * 10^(-i), col = "gray80")}
points(mean~pos, data=SNPFreqs,pch=20,col="#EE667766",cex=0.5)

#add rolling average
roll100<-rollmean(SNPFreqs$mean, k=100, na.rm=T)
SNPFreqs$roll100<-c(rep(NA, times=50),roll100,c(rep(NA, times=49)))
lines(roll100~pos,data=SNPFreqs, col="#AA3377")

dev.off()

#######  no filtering of low mut freq  #######
SNPFreq<-TsMutFreq
SNPFreqs<-merge(n,SNPFreq,by="pos",all.x=T)
pdf(paste0("Output1A/Q35Compare/Ts.mutfreq.NoFilter.pdf"),width=15,height=7.5)
plot(mean~pos, data=SNPFreqs,t="n",log='y',yaxt='n',xlab='Genome position (H77)',ylab="Average transition mutation frequency",
     main="Transition muttaion frequency ",ylim=c(0.0001,0.1),xlim=c(340,8500))
eaxis(side = 2, at = 10^((-1):(-(4))), cex=2)
for(i in 2:4){abline(h = 1:10 * 10^(-i), col = "gray80")}
points(mean~pos, data=SNPFreqs,pch=20,col="#EE667766",cex=0.5)

#add rolling average
roll100<-rollmean(SNPFreqs$mean, k=100, na.rm=T)
SNPFreqs$roll100<-c(rep(NA, times=50),roll100,c(rep(NA, times=49)))
lines(roll100~pos,data=SNPFreqs, col="#AA3377")

dev.off()



###################################################################
## compare 1.Q35 reads>1000, 2.mf<0.001=NA, 3.mf<0.001=0 4. Q30 reads>1000

TsMutQ30<-read.csv("Output1A/Q35Compare/Summary_TsMutFreq_Q30.csv",stringsAsFactors = F)
TsMutQ30<-TsMutQ30[,-1]
mean(TsMutQ30$mean, na.rm=T)  #0.005426269


SumTab<-data.frame("data"=c("Ts.Q35","Ts_NA","Ts_zero","Ts.Q30"))

t<-list()


t[[1]]<-TsMutFreq
t[[2]]<-Ts_NA
t[[3]]<-Ts_zero
t[[4]]<-TsMutQ30
names(t)<-c("Ts.Q35","Ts_NA","Ts_zero","Ts.Q30")

for (i in 1:4){
        Ts<-t[[i]]
        SumTab$mean[i]<-mean(Ts$mean, na.rm=T)
        SumTab$Syn[i]<-mean(Ts$mean[Ts$Type=="syn"],na.rm=T) 
        SumTab$Nonsyn[i]<-mean(Ts$mean[Ts$Type=="nonsyn"],na.rm=T)
        SumTab$Stop[i]<-mean(Ts$mean[Ts$Type=="stop"],na.rm=T)
        
        T2<-Ts[Ts$ref=="a"|Ts$ref=="t",]
        SumTab$CpGmaking_Syn[i]<-mean(Ts$mean[Ts$Type=="syn"&Ts$makesCpG==1],na.rm=T)
        SumTab$NoCpGmaking_Syn[i]<-mean(T2$mean[T2$Type=="syn"&T2$makesCpG==0],na.rm=T)
        r3<-wilcox.test(T2$mean[T2$Type=="syn"&T2$makesCpG==0], T2$mean[T2$Type=="syn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 
        SumTab$P_value.Wilcoxon.Test_syn[i]<-r3[[3]]
        SumTab$CpGmaking_Nonsyn[i]<-mean(Ts$mean[Ts$Type=="nonsyn"&Ts$makesCpG==1],na.rm=T)
        SumTab$NoCpGmaking_Nonsyn[i]<-mean(T2$mean[T2$Type=="nonsyn"&T2$makesCpG==0],na.rm=T)
        r4<-wilcox.test(T2$mean[T2$Type=="nonsyn"&T2$makesCpG==0], T2$mean[T2$Type=="nonsyn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 
        SumTab$P_value.Wilcoxon.Test_nonsyn[i]<-r4[[3]]
}

write.csv(SumTab, "Output1A/Q35Compare/Comparison_summary.Ts.csv")

## Test whether 3 are different 
SumTab<-data.frame("data"=c("Ts.Q35","Ts_NA","Ts_zero"))

t<-list()
t[[1]]<-TsMutFreq
t[[2]]<-Ts_NA
t[[3]]<-Ts_zero
names(t)<-c("Ts.Q35","Ts_NA","Ts_zero")

T1<-t[[1]]
T2<-t[[2]]
T3<-t[[3]]
wilcox.test(T1$mean, T2$mean, alternative = "less", paired = FALSE) 
#W = 28830000, p-value < 2.2e-16
wilcox.test(T1$mean, T3$mean, alternative = "greater", paired = FALSE) 
#W = 32347000, p-value = 0.009262
wilcox.test(T2$mean, T3$mean, alternative = "greater", paired = FALSE) 
#W = 35046000, p-value < 2.2e-16



for (i in 1:length(t)){
        Ts<-t[[i]]
        SumTab$mean[i]<-mean(Ts$mean, na.rm=T)
        SumTab$Syn[i]<-mean(Ts$mean[Ts$Type=="syn"],na.rm=T) 
        SumTab$Nonsyn[i]<-mean(Ts$mean[Ts$Type=="nonsyn"],na.rm=T)
        SumTab$Stop[i]<-mean(Ts$mean[Ts$Type=="stop"],na.rm=T)
        
        r1<-wilcox.test(Ts$mean, Ts2$mean, alternative = "greater", paired = FALSE) 
        
        
        T2<-Ts[Ts$ref=="a"|Ts$ref=="t",]
        SumTab$CpGmaking_Syn[i]<-mean(Ts$mean[Ts$Type=="syn"&Ts$makesCpG==1],na.rm=T)
        SumTab$NoCpGmaking_Syn[i]<-mean(T2$mean[T2$Type=="syn"&T2$makesCpG==0],na.rm=T)
        r3<-wilcox.test(T2$mean[T2$Type=="syn"&T2$makesCpG==0], T2$mean[T2$Type=="syn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 
        SumTab$P_value.Wilcoxon.Test_syn[i]<-r3[[3]]
        SumTab$CpGmaking_Nonsyn[i]<-mean(Ts$mean[Ts$Type=="nonsyn"&Ts$makesCpG==1],na.rm=T)
        SumTab$NoCpGmaking_Nonsyn[i]<-mean(T2$mean[T2$Type=="nonsyn"&T2$makesCpG==0],na.rm=T)
        r4<-wilcox.test(T2$mean[T2$Type=="nonsyn"&T2$makesCpG==0], T2$mean[T2$Type=="nonsyn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 
        SumTab$P_value.Wilcoxon.Test_nonsyn[i]<-r4[[3]]
}

write.csv(SumTab, "Output1A/Q35Compare/Comparison_summary.Ts.csv")
