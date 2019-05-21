library(reshape)
library(tidyverse)
library(zoo)
library(purrr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(sfsmisc)
source("Rscripts/baseRscript.R")


#Load the overview files (summarized using the ref (H77))
#####
if (FALSE){
HCVFiles_overview3<-list.files("Output/Overview_filtered/",pattern="overview3.csv")
FilteredOverview<-list()
for (i in 1:length(HCVFiles_overview3)){ 
     overviews<-read.csv(paste0("Output/Overview_filtered/",HCVFiles_overview3[i]),stringsAsFactors=FALSE)
     overviews<-overviews[,-1]
     FilteredOverview[[i]]<-overviews
     names(FilteredOverview)[i]<-substr(paste(HCVFiles_overview3[i]),start=1,stop=7)
     }
}
######################################################
s<-length(FilteredOverview)
##### read the saved summary mutation frequency files
TsMutFreq<-read.csv("Output/MutFreq/Filtered_Transition_MutFreqSummary.csv",stringsAsFactors = F)
Tv1.MutFreq<-read.csv("Output/MutFreq/Filtered_Transv1_MutFreqSummary.csv",stringsAsFactors = F)
Tv2.MutFreq<-read.csv("Output/MutFreq/Filtered_Transv2_MutFreqSummary.csv",stringsAsFactors = F)
Tvs.MutFreq<-read.csv("Output/MutFreq/Filtered_Transv_MutFreqSummary.csv",stringsAsFactors = F)
AllMutFreq<-read.csv("Output/MutFreq/Filtered_All_MutFreqSummary.csv",stringsAsFactors = F)

TsMutFreq<-TsMutFreq[,-1]
Tv1.MutFreq<-Tv1.MutFreq[,-1]
Tv2.MutFreq<-Tv2.MutFreq[,-1]
Tvs.MutFreq<-Tvs.MutFreq[,-1]
AllMutFreq<-AllMutFreq[,-1]

#remove the sites before the coding region
TsMutFreq2<-TsMutFreq[TsMutFreq$pos>=342,]
Tv1.MutFreq2<-Tv1.MutFreq[Tv1.MutFreq$pos>=342,]
Tv2.MutFreq2<-Tv2.MutFreq[Tv2.MutFreq$pos>=342,]
Tvs.MutFreq2<-Tvs.MutFreq[Tvs.MutFreq$pos>=342,]
AllMutFreq2<-AllMutFreq[AllMutFreq$pos>=342,]

TsMutFreq3<-TsMutFreq2
Tv1.MutFreq3<-Tv1.MutFreq2
Tv2.MutFreq3<-Tv2.MutFreq2
Tvs.MutFreq3<-Tvs.MutFreq2
AllMutFreq3<-AllMutFreq2


#remove the low viral load samples
info<-read.csv("Data/CoInfectionstudy_HCV_Genomes.csv",stringsAsFactors = F)

low_ids<-info$Enumber[info$Sequence.Genotype=="1A"&info$HCV_load<100000]
low_ids<-low_ids[!is.na(low_ids)]
low_ids2<- paste0(low_ids, ".")
low_ids2<-strtrim(low_ids2,7) #18

Trans<-TsMutFreq2[,- which(names(TsMutFreq2) %in% low_ids2)]  # 184 samples (12 removed)
LowViralCount<-TsMutFreq2[, which(names(TsMutFreq2) %in% low_ids2)]   #12 samples


# compare the low viral count samples vs. those that are not:
mean(rowMeans(Trans[2:184],na.rm=T)) #0.005441847
LowViralCount$means<-rowMeans(LowViralCount,na.rm=T)
mean(LowViralCount$means, na.rm=T)  #0.005171185
 ### Difference is very small.
 ### keep all of them for now.


# 1.Filter out Mut Freq < 0.001  -> replaced with NA
mf<-0.001
TsMutFreq2[TsMutFreq2<mf]<-NA
Tv1.MutFreq2[Tv1.MutFreq2<mf]<-NA
Tv2.MutFreq2[Tv2.MutFreq2<mf]<-NA
Tvs.MutFreq2[Tvs.MutFreq2<mf]<-NA
AllMutFreq2[AllMutFreq2<mf]<-NA

#Filter out Mut Freq < 0.001  -> replaced with 0
mf<-0.001
TsMutFreq3[TsMutFreq3<mf]<-0
Tv1.MutFreq3[Tv1.MutFreq3<mf]<-0
Tv2.MutFreq3[Tv2.MutFreq3<mf]<-0
Tvs.MutFreq3[Tvs.MutFreq3<mf]<-0
AllMutFreq3[AllMutFreq3<mf]<-0


### Look at Transition mutations
dat<-FilteredOverview[[3]]
muttypes<-dat[,c("pos","ref","Type","Type.tv1","Type.tv2","WTAA","MUTAA","TVS1_AA","TVS2_AA","makesCpG","makesCpG.tv1","makesCpG.tv2","bigAAChange","bigAAChange.tv1","bigAAChange.tv2")]
Ts_NA<-merge(TsMutFreq2,muttypes,by="pos")
Ts_zero<-merge(TsMutFreq3,muttypes, by="pos")
TsMutFreq<-merge(TsMutFreq,muttypes, by="pos")

#calculate the mean
TsMutFreq<-TsMutFreq[TsMutFreq$pos>=342,]
TsMutFreq$mean<-rowMeans(TsMutFreq[2:(s+1)], na.rm=T)
mean(TsMutFreq$mean) #0.005426269

Ts_NA$mean<-rowMeans(Ts_NA[2:(s+1)], na.rm=T)
mean(Ts_NA$mean) #0.005679108

Ts_zero$mean<-rowMeans(Ts_zero[2:(s+1)],na.rm=T)
mean(Ts_zero$mean) #0.005376377


write.csv(TsMutFreq, "Output/Mut.freq.filtered/Summary_TsMutFreq.csv")
write.csv(Ts_NA, "Output/Mut.freq.filtered/Summary_Ts_NA.csv")
write.csv(Ts_zero, "Output/Mut.freq.filtered/Summary_Ts_zero.csv")


#########################################################################################
##Plot across the genome

endnuc<-TsMutFreq2$pos[nrow(TsMutFreq2)]
SNPFreq<-Ts_NA

n<-data.frame("pos"=c(342:endnuc))
SNPFreqs<-merge(n,SNPFreq,by="pos",all.x=T)
pdf(paste0("Output/SummaryFig.Filtered/Ts.MutFreq.NAreplaced.pdf"),width=15,height=7.5)
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
pdf(paste0("Output/SummaryFig.Filtered/Ts.MutFreq.Zero.replaced.pdf"),width=15,height=7.5)
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
pdf(paste0("Output/SummaryFig.Filtered/Ts.MutFreq.NoFilter_over1000.pdf"),width=15,height=7.5)
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
## compare 1.reads>100, 2. reads>1000, 3.mf<0.001=NA, 4.mf<0.001=0

TsMut<-read.csv("Output/MutFreq/Summary_for_Reads_100Removed/Filtered_Transition_MutFreqSummary.csv",stringsAsFactors = F)
TsMut<-TsMut[,-1]
TsMut<-TsMut[TsMut$pos>=342,]
TsMut$mean<-rowMeans(TsMut[2:196], na.rm=T)
TsMut<-merge(TsMut,muttypes, by="pos")
SumTab<-data.frame("data"=c("Ts_over100","Ts_over1000","Ts_NA","Ts_zero"))

t<-list()

t[[1]]<-TsMut
t[[2]]<-TsMutFreq
t[[3]]<-Ts_NA
t[[4]]<-Ts_zero
names(t)<-c("Ts_over100","Ts_over1000","Ts_NA","Ts_zero")

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

write.csv(SumTab, "Output/ReadsVs.MutFreq/Comparison_summary.csv")

# 1) NA replaced Transition
Ts<-Ts_NA
## Compare Syn vs Nonsyn, Nonsyn vs. Stop
r1<-wilcox.test(Ts$mean[Ts$Type=="syn"], Ts$mean[Ts$Type=="nonsyn"], alternative = "greater", paired = FALSE) #W = 11795000, p-value < 2.2e-16
r2<-wilcox.test(Ts$mean[Ts$Type=="nonsyn"], Ts$mean[Ts$Type=="stop"], alternative = "greater", paired = FALSE) 

# syn vs nonsyn
r1[[3]]  #0
r2[[3]]  #6.72538e-37

## Compare CpG making vs. Non-CpGmaking 
T2<-Ts[Ts$ref=="a"|Ts$ref=="t",]
r3<-wilcox.test(T2$mean[T2$Type=="syn"&T2$makesCpG==0], T2$mean[T2$Type=="syn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 
r4<-wilcox.test(T2$mean[T2$Type=="nonsyn"&T2$makesCpG==0], T2$mean[T2$Type=="nonsyn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 
#syn -CpG
r3[[3]] #0.01079874
#nonsyn CpG
r4[[3]] #5.896308e-10


##################
# 1.2) run Wilcoxin Test on means based on mutation types

WilcoxTest.results.nt<-data.frame(matrix(ncol=3,nrow=6))
colnames(WilcoxTest.results.nt)<-c("nt","test","P.value")

# 1) transition Mutations, using 'mean'
dat<-Ts_NA
ty<-which(colnames(dat)=="Type");fname="Transition"
m<-data.frame()
se<-data.frame()
m_CpG<-data.frame()
se_CpG<-data.frame()
m_nonCpG<-data.frame()
se_nonCpG<-data.frame()

for (typeofsite in c("syn", "nonsyn","stop")){
        for (wtnt in c("a", "t", "c", "g")){
                mutrate<- dat$mean[dat[,ty]==typeofsite & dat$ref==wtnt]
                m[typeofsite,wtnt]<-mean(mutrate[!is.na(mutrate)])
                se[typeofsite,wtnt]<-std.error(mutrate[!is.na(mutrate)])
                
                m_NonCpG<-dat$mean[dat$Type==typeofsite & dat$ref==wtnt & dat$makesCpG==0]
                m_nonCpG[typeofsite,wtnt]<-mean(m_NonCpG[!is.na(m_NonCpG)])
                se_nonCpG[typeofsite,wtnt]<-std.error(m_NonCpG[!is.na(m_NonCpG)])
                
                mu_CpG<-dat$mean[dat$Type==typeofsite & dat$ref==wtnt & dat$makesCpG==1]
                m_CpG[typeofsite,wtnt]<-mean(mu_CpG[!is.na(mu_CpG)])
                se_CpG[typeofsite,wtnt]<-std.error(mu_CpG[!is.na(mu_CpG)])
                
                vectorname<-paste0(typeofsite,"_",wtnt)
                assign(vectorname, mutrate)
                vname1<<-paste0(typeofsite,"_",wtnt,"_noncpg")
                assign(vname1, m_NonCpG)
                vname2<<-paste0(typeofsite,"_",wtnt,"_cpg")
                assign(vname2, mu_CpG)
        }
}
#rownames(m)<-c("syn","nonsyn","stop")
rownames(se)<-c("syn_se","nonsyn_se","stop_se")
rownames(m_nonCpG)<-c("syn_noncpg","nonsyn_noncpg","stop_noncpg")
rownames(se_nonCpG)<-c("syn_noncpg_se","nonsyn_noncpg_se","stop_noncpg_se")
rownames(m_CpG)<-c("syn_cpg","nonsyn_cpg","stop_cpg")
rownames(se_CpG)<-c("syn_cpg_se","nonsyn_cpg_se","stop_cpg_se")

MFbyType<-rbind(m,se,m_nonCpG,se_nonCpG,m_CpG,se_CpG)
MFbyType2<-t(MFbyType)
MFbyType2<-data.frame(MFbyType2)
write.csv(MFbyType2,"Output/SummaryStats/TsMF_byNt_byType_mean_NAreplaced.csv")

#run Wilcoxin Test  
for (i in c("a","t","c","g")) {
        if (i=="a"|i=="t"){
                syncpg<-get(paste0("syn_",i,"_cpg"))
                synnoncpg<-get(paste0("syn_",i,"_noncpg"))
                nonsyncpg<-get(paste0("nonsyn_",i,"_cpg"))
                nonsynnoncpg<-get(paste0("nonsyn_",i,"_noncpg"))
                if (i=="a"){
                        result1<-wilcox.test(syncpg, synnoncpg, alternative = "less", paired = FALSE) 
                        result2<-wilcox.test(nonsyncpg,nonsynnoncpg,alternative = "less", paired = FALSE)
                        for (r in 1:2){
                                result<-get(paste0('result',r))
                                WilcoxTest.results.nt$nt[r]<-i
                                WilcoxTest.results.nt$test[r]<-result[[7]]
                                WilcoxTest.results.nt$P.value[r]<-result[[3]]}
                }
                if (i=="t"){
                        result3<-wilcox.test(syncpg, synnoncpg, alternative = "less", paired = FALSE) 
                        result4<-wilcox.test(nonsyncpg,nonsynnoncpg,alternative = "less", paired = FALSE)
                        for (r in 3:4){
                                result<-get(paste0('result',r))
                                WilcoxTest.results.nt$nt[r]<-i
                                WilcoxTest.results.nt$test[r]<-result[[7]]
                                WilcoxTest.results.nt$P.value[r]<-result[[3]]}
                }
        }
        else {  synnoncpg<-get(paste0("syn_",i,"_noncpg"))
        nonsynnoncpg<-get(paste0("nonsyn_",i,"_noncpg"))
        if (i =="c") {
                result5<-wilcox.test(synnoncpg,nonsynnoncpg, alternative = "greater", paired = FALSE) 
                WilcoxTest.results.nt$nt[5]<-i
                WilcoxTest.results.nt$test[5]<-result5[[7]]
                WilcoxTest.results.nt$P.value[5]<-result5[[3]]     }           
        if (i =="g") {
                result6<-wilcox.test(synnoncpg,nonsynnoncpg, alternative = "greater", paired = FALSE) 
                WilcoxTest.results.nt$nt[6]<-i
                WilcoxTest.results.nt$test[6]<-result6[[7]]
                WilcoxTest.results.nt$P.value[6]<-result6[[3]]    }
        }
} 
write.csv(WilcoxTest.results.nt,paste0("Output/SummaryStats/WilcoxTestResults_Ts_byNt_mean_NAreplaced.csv"))

########
# 2) zero replaced
Ts<-Ts_zero

mean(Ts$mean[Ts$Type=="syn"],na.rm=T) #0.00911283  vs. #0.009074511 (replaced <0.001)
mean(Ts$mean[Ts$Type=="nonsyn"],na.rm=T) #0.003973061  vs # 0.003928531 (replaced <0.001)
mean(Ts$mean[Ts$Type=="stop"],na.rm=T) #0.002562846  vs. 0.002501161 (replaced <0.001)

## Compare Syn vs Nonsyn, Nonsyn vs. Stop
r1<-wilcox.test(Ts$mean[Ts$Type=="syn"], Ts$mean[Ts$Type=="nonsyn"], alternative = "greater", paired = FALSE) #W = 11795000, p-value < 2.2e-16
r2<-wilcox.test(Ts$mean[Ts$Type=="nonsyn"], Ts$mean[Ts$Type=="stop"], alternative = "greater", paired = FALSE) 

# syn vs nonsyn
r1[[3]]  #0
r2[[3]]  #4.004883e-30

## Compare CpG making vs. Non-CpGmaking 
T2<-Ts[Ts$ref=="a"|Ts$ref=="t",]
r3<-wilcox.test(T2$mean[T2$Type=="syn"&T2$makesCpG==0], T2$mean[T2$Type=="syn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 
r4<-wilcox.test(T2$mean[T2$Type=="nonsyn"&T2$makesCpG==0], T2$mean[T2$Type=="nonsyn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 

#syn -CpG
r3[[3]] #0.007407342
#nonsyn CpG
r4[[3]] #3.689598e-09

##################
# 1.2) run Wilcoxin Test on means based on mutation types

WilcoxTest.results.nt<-data.frame(matrix(ncol=3,nrow=6))
colnames(WilcoxTest.results.nt)<-c("nt","test","P.value")

# 1) transition Mutations, using 'mean'
dat<-Ts_zero
ty<-which(colnames(dat)=="Type");fname="Transition"
m<-data.frame()
se<-data.frame()
m_CpG<-data.frame()
se_CpG<-data.frame()
m_nonCpG<-data.frame()
se_nonCpG<-data.frame()

for (typeofsite in c("syn", "nonsyn","stop")){
        for (wtnt in c("a", "t", "c", "g")){
                mutrate<- dat$mean[dat[,ty]==typeofsite & dat$ref==wtnt]
                m[typeofsite,wtnt]<-mean(mutrate[!is.na(mutrate)])
                se[typeofsite,wtnt]<-std.error(mutrate[!is.na(mutrate)])
                
                m_NonCpG<-dat$mean[dat$Type==typeofsite & dat$ref==wtnt & dat$makesCpG==0]
                m_nonCpG[typeofsite,wtnt]<-mean(m_NonCpG[!is.na(m_NonCpG)])
                se_nonCpG[typeofsite,wtnt]<-std.error(m_NonCpG[!is.na(m_NonCpG)])
                
                mu_CpG<-dat$mean[dat$Type==typeofsite & dat$ref==wtnt & dat$makesCpG==1]
                m_CpG[typeofsite,wtnt]<-mean(mu_CpG[!is.na(mu_CpG)])
                se_CpG[typeofsite,wtnt]<-std.error(mu_CpG[!is.na(mu_CpG)])
                
                vectorname<-paste0(typeofsite,"_",wtnt)
                assign(vectorname, mutrate)
                vname1<<-paste0(typeofsite,"_",wtnt,"_noncpg")
                assign(vname1, m_NonCpG)
                vname2<<-paste0(typeofsite,"_",wtnt,"_cpg")
                assign(vname2, mu_CpG)
        }
}
#rownames(m)<-c("syn","nonsyn","stop")
rownames(se)<-c("syn_se","nonsyn_se","stop_se")
rownames(m_nonCpG)<-c("syn_noncpg","nonsyn_noncpg","stop_noncpg")
rownames(se_nonCpG)<-c("syn_noncpg_se","nonsyn_noncpg_se","stop_noncpg_se")
rownames(m_CpG)<-c("syn_cpg","nonsyn_cpg","stop_cpg")
rownames(se_CpG)<-c("syn_cpg_se","nonsyn_cpg_se","stop_cpg_se")

MFbyType<-rbind(m,se,m_nonCpG,se_nonCpG,m_CpG,se_CpG)
MFbyType2<-t(MFbyType)
MFbyType2<-data.frame(MFbyType2)
write.csv(MFbyType2,"Output/SummaryStats/TsMF_byNt_byType_mean_ZEROreplaced.csv")

#run Wilcoxin Test  
for (i in c("a","t","c","g")) {
        if (i=="a"|i=="t"){
                syncpg<-get(paste0("syn_",i,"_cpg"))
                synnoncpg<-get(paste0("syn_",i,"_noncpg"))
                nonsyncpg<-get(paste0("nonsyn_",i,"_cpg"))
                nonsynnoncpg<-get(paste0("nonsyn_",i,"_noncpg"))
                if (i=="a"){
                        result1<-wilcox.test(syncpg, synnoncpg, alternative = "less", paired = FALSE) 
                        result2<-wilcox.test(nonsyncpg,nonsynnoncpg,alternative = "less", paired = FALSE)
                        for (r in 1:2){
                                result<-get(paste0('result',r))
                                WilcoxTest.results.nt$nt[r]<-i
                                WilcoxTest.results.nt$test[r]<-result[[7]]
                                WilcoxTest.results.nt$P.value[r]<-result[[3]]}
                }
                if (i=="t"){
                        result3<-wilcox.test(syncpg, synnoncpg, alternative = "less", paired = FALSE) 
                        result4<-wilcox.test(nonsyncpg,nonsynnoncpg,alternative = "less", paired = FALSE)
                        for (r in 3:4){
                                result<-get(paste0('result',r))
                                WilcoxTest.results.nt$nt[r]<-i
                                WilcoxTest.results.nt$test[r]<-result[[7]]
                                WilcoxTest.results.nt$P.value[r]<-result[[3]]}
                }
        }
        else {  synnoncpg<-get(paste0("syn_",i,"_noncpg"))
        nonsynnoncpg<-get(paste0("nonsyn_",i,"_noncpg"))
        if (i =="c") {
                result5<-wilcox.test(synnoncpg,nonsynnoncpg, alternative = "greater", paired = FALSE) 
                WilcoxTest.results.nt$nt[5]<-i
                WilcoxTest.results.nt$test[5]<-result5[[7]]
                WilcoxTest.results.nt$P.value[5]<-result5[[3]]     }           
        if (i =="g") {
                result6<-wilcox.test(synnoncpg,nonsynnoncpg, alternative = "greater", paired = FALSE) 
                WilcoxTest.results.nt$nt[6]<-i
                WilcoxTest.results.nt$test[6]<-result6[[7]]
                WilcoxTest.results.nt$P.value[6]<-result6[[3]]    }
        }
} 
write.csv(WilcoxTest.results.nt,paste0("Output/SummaryStats/WilcoxTestResults_Ts_byNt_mean_ZEROreplaced.csv"))

