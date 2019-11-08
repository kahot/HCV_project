library(tidyverse)
library(zoo)
library(purrr)

source("Rscripts/baseRscript.R")


HCVFiles3<-list.files("Output1A/Overview3/",pattern="overview3.csv")

#  #### read the filtered MF files from 10.1.2

TsMutFreq  <-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F,row.names=1)
Tv1.MutFreq<-read.csv("Output1A/MutFreq.filtered/Filtered.Tv1.MutFreq.Q35.csv",stringsAsFactors = F,row.names=1)
Tv2.MutFreq<-read.csv("Output1A/MutFreq.filtered/Filtered.Tv2.MutFreq.Q35.csv",stringsAsFactors = F,row.names=1)
Tvs.MutFreq<-read.csv("Output1A/MutFreq.filtered/Filtered.Tvs.MutFreq.Q35.csv",stringsAsFactors = F,row.names=1)
AllMutFreq <-read.csv("Output1A/MutFreq.filtered/Filtered.AllMutFreq.Q35.csv", stringsAsFactors = F,row.names=1)

mf.files<-list()
mf.files[[1]]<-TsMutFreq
mf.files[[2]]<-Tv1.MutFreq
mf.files[[3]]<-Tv2.MutFreq
mf.files[[4]]<-Tvs.MutFreq
mf.files[[5]]<-AllMutFreq
names(mf.files)[1]<-"TransMutFreq"
names(mf.files)[2]<-"Tv1.MutFreq"
names(mf.files)[3]<-"Tv2.MutFreq"
names(mf.files)[4]<-"Tvs.MutFreq"
names(mf.files)[5]<-"AllMutFreq"

s<-length(HCVFiles3)

#### 1) Transition Mutations (mean)
# 1.1) Summary
Ts<-mf.files[[1]]
Ts<-Ts[Ts$pos>=342, ]
mean(Ts$mean[Ts$Type=="syn"]) #0.008093379
mean(Ts$mean[Ts$Type=="nonsyn"]) #0.003348242
mean(Ts$mean[Ts$Type=="stop"]) #0.002022381

table(Ts$Type)
#nonsyn   stop    syn 
#  5193    220   2545  

std.error(Ts$mean[Ts$Type=="syn"]) #9.066418e-05
std.error(Ts$mean[Ts$Type=="nonsyn"]) # 2.768315e-05
std.error(Ts$mean[Ts$Type=="stop"]) # 9.204526e-05

r1<-wilcox.test(Ts$mean[Ts$Type=="syn"], Ts$mean[Ts$Type=="nonsyn"], alternative = "greater", paired = FALSE) 
r2<-wilcox.test(Ts$mean[Ts$Type=="nonsyn"], Ts$mean[Ts$Type=="stop"], alternative = "greater", paired = FALSE) 
r1[[3]]  #P=0
r2[[3]]  #P= 1.933446e-41

<<<<<<< HEAD
=======

>>>>>>> Updated analysis scrits
#CpG creating vs Non-CpG creating
T2<-Ts[Ts$ref=="a"|Ts$ref=="t",]
r3<-wilcox.test(T2$mean[T2$Type=="syn"&T2$makesCpG==0], T2$mean[T2$Type=="syn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 
r4<-wilcox.test(T2$mean[T2$Type=="nonsyn"&T2$makesCpG==0], T2$mean[T2$Type=="nonsyn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 

r3[[3]]  #P=0.005231179
r4[[3]]  #P=3.893716e-12


##################
# 1.2) run Wilcoxin Test on means based on mutation types

WilcoxTest.results.nt<-data.frame(matrix(ncol=3,nrow=6))
colnames(WilcoxTest.results.nt)<-c("nt","test","P.value")

# 1) transition Mutations, using 'mean'
dat<-mf_files[[1]]
#dat<-TransMutFreq_filtered
dat<-dat[dat$pos>=342, ]
which(is.na(dat$freq.Ts))
dat<-dat[1:8236, ]#7895
ty<-which(colnames(dat)=="Type");fname="Transition"
m<-data.frame()
se<-data.frame()
m_CpG<-data.frame()
se_CpG<-data.frame()
m_nonCpG<-data.frame()
se_nonCpG<-data.frame()

table(dat$ref)
#a    c    g    t 
#1556 2431 2310 1661  


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
write.csv(MFbyType2,"Output1A/SummaryStats/TransitionMF_byNt_byType_mean.csv")

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
write.csv(WilcoxTest.results.nt,paste0("Output1A/SummaryStats/WilcoxTestResults_Ts_eachNT_mean.csv"))

al<-mf.files[[5]]
mean(al$mean) #0.005771959

tvs<-mf.files[[4]]
tvs<-tvs[tvs$pos>=342,]
mean(tvs$mean) #0.0009650088

# 2) Transversion with mean ###
## 2.1) Summary Stats
Tv1<-mf.files[[2]]
Tv1<-Tv1[Tv1$pos>=342, ]
Tv2<-mf.files[[3]]
Tv2<-Tv2[Tv2$pos>=342, ]

All<-c(Tv1$mean,Tv2$mean) #0.0004825044
Syn<-c(Tv1$mean[Tv1$Type.tv1=="syn"],Tv2$mean[Tv2$Type.tv2=="syn"])
Nonsyn <- c(Tv1$mean[Tv1$Type.tv1=="nonsyn"],Tv2$mean[Tv2$Type.tv2=="nonsyn"])
Stop <- c(Tv1$mean[Tv1$Type.tv1=="stop"],Tv2$mean[Tv2$Type.tv2=="stop"])

mean(All) #0.0004825044
mean(Syn) #0.0007535948
mean(Nonsyn) #0.0004099817
mean(Stop) #0.0005612278

c(length(Nonsyn),length(Stop),length(Syn))
#nonsyn   stop  syn 
# 12184   625  3107 

std.error(Syn) # 1.508687e-05
std.error(Nonsyn) # 4.921358e-06
std.error(Stop) # 1.833373e-05


r1<-wilcox.test(Syn, Nonsyn, alternative = "greater", paired = FALSE) 
r2<-wilcox.test(Nonsyn, Stop, alternative = "greater", paired = FALSE) 
wilcox.test(Nonsyn, Stop, alternative = "less", paired = FALSE) 

r1[[3]] #2.811567e-175
r2[[3]] #1


Syncpg<-c(Tv1$mean[Tv1$Type.tv1=="syn"&Tv1$makesCpG.tv1==1],Tv2$mean[Tv2$Type.tv2=="syn"&Tv2$makesCpG.tv2==1])
SynNoncpg<-c(Tv1$mean[Tv1$Type.tv1=="syn"&Tv1$makesCpG.tv1==0],Tv2$mean[Tv2$Type.tv2=="syn"&Tv2$makesCpG.tv2==0])
Nonsyncpg <- c(Tv1$mean[Tv1$Type.tv1=="nonsyn"&Tv1$makesCpG.tv1==1],Tv2$mean[Tv2$Type.tv2=="nonsyn"&Tv2$makesCpG.tv2==1])
NNcpg <- c(Tv1$mean[Tv1$Type.tv1=="nonsyn"&Tv1$makesCpG.tv1==0],Tv2$mean[Tv2$Type.tv2=="nonsyn"&Tv2$makesCpG.tv2==0])

r3<-wilcox.test(SynNoncpg,Syncpg, alternative = "greater", paired = FALSE) 
r4<-wilcox.test(NNcpg,Nonsyncpg, alternative = "greater", paired = FALSE) 

r3[[3]] #[1] 5.038601e-143
r4[[3]] #[1] 1.390009e-237


## 2.2) Run Wilcoxon test on each NT, transversion

dat1<-mf_filtered[[2]]
dat1<-dat1[dat1$pos>=342, ]
dat2<-mf_filtered[[3]]
dat2<-dat2[dat2$pos>=342, ]

m<-data.frame()
se<-data.frame()
m_nonCpG<-data.frame()
se_nonCpG<-data.frame()  
mr_CpG<-data.frame()
se_CpG<-data.frame()


for (typeofsite in c("syn", "nonsyn","stop")){
        for (wtnt in c("a", "t", "c", "g")){
                mr1<- dat1$mean[dat1$Type.tv1==typeofsite & dat1$ref==wtnt]
                mr2<- dat2$mean[dat2$Type.tv2==typeofsite & dat1$ref==wtnt]
                
                m_NonCpG1<-dat1$mean[dat1$Type.tv1==typeofsite & dat1$ref==wtnt & dat1$makesCpG.tv1==0]
                m_NonCpG2<-dat2$mean[dat2$Type.tv2==typeofsite & dat2$ref==wtnt & dat2$makesCpG.tv2==0]
                
                m_CpG1<-dat1$mean[dat1$Type.tv1==typeofsite & dat1$ref==wtnt & dat1$makesCpG.tv1==1]
                m_CpG2<-dat2$mean[dat2$Type.tv2==typeofsite & dat2$ref==wtnt & dat2$makesCpG.tv2==1]
                
                mr<-c(mr1,mr2)
                m_NonCpG<-c(m_NonCpG1,m_NonCpG2)
                m_CpG<-c(m_CpG1,m_CpG2)
                
                m[typeofsite,wtnt]<-mean(mr[!is.na(mr)])
                se[typeofsite,wtnt]<-std.error(mr[!is.na(mr)])
                m_nonCpG[typeofsite,wtnt]<-mean(m_NonCpG[!is.na(m_NonCpG)])
                se_nonCpG[typeofsite,wtnt]<-std.error(m_NonCpG[!is.na(m_NonCpG)])
                mr_CpG[typeofsite,wtnt]<-mean(m_CpG[!is.na(m_CpG)])
                se_CpG[typeofsite,wtnt]<-std.error(m_CpG[!is.na(m_CpG)])
                
                vectorname<-paste0(typeofsite,"_",wtnt)
                assign(vectorname, mr)
                vname1<<-paste0(typeofsite,"_",wtnt,"_noncpg")
                assign(vname1, m_NonCpG)
                vname2<<-paste0(typeofsite,"_",wtnt,"_cpg")
                assign(vname2, m_CpG)
        }
}


rownames(se)<-c("syn_se","nonsyn_se","stop_se")
rownames(m_nonCpG)<-c("syn_noncpg","nonsyn_noncpg","stop_noncpg")
rownames(se_nonCpG)<-c("syn_noncpg_se","nonsyn_noncpg_se","stop_noncpg_se")
rownames(mr_CpG)<-c("syn_cpg","nonsyn_cpg","stop_cpg")
rownames(se_CpG)<-c("syn_cpg_se","nonsyn_cpg_se","stop_cpg_se")

MFbyType.tv<-rbind(m,se,mr_CpG,se_CpG,m_nonCpG,se_nonCpG)
MFbyType.tv2<-t(MFbyType.tv)
write.csv(MFbyType.tv2,paste0("Output1A/SummaryStats/Transversion_MF_byNt_byType_mean.csv"))

# CpG vs nonCpG
WilcoxTest.results.nt.tv<-data.frame(matrix(ncol=3,nrow=8))
colnames(WilcoxTest.results.nt.tv)<-c("nt","test","P.value")

for (i in c("a","t","c","g")) {  
   syncpg<-get(paste0("syn_",i,"_cpg"))
   synnoncpg<-get(paste0("syn_",i,"_noncpg"))
   nonsyncpg<-get(paste0("nonsyn_",i,"_cpg"))
   nonsynnoncpg<-get(paste0("nonsyn_",i,"_noncpg"))
   if (i=="a"){
           if (length(syncpg)==0) next
           else{result1<-wilcox.test(syncpg, synnoncpg, alternative = "less", paired = FALSE) 
                result2<-wilcox.test(nonsyncpg,nonsynnoncpg,alternative = "less", paired = FALSE)
                   for (r in 1:2){
                           result<-get(paste0('result',r))
                           WilcoxTest.results.nt.tv$nt[r]<-i
                           WilcoxTest.results.nt.tv$test[r]<-result[[7]]
                           WilcoxTest.results.nt.tv$P.value[r]<-result[[3]]}
                   }}
   if (i=="t"){
           if (length(syncpg)==0) next
           else{ result3<-wilcox.test(syncpg, synnoncpg, alternative = "less", paired = FALSE) 
                 result4<-wilcox.test(nonsyncpg,nonsynnoncpg,alternative = "less", paired = FALSE)
                   for (r in 3:4){
                           result<-get(paste0('result',r))
                           WilcoxTest.results.nt.tv$nt[r]<-i
                           WilcoxTest.results.nt.tv$test[r]<-result[[7]]
                           WilcoxTest.results.nt.tv$P.value[r]<-result[[3]]}  
                   }}
   if (i=="c"){
           if (length(syncpg)==0) next
           else{ result5<-wilcox.test(syncpg, synnoncpg, alternative = "less", paired = FALSE) 
                 result6<-wilcox.test(nonsyncpg,nonsynnoncpg,alternative = "less", paired = FALSE)
                   for (r in 5:6){
                           result<-get(paste0('result',r))
                           WilcoxTest.results.nt.tv$nt[r]<-i
                           WilcoxTest.results.nt.tv$test[r]<-result[[7]]
                           WilcoxTest.results.nt.tv$P.value[r]<-result[[3]]}
                   }}
   if (i=="g"){
           if (length(syncpg)==0) next
           else{ result7<-wilcox.test(syncpg, synnoncpg, alternative = "less", paired = FALSE) 
                 result8<-wilcox.test(nonsyncpg,nonsynnoncpg,alternative = "less", paired = FALSE)
                   for (r in 7:8){
                           result<-get(paste0('result',r))
                           WilcoxTest.results.nt.tv$nt[r]<-i
                           WilcoxTest.results.nt.tv$test[r]<-result[[7]]
                           WilcoxTest.results.nt.tv$P.value[r]<-result[[3]]}         
                   }}     
   }
write.csv(WilcoxTest.results.nt.tv,paste0("Output1A/SummaryStats/WilcoxTestResults_TV_eachNT_mean.csv"))
 

#####################
## Not sure if we need to run these:
# 3) use all data to create summary / run Wilcoxon Test (using all values,not just the mean of 195 files) on Transition

source("Rscripts/MutationFreqSum.filtered.R")

##3.1 Transition 
nuc.mf2<-data.frame("s.cpg"=matrix(nrow=4))
rownames(nuc.mf2)<-c("a","t","c","g")

WilcoxTest.nt.mf2<-data.frame(matrix(ncol=3,nrow=6))
colnames(WilcoxTest.nt.mf2)<-c("nt","test","P.value")

k=1
for (i in c("A","T","C","G")) {
        if (i=="A"|i=="T"){
                syncpg1<-get(paste0(i,"_syn_cpg"))
                synnoncpg1<-get(paste0(i,"_syn_noncpg"))
                nonsyncpg1<-get(paste0(i,"_nonsyn_cpg"))
                nonsynnoncpg1<-get(paste0(i,"_nonsyn_noncpg"))
                
                nuc.mf2$s.cpg[k]<-mean(syncpg1, na.rm=T)
                nuc.mf2$s.cpg.se[k]<-std.error(syncpg1, na.rm=T)
                nuc.mf2$s.ncpg[k]<-mean(synnoncpg1, na.rm=T)
                nuc.mf2$s.ncpg.se[k]<-std.error(synnoncpg1, na.rm=T)
                nuc.mf2$ns.cpg[k]<-mean(nonsyncpg1, na.rm=T)
                nuc.mf2$ns.cpg.se[k]<-std.error(nonsyncpg1, na.rm=T)
                nuc.mf2$ns.ncpg[k]<-mean(nonsynnoncpg1, na.rm=T)
                nuc.mf2$ns.ncpg.se[k]<-std.error(nonsynnoncpg1, na.rm=T)
                #nuc.mf2$stop.ncpg[k]<-NA
                #nuc.mf2$stop.ncpg.se[k]<-NA
                
                
                if (i=="A"){
                        result1<-wilcox.test(syncpg1, synnoncpg1, alternative = "less", paired = FALSE) 
                        result2<-wilcox.test(nonsyncpg1,nonsynnoncpg1,alternative = "less", paired = FALSE)
                        for (r in 1:2){
                                result<-get(paste0('result',r))
                                WilcoxTest.nt.mf2$nt[r]<-i
                                WilcoxTest.nt.mf2$test[r]<-result[[7]]
                                WilcoxTest.nt.mf2$P.value[r]<-result[[3]]}
                }
                if (i=="T"){
                        result3<-wilcox.test(syncpg1, synnoncpg1, alternative = "less", paired = FALSE) 
                        result4<-wilcox.test(nonsyncpg1,nonsynnoncpg1,alternative = "less", paired = FALSE)
                        for (r in 3:4){
                                result<-get(paste0('result',r))
                                WilcoxTest.nt.mf2$nt[r]<-i
                                WilcoxTest.nt.mf2$test[r]<-result[[7]]
                                WilcoxTest.nt.mf2$P.value[r]<-result[[3]]}
                }
        }
        else { 
                
                synnoncpg1<-get(paste0(i,"_syn_noncpg"))
                nonsynnoncpg1<-get(paste0(i,"_nonsyn_noncpg"))
                
                nuc.mf2$s.cpg[k]<-NA
                nuc.mf2$s.cpg.se[k]<-NA
                nuc.mf2$s.ncpg[k]<-mean(synnoncpg1, na.rm=T)
                nuc.mf2$s.ncpg.se[k]<-std.error(synnoncpg1, na.rm=T)
                nuc.mf2$ns.cpg[k]<-NA
                nuc.mf2$ns.cpg.se[k]<-NA
                nuc.mf2$ns.ncpg[k]<-mean(nonsynnoncpg1, na.rm=T)
                nuc.mf2$ns.ncpg.se[k]<-std.error(nonsynnoncpg1, na.rm=T)
                 
                if (i =="C") {
                        result5<-wilcox.test(synnoncpg1,nonsynnoncpg1, alternative = "greater", paired = FALSE) 
                        WilcoxTest.nt.mf2$nt[5]<-i
                        WilcoxTest.nt.mf2$test[5]<-result5[[7]]
                        WilcoxTest.nt.mf2$P.value[5]<-result5[[3]]     }           
                if (i =="G") {
                        result6<-wilcox.test(synnoncpg1,nonsynnoncpg1, alternative = "greater", paired = FALSE) 
                        WilcoxTest.nt.mf2$nt[6]<-i
                        WilcoxTest.nt.mf2$test[6]<-result6[[7]]
                        WilcoxTest.nt.mf2$P.value[6]<-result6[[3]]    }
        }
        k=k+1
} 
write.csv(WilcoxTest.nt.mf2,paste0("Output/SummaryStats/WilcoxTestResults_Transition_eachNT_all.csv"))
write.csv(nuc.mf2,paste0("Output/SummaryStats/TransitionMF_byNT_byType_all.csv"))



wilcoxtest2<-data.frame("test"=matrix(nrow=4))
Typelist2<-list()
syn_all<-c(A_syn,T_syn,C_syn,G_syn)
Typelist2[[1]]<-syn_all; names(Typelist2)[1]<-"syn_all"
nonsyn_all<-c(A_nonsyn,T_nonsyn,C_nonsyn,G_nonsyn)
Typelist2[[2]]<-nonsyn_all; names(Typelist2)[2]<-"nonsyn_all"
stop_all<-c(A_stop,T_stop,C_stop,G_stop)
Typelist2[[3]]<-stop_all; names(Typelist2)[3]<-"stop_all"
syn_allCpG<-c(A_syn_cpg,T_syn_cpg)
Typelist2[[4]]<-syn_allCpG; names(Typelist2)[4]<-"syn_allCpG"
syn_allnonCpG<-c(A_syn_noncpg,T_syn_noncpg)
Typelist2[[5]]<-syn_allnonCpG; names(Typelist2)[5]<-"syn_allnonCpG"
nonsyn_allCpG<-c(A_nonsyn_cpg,T_nonsyn_cpg)
Typelist2[[6]]<-nonsyn_allCpG; names(Typelist2)[6]<-"nonsyn_allCpG"
nonsyn_allnonCpG<-c(A_nonsyn_noncpg,T_nonsyn_noncpg)
Typelist2[[7]]<-nonsyn_allnonCpG; names(Typelist2)[7]<-"nonsyn_allnonCpG"


re1<-wilcox.test(syn_all,nonsyn_all, alternative = "greater", paired = FALSE) 
wilcoxtest2$test[1]<-re1[[7]]
wilcoxtest2$P.value[1]<-re1[[3]]

re2<-wilcox.test(nonsyn_all,stop_all, alternative = "greater", paired = FALSE) 
wilcoxtest2$test[2]<-re2[[7]]
wilcoxtest2$P.value[2]<-re2[[3]]

re3<-wilcox.test(syn_allCpG,syn_allnonCpG, alternative = "less", paired = FALSE) 
wilcoxtest2$test[3]<-re3[[7]]
wilcoxtest2$P.value[3]<-re3[[3]]

re4<-wilcox.test(nonsyn_allCpG,nonsyn_allnonCpG, alternative = "less", paired = FALSE) 
wilcoxtest2$test[4]<-re4[[7]]
wilcoxtest2$P.value[4]<-re4[[3]]

write.csv(wilcoxtest2,"Output/SummaryStats/WilcoxTestResults_by_Type_MF_Ts_all.csv")


Type.mf1<-data.frame("mean"=matrix(nrow=7))

for (i in 1:7){
        rownames(Type.mf1)[i]<-names(Typelist2)[i]
        Type.mf1$mean[i]<-mean(Typelist2[[i]],na.rm=T)
        Type.mf1$se[i]<-std.error(Typelist2[[i]], na.rm=T)
        
}

write.csv(Type.mf1,paste0("Output/SummaryStats/MutFreq_byType_Summary_all.csv"))


##################################
##3.2 Transversion (all)
nuc.mft<-data.frame("s.cpg"=matrix(nrow=4))
rownames(nuc.mft)<-c("a","t","c","g")

WilcoxTest.nt.mft<-data.frame(matrix(ncol=3,nrow=8))
colnames(WilcoxTest.nt.mft)<-c("nt","test","P.value")

k=1
for (i in c("A","T","C","G")) {
        if (i=="A"|i=="G"){
                syncpg1<-get(paste0(i,"_tv1_syn_cpg"))
                #syncpg2<-get(paste0(i,"_tv2_syn_cpg"))
                synnoncpg1<-get(paste0(i,"_tv1_syn_noncpg"))
                synnoncpg2<-get(paste0(i,"_tv2_syn_noncpg"))
                nonsyncpg1<-get(paste0(i,"_tv1_nonsyn_cpg"))
                #nonsyncpg2<-get(paste0(i,"_tv2_nonsyn_cpg"))
                nonsynnoncpg1<-get(paste0(i,"_tv1_nonsyn_noncpg"))
                nonsynnoncpg2<-get(paste0(i,"_tv2_nonsyn_noncpg"))
                
                synnoncpg<-c(synnoncpg1,synnoncpg2)
                nonsynnoncpg<-c(nonsynnoncpg1,nonsynnoncpg2)
                
                nuc.mft$s.cpg[k]<-mean(syncpg1, na.rm=T)
                nuc.mft$s.cpg.se[k]<-std.error(syncpg1, na.rm=T)
                nuc.mft$s.ncpg[k]<-mean(synnoncpg, na.rm=T)
                nuc.mft$s.ncpg.se[k]<-std.error(synnoncpg, na.rm=T)
                nuc.mft$ns.cpg[k]<-mean(nonsyncpg1, na.rm=T)
                nuc.mft$ns.cpg.se[k]<-std.error(nonsyncpg1, na.rm=T)
                nuc.mft$ns.ncpg[k]<-mean(nonsynnoncpg, na.rm=T)
                nuc.mft$ns.ncpg.se[k]<-std.error(nonsynnoncpg, na.rm=T)
                #nuc.mft$stop.ncpg[k]<-NA
                #nuc.mft$stop.ncpg.se[k]<-NA
                
                
                if (i=="A"){
                        result1<-wilcox.test(syncpg1, synnoncpg, alternative = "less", paired = FALSE) 
                        result2<-wilcox.test(nonsyncpg1,nonsynnoncpg,alternative = "less", paired = FALSE)
                        for (r in 1:2){
                                result<-get(paste0('result',r))
                                WilcoxTest.nt.mft$nt[r]<-i
                                WilcoxTest.nt.mft$test[r]<-result[[7]]
                                WilcoxTest.nt.mft$P.value[r]<-result[[3]]}
                }
                if (i=="G"){
                        result7<-wilcox.test(syncpg1, synnoncpg, alternative = "less", paired = FALSE) 
                        result8<-wilcox.test(nonsyncpg1,nonsynnoncpg,alternative = "less", paired = FALSE)
                        for (r in 7:8){
                                result<-get(paste0('result',r))
                                WilcoxTest.nt.mft$nt[r]<-i
                                WilcoxTest.nt.mft$test[r]<-result[[7]]
                                WilcoxTest.nt.mft$P.value[r]<-result[[3]]}
                }
        }
        if (i=="T"|i=="C") { 
                
                syncpg2<-get(paste0(i,"_tv2_syn_cpg"))
                synnoncpg1<-get(paste0(i,"_tv1_syn_noncpg"))
                synnoncpg2<-get(paste0(i,"_tv2_syn_noncpg"))
                nonsyncpg2<-get(paste0(i,"_tv2_nonsyn_cpg"))
                nonsynnoncpg1<-get(paste0(i,"_tv1_nonsyn_noncpg"))
                nonsynnoncpg2<-get(paste0(i,"_tv2_nonsyn_noncpg"))
                
                synnoncpg<-c(synnoncpg1,synnoncpg2)
                nonsynnoncpg<-c(nonsynnoncpg1,nonsynnoncpg2)
                
                nuc.mft$s.cpg[k]<-mean(syncpg2, na.rm=T)
                nuc.mft$s.cpg.se[k]<-std.error(syncpg2, na.rm=T)
                nuc.mft$s.ncpg[k]<-mean(synnoncpg, na.rm=T)
                nuc.mft$s.ncpg.se[k]<-std.error(synnoncpg, na.rm=T)
                nuc.mft$ns.cpg[k]<-mean(nonsyncpg2, na.rm=T)
                nuc.mft$ns.cpg.se[k]<-std.error(nonsyncpg2, na.rm=T)
                nuc.mft$ns.ncpg[k]<-mean(nonsynnoncpg, na.rm=T)
                nuc.mft$ns.ncpg.se[k]<-std.error(nonsynnoncpg, na.rm=T)
                
                if (i =="T") {
                        result3<-wilcox.test(syncpg2, synnoncpg, alternative = "less", paired = FALSE) 
                        result4<-wilcox.test(nonsyncpg2,nonsynnoncpg,alternative = "less", paired = FALSE)
                        for (r in 3:4){
                                result<-get(paste0('result',r))
                                WilcoxTest.nt.mft$nt[r]<-i
                                WilcoxTest.nt.mft$test[r]<-result[[7]]
                                WilcoxTest.nt.mft$P.value[r]<-result[[3]]}
                         }           
                if (i =="C") {
                        result5<-wilcox.test(syncpg2, synnoncpg, alternative = "less", paired = FALSE) 
                        result6<-wilcox.test(nonsyncpg2,nonsynnoncpg,alternative = "less", paired = FALSE)
                        for (r in 5:6){
                                result<-get(paste0('result',r))
                                WilcoxTest.nt.mft$nt[r]<-i
                                WilcoxTest.nt.mft$test[r]<-result[[7]]
                                WilcoxTest.nt.mft$P.value[r]<-result[[3]]}   
                        }
        }
        k=k+1
} 

write.csv(WilcoxTest.nt.mft,paste0("Output/SummaryStats/WilcoxTestResults_TV_eachNT_all.csv"))
write.csv(nuc.mft,paste0("Output/SummaryStats/TransversionMF_byNT_byType_all.csv"))



wilcoxtest3<-data.frame("test"=matrix(nrow=4))
Typelist3<-list()
syn_all<-c(A_tv1_syn,T_tv1_syn,C_tv1_syn,G_tv1_syn,A_tv2_syn,T_tv2_syn,C_tv2_syn,G_tv2_syn)
Typelist3[[1]]<-syn_all; names(Typelist3)[1]<-"syn_all"
nonsyn_all<-c(A_tv1_nonsyn,T_tv1_nonsyn,C_tv1_nonsyn,G_tv1_nonsyn,A_tv2_nonsyn,T_tv2_nonsyn,C_tv2_nonsyn,G_tv2_nonsyn)
Typelist3[[2]]<-nonsyn_all; names(Typelist3)[2]<-"nonsyn_all"
syn_allCpG<-c(A_tv1_syn_cpg,T_tv2_syn_cpg,G_tv1_syn_cpg,C_tv2_syn_cpg )
Typelist3[[3]]<-syn_allCpG; names(Typelist3)[3]<-"syn_allCpG"
syn_allnonCpG<-c(A_tv1_syn_noncpg,T_tv1_syn_noncpg,G_tv1_syn_noncpg,C_tv1_syn_noncpg,A_tv2_syn_noncpg,T_tv2_syn_noncpg,G_tv2_syn_noncpg,C_tv2_syn_noncpg)
Typelist3[[4]]<-syn_allnonCpG; names(Typelist3)[4]<-"syn_allnonCpG"
nonsyn_allCpG<-c(A_tv1_nonsyn_cpg,T_tv2_nonsyn_cpg,G_tv1_nonsyn_cpg,C_tv2_nonsyn_cpg)
Typelist3[[5]]<-nonsyn_allCpG; names(Typelist3)[5]<-"nonsyn_allCpG"
nonsyn_allnonCpG<-c(A_tv1_nonsyn_noncpg,T_tv1_nonsyn_noncpg,G_tv1_nonsyn_noncpg,C_tv1_nonsyn_noncpg,A_tv2_nonsyn_noncpg,T_tv2_nonsyn_noncpg,G_tv2_nonsyn_noncpg,C_tv2_nonsyn_noncpg)
Typelist3[[6]]<-nonsyn_allnonCpG; names(Typelist3)[6]<-"nonsyn_allnonCpG"

stop_all<-c(A_tv1_stop,T_tv1_stop,C_tv1_stop,G_tv1_stop,A_tv2_stop,T_tv2_stop,C_tv2_stop,G_tv2_stop)
Typelist3[[7]]<-stop_all; names(Typelist3)[7]<-"stop_all"



re1<-wilcox.test(syn_all,nonsyn_all, alternative = "greater", paired = FALSE) 
wilcoxtest3$test[1]<-re1[[7]]
wilcoxtest3$P.value[1]<-re1[[3]]

re2<-wilcox.test(nonsyn_all,stop_all, alternative = "greater", paired = FALSE) 
wilcoxtest3$test[2]<-re2[[7]]
wilcoxtest3$P.value[2]<-re2[[3]]

re3<-wilcox.test(syn_allCpG,syn_allnonCpG, alternative = "less", paired = FALSE) 
wilcoxtest3$test[3]<-re3[[7]]
wilcoxtest3$P.value[3]<-re3[[3]]

re4<-wilcox.test(nonsyn_allCpG,nonsyn_allnonCpG, alternative = "less", paired = FALSE) 
wilcoxtest3$test[4]<-re4[[7]]
wilcoxtest3$P.value[4]<-re4[[3]]

write.csv(wilcoxtest3,"Output/SummaryStats/WilcoxTestResults_by_Type_MF_Transversion_all.csv")


Type.mf2<-data.frame("mean"=matrix(nrow=7))

for (i in 1:7){
        rownames(Type.mf2)[i]<-names(Typelist3)[i]
        Type.mf2$mean[i]<-mean(Typelist3[[i]],na.rm=T)
        Type.mf2$se[i]<-std.error(Typelist3[[i]], na.rm=T)
}

write.csv(Type.mf2,paste0("Output/SummaryStats/MutFreq_byType_Summary_Transversion_all.csv"))


####### plot to visually assure ###

Type.mf2$names<-c("Syn","Nonsyn","Syn_CpG", "Syn_nonCpG","Nonsyn_CpG","Nonsyn_nonCpG","Stop")
x=seq(1,7, by=1)

plot(x, Type.mf2$mean, xaxt="n",main="", xlab ="",ylim=c(0,0.0012),pch=".")
#axis(1, at=1:7, labels=Type.mf2$names, las = 1, cex.axis = 0.8)
text(cex=1,x=x, y=-0.0001, labels=Type.mf2$names, xpd=TRUE, srt=35,adj= 1)
bar<-0.02

segments(x,(Type.mf2$mean-Type.mf2$se),x,(Type.mf2$mean+Type.mf2$se))
segments(x-bar,(Type.mf2$mean-Type.mf2$se),x+bar,(Type.mf2$mean-Type.mf2$se))
segments(x-bar,(Type.mf2$mean+Type.mf2$se),x+bar,(Type.mf2$mean+Type.mf2$se))

##################################




