library(reshape)
library(tidyverse)
library(zoo)
library(purrr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(sfsmisc)
source("Rscripts/baseRscript.R")


HCVFiles_overview3<-list.files("Output1A/Overview3/",pattern="overview3.csv")

##### read the summary mutation frequency files before removing the sites with >66% NA
TsMutFreq  <-read.csv("Output1A/MutFreq.filtered/Ts_MutFreq_1A_all.csv",stringsAsFactors = F,row.names=1)
#remove the sites before the coding region
Ts<-TsMutFreq[TsMutFreq$pos>=342,]
#attach the metadata
M<-read.csv("Output1A/Overview3/D75002-_overview3.csv", stringsAsFactors = F,row.names=1)
mut<-M[,c("pos","ref","Type","WTAA","MUTAA","bigAAChange","makesCpG")]
Ts<-merge(TsMutFreq,mut, by="pos")
Ts_NA<-Ts
Ts_zero<-Ts
# 1.Filter out Mut Freq < 0.01  -> replaced with NA
mf<-0.01
Ts_NA[Ts_NA<mf]<-NA

#Filter out Mut Freq < 0.001  -> replaced with 0
Ts_zero[Ts_zero<mf]<-0

#calculate the mean
Ts_NA$mean<-rowMeans(Ts_NA[2:(s+1)], na.rm=T)
mean(Ts_NA$mean,na.rm=T) #0.02755865

Ts_zero$mean<-rowMeans(Ts_zero[2:(s+1)],na.rm=T)
mean(Ts_zero$mean) #0.002039763

# Ts 
Ts$mean<-rowMeans(Ts[2:(s+1)], na.rm=T)


write.csv(Ts_NA,   "Output1A/Q35Compare/Summary_Ts_NA0.01.csv")
write.csv(Ts_zero, "Output1A/Q35Compare/Summary_Ts_zero0.01.csv")

#########################################################################################
##Plot across the genome

endnuc<-Ts$pos[nrow(Ts)]
SNPFreq<-Ts_NA

n<-data.frame("pos"=c(342:endnuc))
SNPFreqs<-merge(n,SNPFreq,by="pos",all.x=T)
pdf(paste0("Output1A/Q35Compare/Ts.MutFreq.NAreplaced0.01.pdf"),width=15,height=7.5)
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

n<-data.frame("pos"=c(342:endnuc))
SNPFreqs<-merge(n,SNPFreq,by="pos",all.x=T)
pdf(paste0("Output1A/Q35Compare/Ts.MutFreq.Zero.replaced0.01.pdf"),width=15,height=7.5)
plot(mean~pos, data=SNPFreqs,t="n",log='y',yaxt='n',xlab='Genome position (H77)',ylab="Average transition mutation frequency",
     main="Transition muttaion frequency ",ylim=c(0.00001,0.1),xlim=c(340,8500))
eaxis(side = 2, at = 10^((-1):(-(5))), cex=2)
for(i in 2:5){abline(h = 1:10 * 10^(-i), col = "gray80")}
points(mean~pos, data=SNPFreqs,pch=20,col="#EE667766",cex=0.5)

#add rolling average
roll100<-rollmean(SNPFreqs$mean, k=100, na.rm=T)
SNPFreqs$roll100<-c(rep(NA, times=50),roll100,c(rep(NA, times=49)))
lines(roll100~pos,data=SNPFreqs, col="#AA3377")

dev.off()

###########################################################################
SumTab<-data.frame("data"=c("Ts","Ts_zero","Ts_NA"))
t<-list()
t[[1]]<-Ts
t[[2]]<-Ts_zero
t[[3]]<-Ts_NA
names(t)<-c("Ts","Ts_zero","Ts_NA")

for (i in 1:3){
        dat<-t[[i]]
        SumTab$mean[i]<-mean(dat$mean, na.rm=T)
        SumTab$Syn[i]<-mean(dat$mean[dat$Type=="syn"],na.rm=T) 
        SumTab$Nonsyn[i]<-mean(dat$mean[dat$Type=="nonsyn"],na.rm=T)
        SumTab$Stop[i]<-mean(dat$mean[dat$Type=="stop"],na.rm=T)
        
        T2<-dat[dat$ref=="a"|dat$ref=="t",]
        SumTab$CpGmaking_Syn[i]  <-mean(T2$mean[T2$Type=="syn"&T2$makesCpG==1],na.rm=T)
        SumTab$NoCpGmaking_Syn[i]<-mean(T2$mean[T2$Type=="syn"&T2$makesCpG==0],na.rm=T)
        if (is.nan(SumTab$CpGmaking_Syn[i])|is.nan(SumTab$NoCpGmaking_Syn[i])) SumTab$P_value.Wilcoxon.Test_syn[i]<-NA
        else {
                r3<-wilcox.test(T2$mean[T2$Type=="syn"&T2$makesCpG==0], T2$mean[T2$Type=="syn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 
                SumTab$P_value.Wilcoxon.Test_syn[i]<-r3[[3]]
                
        }          
        SumTab$CpGmaking_Nonsyn[i]<-mean(Ts$mean[T2$Type=="nonsyn"&T2$makesCpG==1],na.rm=T)
        SumTab$NoCpGmaking_Nonsyn[i]<-mean(T2$mean[T2$Type=="nonsyn"&T2$makesCpG==0],na.rm=T)
        if (is.nan(SumTab$CpGmaking_Nonsyn[i])|is.nan(SumTab$NoCpGmaking_Nonsyn[i])) SumTab$P_value.Wilcoxon.Test_nonsyn[i]<-NA
        else {
                r4<-wilcox.test(T2$mean[T2$Type=="nonsyn"&T2$makesCpG==0], T2$mean[T2$Type=="nonsyn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 
                SumTab$P_value.Wilcoxon.Test_nonsyn[i]<-r4[[3]]
        }
}

write.csv(SumTab, "Output1A/Q35Compare/Tsfiltered_zero_NA_compare.summary0.01.csv")







###################################################################
# 1) NA replaced Transition
Ts<-Ts_NA

mean(Ts$mean[Ts$Type=="syn"],na.rm=T) #0.00911283  vs. #0.03059845 (removed <0.01)
mean(Ts$mean[Ts$Type=="nonsyn"],na.rm=T) #0.003973061  vs #0.02175058 (removed <0.01)
mean(Ts$mean[Ts$Type=="stop"],na.rm=T) #0.002562846  vs.0.02039543(removed <0.01)

## Compare Syn vs Nonsyn, Nonsyn vs. Stop
r1<-wilcox.test(Ts$mean[Ts$Type=="syn"], Ts$mean[Ts$Type=="nonsyn"], alternative = "greater", paired = FALSE) #W = 11795000, p-value < 2.2e-16
r2<-wilcox.test(Ts$mean[Ts$Type=="nonsyn"], Ts$mean[Ts$Type=="stop"], alternative = "greater", paired = FALSE) 

# syn vs nonsyn
r1[[3]]  #0
r2[[3]]  #0.09414204

## Compare CpG making vs. Non-CpGmaking 
T2<-Ts[Ts$ref=="a"|Ts$ref=="t",]
r3<-wilcox.test(T2$mean[T2$Type=="syn"&T2$makesCpG==0], T2$mean[T2$Type=="syn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 
r4<-wilcox.test(T2$mean[T2$Type=="nonsyn"&T2$makesCpG==0], T2$mean[T2$Type=="nonsyn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 

#syn -CpG
r3[[3]] #0.2197818
#nonsyn CpG
r4[[3]] #0.7039458

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
write.csv(MFbyType2,"Output/SummaryStats/TsMF_byNt_byType_mean_NAreplaced0.01.csv")

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
write.csv(WilcoxTest.results.nt,paste0("Output/SummaryStats/WilcoxTestResults_Ts_byNt_mean_NAreplaced0.01.csv"))

########
# 2) zero replaced
Ts<-Ts_zero

mean(Ts$mean[Ts$Type=="syn"],na.rm=T) #0.00911283  vs. #0.005881737 (replaced <0.01)
mean(Ts$mean[Ts$Type=="nonsyn"],na.rm=T) #0.003973061  vs # 0.000678216 (replaced <0.01)
mean(Ts$mean[Ts$Type=="stop"],na.rm=T) #0.002562846  vs. 0.0001589955 (replaced <0.01)

## Compare Syn vs Nonsyn, Nonsyn vs. Stop
r1<-wilcox.test(Ts$mean[Ts$Type=="syn"], Ts$mean[Ts$Type=="nonsyn"], alternative = "greater", paired = FALSE) #W = 11795000, p-value < 2.2e-16
r2<-wilcox.test(Ts$mean[Ts$Type=="nonsyn"], Ts$mean[Ts$Type=="stop"], alternative = "greater", paired = FALSE) 

# syn vs nonsyn
r1[[3]]  #0
r2[[3]]  #3.186742e-31

## Compare CpG making vs. Non-CpGmaking 
T2<-Ts[Ts$ref=="a"|Ts$ref=="t",]
r3<-wilcox.test(T2$mean[T2$Type=="syn"&T2$makesCpG==0], T2$mean[T2$Type=="syn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 
r4<-wilcox.test(T2$mean[T2$Type=="nonsyn"&T2$makesCpG==0], T2$mean[T2$Type=="nonsyn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 

#syn -CpG
r3[[3]] #0.01818175
#nonsyn CpG
r4[[3]] #1.372421e-06

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
write.csv(MFbyType2,"Output/SummaryStats/TsMF_byNt_byType_mean_ZEROreplaced0.01.csv")

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
write.csv(WilcoxTest.results.nt,paste0("Output/SummaryStats/WilcoxTestResults_Ts_byNt_mean_ZEROreplaced0.01.csv"))


###################################################################
###################################################################
#Plot summary of 1) Mut Freq _compare CpG creating vs. non-CpGcreating   

#1.Zero-replaced (0.01)
TS<-Ts_zero

#1. Transition SYN
k=1
transMF<-list()
for (i in c("a","t")) {
        for (cp in c(0,1)){
                datavector<-TS$mean[TS$Type=="syn" & TS$ref==i &TS$makesCpG==cp]
                nt<-toupper(i)
                if (cp==0) cpg<-''
                if (cp==1) cpg<-"(CpG)"
                vname<-paste0(nt,cpg)
                dat<-data.frame(base=rep(vname,times=length(datavector)), mf=datavector)
                transMF[[k]]<-dat
                names(transMF)[k]<-vname
                k=k+1
        }
}

mfdata1<-do.call(rbind, transMF)

z=c(0.7,0.3,0.7,0.3)
mfplot1<-ggplot(mfdata1,aes(x=base,y=mf,fill=base))+geom_boxplot(alpha=z)+labs(x="Nucleotide",y="Mutation frequency")+
        scale_fill_manual(values=c("#66CCEEB3","#66CCEE4D","#228833B3","#2288334D")) + theme_classic()+
        guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10))+
        ggtitle("CpG vs. nonCpG syn transition") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
        theme(plot.title = element_text(hjust = 0.5))+
        geom_vline(xintercept = 2.5)
ggsave(filename="Output/SummaryFig.Filtered/Zero_replaced/CpGvsNonCpG_synTs_ZeroReplaced0.01.pdf",width=4, height=4, units='in',device='pdf', plot=mfplot1)

### Transition Nonsyn

k=1
transMF2<-list()
for (i in c("a","t")) {
        for (cp in c(0,1)){
                datavector<-TS$mean[TS$Type=="nonsyn" & TS$ref==i &TS$makesCpG==cp]
                nt<-toupper(i)
                if (cp==0) cpg<-''
                if (cp==1) cpg<-"(CpG)"
                vname<-paste0(nt,cpg)
                dat<-data.frame(base=rep(vname,times=length(datavector)), mf=datavector)
                transMF2[[k]]<-dat
                names(transMF2)[k]<-vname
                k=k+1
        }
}

mfdata2<-do.call(rbind, transMF2)

z=c(0.7,0.3,0.7,0.3)
mfplot2<-ggplot(mfdata2,aes(x=base,y=mf,fill=base))+geom_boxplot(alpha=z)+labs(x="Nucleotide",y="Mutatin frequency")+
        scale_fill_manual(values=c("#66CCEEB3","#66CCEE4D","#228833B3","#2288334D")) + theme_classic()+
        guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10))+
        ggtitle("CpG vs. nonCpG nonsyn transition") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
        theme(plot.title = element_text(hjust = 0.5))+
        geom_vline(xintercept = 2.5)
ggsave(filename="Output/SummaryFig.Filtered/Zero_replaced/CpGvsNonCpG_nonsyn_Transition_ZeroReplaced0.01.pdf",width=4, height=4, units='in',device='pdf', plot=mfplot2)


######
##2) NA replaced (0.01)

#1.Zero-replaced (0.01)
TS<-Ts_NA

#1. Transition SYN
k=1
transMF<-list()
for (i in c("a","t")) {
        for (cp in c(0,1)){
                datavector<-TS$mean[TS$Type=="syn" & TS$ref==i &TS$makesCpG==cp]
                nt<-toupper(i)
                if (cp==0) cpg<-''
                if (cp==1) cpg<-"(CpG)"
                vname<-paste0(nt,cpg)
                dat<-data.frame(base=rep(vname,times=length(datavector)), mf=datavector)
                transMF[[k]]<-dat
                names(transMF)[k]<-vname
                k=k+1
        }
}

mfdata1<-do.call(rbind, transMF)

z=c(0.7,0.3,0.7,0.3)
mfplot1<-ggplot(mfdata1,aes(x=base,y=mf,fill=base))+geom_boxplot(alpha=z)+labs(x="Nucleotide",y="Mutation frequency")+
        scale_fill_manual(values=c("#66CCEEB3","#66CCEE4D","#228833B3","#2288334D")) + theme_classic()+
        guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10))+
        ggtitle("CpG vs. nonCpG syn transition") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
        theme(plot.title = element_text(hjust = 0.5))+
        geom_vline(xintercept = 2.5)
ggsave(filename="Output/SummaryFig.Filtered/NA_replaced/CpGvsNonCpG_synTs_NAReplaced0.01.pdf",width=4, height=4, units='in',device='pdf', plot=mfplot1)

### Transition Nonsyn

k=1
transMF2<-list()
for (i in c("a","t")) {
        for (cp in c(0,1)){
                datavector<-TS$mean[TS$Type=="nonsyn" & TS$ref==i &TS$makesCpG==cp]
                nt<-toupper(i)
                if (cp==0) cpg<-''
                if (cp==1) cpg<-"(CpG)"
                vname<-paste0(nt,cpg)
                dat<-data.frame(base=rep(vname,times=length(datavector)), mf=datavector)
                transMF2[[k]]<-dat
                names(transMF2)[k]<-vname
                k=k+1
        }
}

mfdata2<-do.call(rbind, transMF2)

z=c(0.7,0.3,0.7,0.3)
mfplot2<-ggplot(mfdata2,aes(x=base,y=mf,fill=base))+geom_boxplot(alpha=z)+labs(x="Nucleotide",y="Mutatin frequency")+
        scale_fill_manual(values=c("#66CCEEB3","#66CCEE4D","#228833B3","#2288334D")) + theme_classic()+
        guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10))+
        ggtitle("CpG vs. nonCpG nonsyn transition") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
        theme(plot.title = element_text(hjust = 0.5))+
        geom_vline(xintercept = 2.5)
ggsave(filename="Output/SummaryFig.Filtered/NA_replaced/CpGvsNonCpG_nonsyn_Transition_NAReplaced0.01.pdf",width=4, height=4, units='in',device='pdf', plot=mfplot2)






