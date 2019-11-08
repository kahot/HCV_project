library(zoo)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(ggthemes)
library(sfsmisc)
source("Rscripts/baseRscript.R")

#SC<-list()
#fnames<-c("Ts", "Ts_NA", "Ts_zero")
mutrates<-read.csv("Data/Geller.mutation.rates.csv")


df<-read.csv("Output1A/Mutfreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F,row.names = 1 )
df<-df[,c(1,197:203)]
        
df$TSmutrate[df$ref=="a"]<-mutrates$mut.rate[mutrates$mutations=="AG"]
df$TSmutrate[df$ref=="c"]<-mutrates$mut.rate[mutrates$mutations=="CU"]
df$TSmutrate[df$ref=="g"]<-mutrates$mut.rate[mutrates$mutations=="GA"]
df$TSmutrate[df$ref=="t"]<-mutrates$mut.rate[mutrates$mutations=="UC"]
df$EstSC<-""
for (j in 1:nrow(df)){
        df$EstSC[j] <- EstimatedS(df$TSmutrate[j],df$mean[j])
}
df$EstSC<-as.numeric(df$EstSC)
write.csv(df,"Output1A/SelCoeff/SC.csv")


### Wilcoxon Test on SC  (use mean SC)

## Use A & T only for CpG Analysis
ty<-which(colnames(df)=="Type")
dat<-df
#fname=names(SC)[i]
k=1
TypeList<-list()
for (typeofsite in c("syn", "nonsyn","stop")){
        all<-dat$EstSC[dat[,ty]==typeofsite]
        allcpg<-dat$EstSC[dat[,ty]==typeofsite & dat$makesCpG==1]
        allnoncpg1<-dat$EstSC[dat[,ty]==typeofsite & dat$makesCpG==0 & dat$ref=="t"]
        allnoncpg2<-dat$EstSC[dat[,ty]==typeofsite & dat$makesCpG==0 & dat$ref=="a"]
        for (wtnt in c("a", "t", "c", "g")){
                selco<- dat$EstSC[dat[,ty]==typeofsite & dat$ref==wtnt]
                sc_NonCpG<-dat$EstSC[dat[,ty]==typeofsite & dat$ref==wtnt & dat$makesCpG==0]
                sc_CpG<-dat$EstSC[dat[,ty]==typeofsite & dat$ref==wtnt & dat$makesCpG==1]
                vectorname<-paste0(typeofsite,"_",wtnt)
                assign(vectorname, selco)
                vname1<<-paste0(typeofsite,"_",wtnt,"_noncpg")
                assign(vname1, sc_NonCpG)
                vname2<<-paste0(typeofsite,"_",wtnt,"_cpg")
                assign(vname2, sc_CpG)
                }
        typev<-paste0(typeofsite,"_all")
        assign(typev, all)
        TypeList[[k]]<-all
        names(TypeList)[k]<-typev
        k=k+1
        cpgv1<-paste0(typeofsite,"_allCpG")
        assign(cpgv1, allcpg)
        TypeList[[k]]<-allcpg
        names(TypeList)[k]<-cpgv1
        k=k+1
        cpgv0<-paste0(typeofsite,"_allnonCpG")
        allnoncpg<-c(allnoncpg1,allnoncpg2)
        assign(cpgv0, allnoncpg)
        TypeList[[k]]<-allnoncpg
        names(TypeList)[k]<-cpgv0
        k=k+1
        }

#wilcox.test 
wilcoxtest1<-data.frame("test"=matrix(nrow=4))
        
re1<-wilcox.test(syn_all,nonsyn_all, alternative = "less", paired = FALSE) 
wilcoxtest1$test[1]<-re1[[7]]
wilcoxtest1$P.value[1]<-re1[[3]]
        
re2<-wilcox.test(nonsyn_all,stop_all, alternative = "less", paired = FALSE) 
wilcoxtest1$test[2]<-re2[[7]]
wilcoxtest1$P.value[2]<-re2[[3]]
        
re3<-wilcox.test(syn_allCpG,syn_allnonCpG, alternative = "greater", paired = FALSE) 
wilcoxtest1$test[3]<-re3[[7]]
wilcoxtest1$P.value[3]<-re3[[3]]
        
re4<-wilcox.test(nonsyn_allCpG,nonsyn_allnonCpG, alternative = "greater", paired = FALSE) 
wilcoxtest1$test[4]<-re4[[7]]
wilcoxtest1$P.value[4]<-re4[[3]]
        
write.csv(wilcoxtest1,"Output1A/SelCoeff/WilcoxonResults_by_Type.csv")

Type.sc<-data.frame("mean"=matrix(nrow=9))
for (i in 1:9){
        rownames(Type.sc)[i]<-names(TypeList)[i]
        Type.sc$mean[i]<-mean(TypeList[[i]],na.rm=T)
        Type.sc$se[i]<-std.error(TypeList[[i]],na.rm=T)
}
write.csv(Type.sc,"Output1A/SelCoeff/SC_byType_Summary.csv")
      
  
## Test on NT by NT
S<-data.frame()
se<-data.frame()
S_CpG<-data.frame()
se_CpG<-data.frame()
S_nonCpG<-data.frame()
se_nonCpG<-data.frame()
for (typeofsite in c("syn", "nonsyn","stop")){
     for (wtnt in c("a", "t", "c", "g")){
        sc<- dat$EstSC[dat[,ty]==typeofsite & dat$ref==wtnt]
        S[typeofsite,wtnt]<-mean(sc[!is.na(sc)])
        se[typeofsite,wtnt]<-std.error(sc[!is.na(sc)])
        
        S_NonCpG<-dat$EstSC[dat$Type==typeofsite & dat$ref==wtnt & dat$makesCpG==0]
        S_nonCpG[typeofsite,wtnt]<-mean(S_NonCpG[!is.na(S_NonCpG)])
        se_nonCpG[typeofsite,wtnt]<-std.error(S_NonCpG[!is.na(S_NonCpG)])
        
        Sc_CpG<-dat$EstSC[dat$Type==typeofsite & dat$ref==wtnt & dat$makesCpG==1]
        S_CpG[typeofsite,wtnt]<-mean(Sc_CpG[!is.na(Sc_CpG)])
        se_CpG[typeofsite,wtnt]<-std.error(Sc_CpG[!is.na(Sc_CpG)])
        
        vectorname<-paste0(typeofsite,"_",wtnt)
        assign(vectorname, sc)
        vname1<<-paste0(typeofsite,"_",wtnt,"_noncpg")
        assign(vname1, S_NonCpG)
        vname2<<-paste0(typeofsite,"_",wtnt,"_cpg")
        assign(vname2, Sc_CpG)
     }
}
rownames(se)<-c("syn_se","nonsyn_se","stop_se")
rownames(S_nonCpG)<-c("syn_noncpg","nonsyn_noncpg","stop_noncpg")
rownames(se_nonCpG)<-c("syn_noncpg_se","nonsyn_noncpg_se","stop_noncpg_se")
rownames(S_CpG)<-c("syn_cpg","nonsyn_cpg","stop_cpg")
rownames(se_CpG)<-c("syn_cpg_se","nonsyn_cpg_se","stop_cpg_se")

SCbyType<-rbind(S,se,S_nonCpG,se_nonCpG,S_CpG,se_CpG)
SCbyType2<-data.frame(t(SCbyType))
write.csv(SCbyType2,"Output1A/SelCoeff/SC_byNt_byType.csv")

#run Wilcoxin Test  
nuc.sc<-data.frame("syn.cpg"=matrix(nrow=4))
rownames(nuc.sc)<-c("a","t","c","g")


WilcoxTest.nt.sc<-data.frame(matrix(ncol=3,nrow=6))
colnames(WilcoxTest.nt.sc)<-c("nt","test","P.value")
k=1
for (i in c("a","t","c","g")) {
    if (i=="a"|i=="t"){
        syncpg<-get(paste0("syn_",i,"_cpg"))
        synnoncpg<-get(paste0("syn_",i,"_noncpg"))
        nonsyncpg<-get(paste0("nonsyn_",i,"_cpg"))
        nonsynnoncpg<-get(paste0("nonsyn_",i,"_noncpg"))
        
        nuc.sc$syn.cpg[k]<-mean(syncpg, na.rm=T)
        nuc.sc$syn.cpg.se[k]<-std.error(syncpg, na.rm=T)
        nuc.sc$syn.ncpg[k]<-mean(synnoncpg, na.rm=T)
        nuc.sc$syn.ncpg.se[k]<-std.error(synnoncpg, na.rm=T)
        nuc.sc$ns.cpg[k]<-mean(nonsyncpg, na.rm=T)
        nuc.sc$ns.cpg.se[k]<-std.error(nonsyncpg, na.rm=T)
        nuc.sc$ns.ncpg[k]<-mean(nonsynnoncpg, na.rm=T)
        nuc.sc$ns.ncpg.se[k]<-std.error(nonsynnoncpg, na.rm=T)
        nuc.sc$stop.ncpg[k]<-NA
        nuc.sc$stop.ncpg.se[k]<-NA
        
        if (i=="a"){
                result1<-wilcox.test(syncpg, synnoncpg, alternative = "greater", paired = FALSE) 
                result2<-wilcox.test(nonsyncpg,nonsynnoncpg,alternative = "greater", paired = FALSE)
                for (r in 1:2){
                        result<-get(paste0('result',r))
                        WilcoxTest.nt.sc$nt[r]<-i
                        WilcoxTest.nt.sc$test[r]<-result[[7]]
                        WilcoxTest.nt.sc$P.value[r]<-result[[3]]}
        }
        if (i=="t"){
                result3<-wilcox.test(syncpg, synnoncpg, alternative = "greater", paired = FALSE) 
                result4<-wilcox.test(nonsyncpg,nonsynnoncpg,alternative = "greater", paired = FALSE)
                for (r in 3:4){
                        result<-get(paste0('result',r))
                        WilcoxTest.nt.sc$nt[r]<-i
                        WilcoxTest.nt.sc$test[r]<-result[[7]]
                        WilcoxTest.nt.sc$P.value[r]<-result[[3]]}
        }
    }
else {  synnoncpg<-get(paste0("syn_",i,"_noncpg"))
        nonsynnoncpg<-get(paste0("nonsyn_",i,"_noncpg"))
        stopnoncpg<-get(paste0("stop_",i,"_noncpg"))
        
        nuc.sc$syn.cpg[k]<-NA
        nuc.sc$syn.cpg.se[k]<-NA
        nuc.sc$syn.ncpg[k]<-mean(synnoncpg, na.rm=T)
        nuc.sc$syn.ncpg.se[k]<-std.error(synnoncpg, na.rm=T)
        nuc.sc$ns.cpg[k]<-NA
        nuc.sc$ns.cpg.se[k]<-NA
        nuc.sc$ns.ncpg[k]<-mean(nonsynnoncpg, na.rm=T)
        nuc.sc$ns.ncpg.se[k]<-std.error(nonsynnoncpg, na.rm=T)
        nuc.sc$stop.ncpg[k]<-mean(stopnoncpg, na.rm=T)
        nuc.sc$stop.ncpg.se[k]<-std.error(stopnoncpg, na.rm=T)
        
        if (i =="c") {
                result5<-wilcox.test(synnoncpg,nonsynnoncpg, alternative = "less", paired = FALSE) 
                WilcoxTest.nt.sc$nt[5]<-i
                WilcoxTest.nt.sc$test[5]<-result5[[7]]
                WilcoxTest.nt.sc$P.value[5]<-result5[[3]]     }           
        if (i =="g") {
                result6<-wilcox.test(synnoncpg,nonsynnoncpg, alternative = "less", paired = FALSE) 
                WilcoxTest.nt.sc$nt[6]<-i
                WilcoxTest.nt.sc$test[6]<-result6[[7]]
                WilcoxTest.nt.sc$P.value[6]<-result6[[3]]    }
        }
     k=k+1
} 
write.csv(WilcoxTest.nt.sc,"Output1A/SelCoeff/WilcoxTestResults_byNT.csv")
write.csv(nuc.sc,"Output1A/SelCoeff/SC_Summary_byNT.csv")






