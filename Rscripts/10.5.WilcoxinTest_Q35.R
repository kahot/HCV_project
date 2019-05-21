library(tidyverse)
library(zoo)
library(purrr)

source("Rscripts/baseRscript.R")

filtered<-list()
fnames<-c("Ts", "Ts_NA", "Ts_zero")
for (i in 1:3){
        filename<-fnames[i]
        df<-read.csv(paste0("Output/Mut.freq.filtered/Summary_",filename,".Q35.csv"),stringsAsFactors = F)
        df<-df[,-1]
        filtered[[i]]<-df
        names(filtered)[i]<-filename
}

f<-filtered[[1]]
f2<-filtered[[2]]

mean(f$mean[f$Type=="syn"&f$ref=="a"], na.rm=T)
mean(f2$mean[f2$Type=="syn"&f2$ref=="a"], na.rm=T)

##################
# 1) run Wilcoxin Test on means based on mutation types

for (j in 1:3){
        WilcoxTest.results.nt<-data.frame(matrix(ncol=3,nrow=6))
        colnames(WilcoxTest.results.nt)<-c("nt","test","P.value")
        WilcoxTest.results<-data.frame(matrix(ncol=2,nrow=4))
        colnames(WilcoxTest.results)<-c("test","P.value")
        
        dat<-filtered[[j]]
        dat<-dat[dat$pos>=342, ] 
        ty<-which(colnames(dat)=="Type")
        fname=names(filtered)[j]
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
        rownames(se)<-c("syn_se","nonsyn_se","stop_se")
        rownames(m_nonCpG)<-c("syn_noncpg","nonsyn_noncpg","stop_noncpg")
        rownames(se_nonCpG)<-c("syn_noncpg_se","nonsyn_noncpg_se","stop_noncpg_se")
        rownames(m_CpG)<-c("syn_cpg","nonsyn_cpg","stop_cpg")
        rownames(se_CpG)<-c("syn_cpg_se","nonsyn_cpg_se","stop_cpg_se")
        
        MFbyType<-rbind(m,se,m_nonCpG,se_nonCpG,m_CpG,se_CpG)
        write.csv(MFbyType,paste0("Output/SummaryStats/MutFreq_byNt_byType_",fname,".Q35.csv"))

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
        write.csv(WilcoxTest.results.nt,paste0("Output/SummaryStats/WilcoxTestResults_",fname,".csv"))

        #general Test
        ## Compare Syn vs Nonsyn, Nonsyn vs. Stop
        R1<-wilcox.test(dat$mean[dat$Type=="syn"], dat$mean[dat$Type=="nonsyn"], alternative = "greater", paired = FALSE)  
        R2<-wilcox.test(dat$mean[dat$Type=="nonsyn"], dat$mean[dat$Type=="stop"], alternative = "greater", paired = FALSE) 
        
        ## Compare CpG making vs. Non-CpGmaking 
        T2<-dat[dat$ref=="a"|dat$ref=="t",]
        R3<-wilcox.test(T2$mean[T2$Type=="syn"&T2$makesCpG==0], T2$mean[T2$Type=="syn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 
        R4<-wilcox.test(T2$mean[T2$Type=="nonsyn"&T2$makesCpG==0], T2$mean[T2$Type=="nonsyn"&T2$makesCpG==1], alternative = "greater", paired = FALSE) 
        
        for (r in 1:4){
                result<-get(paste0('R',r))
                WilcoxTest.results$test[r]<-result[[7]]
                WilcoxTest.results$P.value[r]<-result[[3]]
        }
        write.csv(WilcoxTest.results, paste0("Output/SummaryStats/WilcoxTest_general.",fname,".csv"))
}

