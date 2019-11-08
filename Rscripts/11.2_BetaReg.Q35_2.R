#GLM/Beta regression preparation:

library(tidyverse)
library(zoo)
library(purrr)
library(MASS)
library(betareg)
library(miscTools)
source("Rscripts/baseRscript.R")

####################

#Read data file
BetaReg   <-read.csv("Output1A/GLM/BetaRegFull.Ts.FilteredData.csv", row.names = 1, stringsAsFactors = F)

#1. Format the data:
### addd gene annotation info for all genes
genes<-read.csv("Data/HCV_annotations2.csv")
genenames<-genes$Gene
gene<-c()
for (i in 2:12){
        gene<-c(gene, rep(i, times=genes$end[i]-genes$start[i]+1))
}

n<-data.frame(pos=342:(length(gene)+341))
g<-cbind(n,gene)


colnames(BetaReg)[1]<-"pos"
betar<-merge(BetaReg, g, by ="pos")

for (i in 2:12){
        gname<-paste(genes$Gene[i])
        n<-ncol(betar)
        betar[,n+1]<-0
        colnames(betar)[n+1]<-gname
        betar[betar$gene==i,n+1]<-1
}
co<-which(colnames(betar)=="NS1(P7)" )
colnames(betar)[co]<-"NS1"

write.csv(betar,paste0("Output1A/GLM/BetaRegFull.Ts.FilteredData.csv"))

###

betar<-read.csv("Output1A/GLM/BetaRegFull.Ts.FilteredData.csv", stringsAsFactors = F, row.names = 1)

mod.g1 <- betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                          Core +E1 +HVR1++E2 +NS1 +NS2+NS3+NS4A+NS5A+NS5B, data = betar[betar$Stop == 0,])
AIC(mod.g1) #-75866.15
#remove the ns genes: NS3, and NS5A        
mod.g2 <- update(mod.g1, ~. -NS3 - NS5A)
AIC(mod.g2) # -75870.06  *** BEST MODEL
summary(mod.g2)

m3<-betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                    Core +E1 +HVR1++E2 +NS1 +NS2+NS4A+NS5B+CpG:Nonsyn, data = betar[betar$Stop == 0,])

AIC(m3) #-75868.22


m4<-betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                    Core +E1 +HVR1++E2 +NS1 +NS2+NS4A+NS5B+, data = betar[betar$Stop == 0,])
AIC(m4)





