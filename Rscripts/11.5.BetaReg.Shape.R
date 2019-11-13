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
BetaReg   <-read.csv("Output1A/GLM/BetaReg.Data.Shape.csv", row.names = 1, stringsAsFactors = F)

#1. Format the data:
### add gene annotation info for all genes
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

#write.csv(betar,paste0("Output1A/GLM/BetaRegFull.Ts.FilteredData.csv"))



###

#betar<-read.csv("Output1A/GLM/BetaRegFull.Ts.FilteredData.csv", stringsAsFactors = F, row.names = 1)

################
#2. run Beta Regression
mod1 <- betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                        Core +E1 +HVR1++E2 +NS1 +NS2+NS3+ NS4A+NS5A+NS5B + Shape, data = betar[betar$Stop == 0,])
AIC(mod1) #-75864.86
summary(mod1)
#Coefficients (mean model with logit link):
#Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.473459   0.022983 -194.643  < 2e-16 ***
#t            0.115463   0.020862    5.535 3.12e-08 ***
#c           -0.532611   0.020508  -25.970  < 2e-16 ***
#g           -0.722394   0.022949  -31.478  < 2e-16 ***
#CpG         -0.072600   0.013443   -5.401 6.64e-08 ***
#Nonsyn      -0.758651   0.022349  -33.945  < 2e-16 ***
#bigAAChange -0.124901   0.014054   -8.887  < 2e-16 ***
#Core        -0.223432   0.025866   -8.638  < 2e-16 ***
#E1           0.037405   0.022843    1.637  0.10154    
#HVR1         0.237318   0.048289    4.914 8.90e-07 ***
#E2           0.084702   0.019787    4.281 1.86e-05 ***
#NS1          0.071473   0.033029    2.164  0.03047 *  
#NS2          0.091152   0.021842    4.173 3.00e-05 ***
#NS3          0.002000   0.018022    0.111  0.91164    
#NS4A        -0.054267   0.037271   -1.456  0.14538    
#NS5A        -0.002029   0.019209   -0.106  0.91590    
#NS5B        -0.164191   0.020978   -7.827 5.00e-15 ***
#Shape        0.012517   0.014699    0.852  0.39444    
#t:Nonsyn    -0.185900   0.026814   -6.933 4.12e-12 ***
#c:Nonsyn    -0.143042   0.027752   -5.154 2.55e-07 ***
#g:Nonsyn     0.092746   0.028682    3.234  0.00122 ** 

#remove the ns genes: NS3, and NS5A        
mod2 <- update(mod1, ~. -NS3 - NS5A)
AIC(mod2) # -75868.79  (lower than without Shape variable)
summary(mod2)
#Coefficients (mean model with logit link):
#Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.47319    0.01842 -242.886  < 2e-16 ***
#t            0.11547    0.02086    5.536 3.10e-08 ***
#c           -0.53259    0.02050  -25.975  < 2e-16 ***
#g           -0.72244    0.02295  -31.480  < 2e-16 ***
#CpG         -0.07251    0.01343   -5.399 6.71e-08 ***
#Nonsyn      -0.75857    0.02235  -33.945  < 2e-16 ***
#bigAAChange -0.12495    0.01404   -8.897  < 2e-16 ***
#Core        -0.22371    0.02126  -10.524  < 2e-16 ***
#E1           0.03711    0.01817    2.042  0.04110 *  
#HVR1         0.23701    0.04635    5.113 3.16e-07 ***
#E2           0.08441    0.01419    5.949 2.69e-09 ***
#NS1          0.07119    0.03010    2.365  0.01803 *  
#NS2          0.09086    0.01708    5.321 1.03e-07 ***
#NS4A        -0.05457    0.03471   -1.572  0.11591    
#NS5B        -0.16449    0.01564  -10.515  < 2e-16 ***
#Shape        0.01251    0.01460    0.857  0.39152    
#t:Nonsyn    -0.18595    0.02681   -6.936 4.03e-12 ***
#c:Nonsyn    -0.14311    0.02775   -5.158 2.50e-07 ***
#g:Nonsyn     0.09282    0.02868    3.237  0.00121 ** 



## Mod1 minus NS3 only
#Coefficients (mean model with logit link):
#             Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.472011   0.019024 -235.066  < 2e-16 ***
#t            0.115482   0.020860    5.536 3.09e-08 ***
#c           -0.532649   0.020508  -25.973  < 2e-16 ***
#g           -0.722415   0.022949  -31.479  < 2e-16 ***
#CpG         -0.072632   0.013440   -5.404 6.52e-08 ***
#Nonsyn      -0.758678   0.022349  -33.947  < 2e-16 ***
#bigAAChange -0.124858   0.014049   -8.887  < 2e-16 ***
#Core        -0.224975   0.021861  -10.291  < 2e-16 ***
#E1           0.035960   0.018754    1.917  0.05519 .  
#HVR1         0.235898   0.046566    5.066 4.06e-07 ***
#E2           0.083262   0.014922    5.580 2.41e-08 ***
#NS1          0.070050   0.030439    2.301  0.02137 *  
#NS2          0.089728   0.017671    5.078 3.82e-07 ***
#NS4A        -0.055687   0.035000   -1.591  0.11160    
#NS5A        -0.003477   0.014083   -0.247  0.80500    
#NS5B        -0.165651   0.016339  -10.139  < 2e-16 ***
#Shape        0.012689   0.014623    0.868  0.38554    
#t:Nonsyn    -0.185956   0.026809   -6.936 4.03e-12 ***
#c:Nonsyn    -0.143043   0.027751   -5.154 2.54e-07 ***
#g:Nonsyn     0.092729   0.028681    3.233  0.00122 ** 


