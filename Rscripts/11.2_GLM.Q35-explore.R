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
GlmData   <-read.csv("Output1A/GLM/GlmData_Ts_FilteredData.csv", row.names = 1, stringsAsFactors = F)

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


colnames(GlmData)[1]<-"pos"
glmData<-merge(GlmData, g, by ="pos")

for (i in 2:12){
        gname<-paste(genes$Gene[i])
        n<-ncol(glmData)
        glmData[,n+1]<-0
        colnames(glmData)[n+1]<-gname
        glmData[glmData$gene==i,n+1]<-1
}
co<-which(colnames(glmData)=="NS1(P7)" )
colnames(glmData)[co]<-"NS1"

write.csv(glmData,paste0("Output1A/GLM/GlmdataFull.Ts.FilteredData.csv"))

glmData<-read.csv("Output1A/GLM/BetaRegFull.Ts.FilteredData.csv", stringsAsFactors = F, row.names = 1)
################
#2. run GLM
mod1 <- betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                        Core +E1 +HVR1++E2 +NS1 +NS2+NS3 +NS4A+NS5B, data = glmData[glmData$Stop == 0,])
summary(mod1)
AIC(mod1) #-75868.15

#              Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.4737177  0.0229838 -194.647  < 2e-16 ***
#t1           0.1157091  0.0208632    5.546 2.92e-08 ***
#c1          -0.5323171  0.0205090  -25.955  < 2e-16 ***
#g1          -0.7220995  0.0229500  -31.464  < 2e-16 ***
#CpG         -0.0725060  0.0134435   -5.393 6.91e-08 ***
#Nonsyn1     -0.7585679  0.0223505  -33.940  < 2e-16 ***
#bigAAChange -0.1249080  0.0140548   -8.887  < 2e-16 ***
#Core        -0.2152590  0.0240976   -8.933  < 2e-16 ***
#E1           0.0390887  0.0227549    1.718   0.0858 .  
#HVR1         0.2373264  0.0482918    4.914 8.90e-07 ***
#E2           0.0859036  0.0197335    4.353 1.34e-05 ***
#NS1          0.0714742  0.0330309    2.164   0.0305 *  
#NS2          0.0911630  0.0218428    4.174 3.00e-05 ***
#NS3          0.0035818  0.0179292    0.200   0.8417    
#NS4A        -0.0542810  0.0372726   -1.456   0.1453    
#NS5A        -0.0003067  0.0190949   -0.016   0.9872    
#NS5B        -0.1615474  0.0207384   -7.790 6.71e-15 ***
#t1:Nonsyn1  -0.1860138  0.0268162   -6.937 4.02e-12 ***
#c1:Nonsyn1  -0.1431568  0.0277536   -5.158 2.49e-07 ***
#g1:Nonsyn1   0.0929057  0.0286821    3.239   0.0012 ** 

m2<-update(mod1,~. -NS3 )
AIC(m2) #-75870.06
summary(m2)


m3<-update(m2,~. -NS4A)
AIC(m3)#-75869.38
summary(m3)
m4<-update(m3,~. -NS4A)
AIC(m4) #-75869.38  #best values
summary(m4)

m5<-update(m4, ~. -NS1)
AIC(m5)
###
#Best Model
mod <- betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                        Core +E1 +HVR1++E2 +NS1 +NS2++NS4A+NS5B, data = glmData[glmData$Stop == 0,])
summary(mod)
AIC(mod)
#-75870.06

summary<-summary(mod)
modcoef<-summary$coefficients
coef1<-modcoef[[1]]
write.csv(coef1,"Output1A/GLM/BetaReg_BestModel.Q35.csv")

#### include more interaction

mod1 <- betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                        Core +E1 +HVR1++E2 +NS1 +NS2++NS3 +NS4A+NS5A+NS5B , data = glmData[glmData$Stop == 0,])
summary(mod1)
AIC(mod1)

m2<-update(mod1,~. -NS5A )
AIC(m2)
summary(m2)
m3<-update(m2,~. -NS3)
AIC(m3)
summary(m3)
m4<-update(m3,~. -NS4A)
AIC(m4) #-75869.38  #best values
summary(m4)

m5<-update(m4, ~. -NS1)
AIC(m5)
AIC(mod1,m2,m3,m4)


####include Stop

mod1s <- betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + Stop+
                         Core +E1 +HVR1++E2 +NS1 +NS2++NS3+ NS4A+NS5A+NS5B, data = glmData)
summary(mod1s)
AIC(mod1s)
#-78305.91

m2<-update(mod1s,~. -NS3 )
AIC(m2)
summary(m2)

m3<-update(m2,~. -NS5A)
AIC(m3)
summary(m3)
m4<-update(m3,~. -NS5B)
AIC(m4) #-75869.38  #best values
summary(m4)

### Best Model ###
mod <- betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + Stop+
                         Core +E1 +HVR1++E2 +NS1 +NS2+NS4A+NS5B, data = glmData)
summary(mod)

#Coefficients (mean model with logit link):
#             Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.47336    0.01823 -245.430  < 2e-16 ***
#t            0.11564    0.02072    5.582 2.37e-08 ***
#c           -0.53469    0.02032  -26.313  < 2e-16 ***
#g           -0.71884    0.02262  -31.783  < 2e-16 ***
#CpG         -0.07237    0.01334   -5.423 5.86e-08 ***
#Nonsyn      -0.75904    0.02220  -34.190  < 2e-16 ***
#bigAAChange -0.12516    0.01396   -8.964  < 2e-16 ***
#Stop        -0.86570    0.04146  -20.879  < 2e-16 ***
#Core        -0.21264    0.01949  -10.909  < 2e-16 ***
#E1           0.03663    0.01791    2.045  0.04089 *  
#HVR1         0.23281    0.04591    5.071 3.95e-07 ***
#E2           0.08541    0.01398    6.111 9.90e-10 ***
#NS1          0.07153    0.02958    2.418  0.01562 *  
#NS2          0.08990    0.01674    5.372 7.80e-08 ***
#NS4A        -0.04988    0.03420   -1.458  0.14473    
#NS5B        -0.16186    0.01538  -10.524  < 2e-16 ***
#t:Nonsyn    -0.18611    0.02663   -6.988 2.79e-12 ***
#c:Nonsyn    -0.14192    0.02754   -5.153 2.57e-07 ***
#g:Nonsyn     0.08861    0.02836    3.125  0.00178 ** 

AIC(mod)
# -78309.87

summary<-summary(mod)
modcoef1<-summary$coefficients
coef1<-modcoef1[[1]]
write.csv(coef1,"Output1A/GLM/BetaReg_BestModel.Stop.Q35.csv")


