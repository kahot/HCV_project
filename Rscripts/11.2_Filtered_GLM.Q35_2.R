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
GlmData<-read.csv("Output/GLM/GlmData_Ts.csv")
GlmData.na<-read.csv("Output/GLM/GlmData_Ts_NA.csv")
GlmData.0<-read.csv("Output/GLM/GlmData_Ts_zero.csv")

GlmData<-GlmData[,-1]
GlmData.na<-GlmData.na[,-1]
GlmData.0<-GlmData.0[,-1]

glmdata<-list()
glmdata[[1]]<-GlmData
glmdata[[2]]<-GlmData.na
glmdata[[3]]<-GlmData.0
names(glmdata)[1]<-"Ts"
names(glmdata)[2]<-"Ts_NA"
names(glmdata)[2]<-"Ts_zero"

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

glmdata2<-list()
for (k in 1:length(glmdata)){
        GlmData1<-glmdata[[k]]
        colnames(GlmData1)[1]<-"pos"
        glmData<-merge(GlmData1, g, by ="pos")
        
        for (i in 2:12){
                gname<-paste(genes$Gene[i])
                n<-ncol(glmData)
                glmData[,n+1]<-0
                colnames(glmData)[n+1]<-gname
                glmData[glmData$gene==i,n+1]<-1
        }

        co<-which(colnames(glmData)=="NS1(P7)" )
        colnames(glmData)[co]<-"NS1"
        
        glmdata2[[k]]<-glmData
        names(glmdata2)[k]<-names(glmdata)[k]
        #write.csv(glmData,paste0("Output/GLM/GlmdataFull.",names(glmdata)[k],".Q35.csv"))
}

################
#2. run GLM for each dataset

# 2.1  Ts.Q35
glmData<-glmdata2[[1]]
# not including NS4B because betareg won't run due to Rank deficiency
mod.g1 <- betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                          Core+E1 +HVR1 + E2 +NS1 +NS2+NS3 +NS4A + NS5A+NS5B, data = glmData[glmData$Stop == 0,])
summary(mod.g1)
#nofilter Q35
#             Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.4737177  0.0229838 -194.647  < 2e-16 ***
#t            0.1157091  0.0208632    5.546 2.92e-08 ***
#c           -0.5323171  0.0205090  -25.955  < 2e-16 ***
#g           -0.7220995  0.0229500  -31.464  < 2e-16 ***
#CpG         -0.0725060  0.0134435   -5.393 6.91e-08 ***
#Nonsyn      -0.7585679  0.0223505  -33.940  < 2e-16 ***
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
#t:Nonsyn    -0.1860138  0.0268162   -6.937 4.02e-12 ***
#c:Nonsyn    -0.1431568  0.0277536   -5.158 2.49e-07 ***
#g:Nonsyn     0.0929057  0.0286821    3.239   0.0012 **   
#remove the ns genes: NS3, and NS5A        

mod.g2 <- update(mod.g1, ~. -NS3 )
AIC(mod.g2)
mod.g3 <- update(mod.g2, ~. -NS4A)
mod.g4 <- update(mod.g3, ~. -NS2)

AIC(mod.g1,mod.g2, mod.g3)

mod.g2 <- update(mod.g1, ~. -NS3 - NS5A)
summary(mod.g3)
#(Intercept) -4.47210    0.01836 -243.547  < 2e-16 ***
#t            0.11574    0.02086    5.548 2.89e-08 ***
#c           -0.53233    0.02050  -25.962  < 2e-16 ***
#g           -0.72216    0.02295  -31.467  < 2e-16 ***
#CpG         -0.07246    0.01343   -5.395 6.86e-08 ***
#Nonsyn      -0.75852    0.02235  -33.941  < 2e-16 ***
#bigAAChange -0.12491    0.01404   -8.894  < 2e-16 ***
#Core        -0.21688    0.01979  -10.960  < 2e-16 ***
#E1           0.03746    0.01816    2.063  0.03916 *  
#HVR1         0.23570    0.04633    5.088 3.63e-07 ***
#E2           0.08428    0.01419    5.940 2.85e-09 ***
#NS1          0.06986    0.03006    2.324  0.02013 *  
#NS2          0.08954    0.01700    5.266 1.40e-07 ***
#NS4A        -0.05591    0.03467   -1.612  0.10689    
#NS5B        -0.16317    0.01556  -10.484  < 2e-16 ***
#t:Nonsyn    -0.18612    0.02681   -6.942 3.87e-12 ***
#c:Nonsyn    -0.14322    0.02775   -5.161 2.46e-07 ***
#g:Nonsyn     0.09295    0.02868    3.241  0.00119 ** 


#remove the ns genes: NS3, NS4A, and MS5A
mod.g2.2 <- update(mod.g1, ~. -NS3 -NS4A - NS5A)
summary(mod.g2.2)
#            Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.47368    0.01834 -243.971  < 2e-16 ***
#t            0.11525    0.02086    5.524 3.32e-08 ***
#c           -0.53291    0.02051  -25.987  < 2e-16 ***
#g           -0.72290    0.02295  -31.498  < 2e-16 ***
#CpG         -0.07259    0.01343   -5.403 6.54e-08 ***
#Nonsyn      -0.75895    0.02235  -33.958  < 2e-16 ***
#bigAAChange -0.12471    0.01404   -8.879  < 2e-16 ***
#Core        -0.21476    0.01975  -10.873  < 2e-16 ***
#E1           0.03963    0.01812    2.187  0.02878 *  
#HVR1         0.23781    0.04632    5.134 2.83e-07 ***
#E2           0.08641    0.01414    6.113 9.78e-10 ***
#NS1          0.07202    0.03004    2.398  0.01649 *  
#NS2          0.09168    0.01696    5.405 6.47e-08 ***
#NS5B        -0.16105    0.01552  -10.380  < 2e-16 ***
#t:Nonsyn    -0.18616    0.02682   -6.942 3.86e-12 ***
#c:Nonsyn    -0.14228    0.02775   -5.127 2.94e-07 ***
#g:Nonsyn     0.09324    0.02868    3.251  0.00115 ** 
        

AIC(mod.g1,mod.g2,mod.g2.2)
#         df       AIC
#mod.g1   21 -75866.15
#mod.g2   19 -75870.06
#mod.g2.2 18 -75869.38

sum.g1<-summary(mod.g1)
sum.g2<-summary(mod.g2)

modcoef1<-sum.g1$coefficients
coef1<-modcoef1[[1]]
write.csv(coef1,"Output/GLM/BetaReg_mod.g1_Ts.Q35.csv")
modcoef2<-sum.g2$coefficients
coef2<-modcoef2[[1]]
write.csv(coef2,"Output/GLM/BetaReg_mod.g2_Ts.Q35.csv")


#######################################################################################################

#2.2. NA replaced (mf<0.001)
GlmData1<-glmdata2[[2]]

mod.g1 <- betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                          Core +E1 +HVR1++E2 +NS1 +NS2++NS3 +NS4A+NS5A+NS5B, data = glmData[glmData$Stop == 0,])
summary(mod.g1)

#Ts_NA Q35
#             Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.4737177  0.0229838 -194.647  < 2e-16 ***
#t            0.1157091  0.0208632    5.546 2.92e-08 ***
#c           -0.5323171  0.0205090  -25.955  < 2e-16 ***
#g           -0.7220995  0.0229500  -31.464  < 2e-16 ***
#CpG         -0.0725060  0.0134435   -5.393 6.91e-08 ***
#Nonsyn      -0.7585679  0.0223505  -33.940  < 2e-16 ***
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
#t:Nonsyn    -0.1860138  0.0268162   -6.937 4.02e-12 ***
#c:Nonsyn    -0.1431568  0.0277536   -5.158 2.49e-07 ***
#g:Nonsyn     0.0929057  0.0286821    3.239   0.0012 ** 
        

#remove the two ns genes in Ts-no.filter
mod.g2 <- update(mod.g1, ~. -NS3 - NS5A)
summary(mod.g2)
#            Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.47210    0.01836 -243.547  < 2e-16 ***
t            0.11574    0.02086    5.548 2.89e-08 ***
#c           -0.53233    0.02050  -25.962  < 2e-16 ***
#g           -0.72216    0.02295  -31.467  < 2e-16 ***
#CpG         -0.07246    0.01343   -5.395 6.86e-08 ***
#Nonsyn      -0.75852    0.02235  -33.941  < 2e-16 ***
#bigAAChange -0.12491    0.01404   -8.894  < 2e-16 ***
#Core        -0.21688    0.01979  -10.960  < 2e-16 ***
#E1           0.03746    0.01816    2.063  0.03916 *  
#HVR1         0.23570    0.04633    5.088 3.63e-07 ***
#E2           0.08428    0.01419    5.940 2.85e-09 ***
#NS1          0.06986    0.03006    2.324  0.02013 *  
#NS2          0.08954    0.01700    5.266 1.40e-07 ***
#NS4A        -0.05591    0.03467   -1.612  0.10689    
#NS5B        -0.16317    0.01556  -10.484  < 2e-16 ***
#t:Nonsyn    -0.18612    0.02681   -6.942 3.87e-12 ***
#c:Nonsyn    -0.14322    0.02775   -5.161 2.46e-07 ***
#g:Nonsyn     0.09295    0.02868    3.241  0.00119 ** 
        
mod.g3 <- update(mod.g2, ~. -NS4A)
summary(mod.g3)        

#            Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.47368    0.01834 -243.971  < 2e-16 ***
#t            0.11525    0.02086    5.524 3.32e-08 ***
#c           -0.53291    0.02051  -25.987  < 2e-16 ***
#g           -0.72290    0.02295  -31.498  < 2e-16 ***
#CpG         -0.07259    0.01343   -5.403 6.54e-08 ***
#Nonsyn      -0.75895    0.02235  -33.958  < 2e-16 ***
#bigAAChange -0.12471    0.01404   -8.879  < 2e-16 ***
#Core        -0.21476    0.01975  -10.873  < 2e-16 ***
#E1           0.03963    0.01812    2.187  0.02878 *  
#HVR1         0.23781    0.04632    5.134 2.83e-07 ***
#E2           0.08641    0.01414    6.113 9.78e-10 ***
#NS1          0.07202    0.03004    2.398  0.01649 *  
#NS2          0.09168    0.01696    5.405 6.47e-08 ***
#NS5B        -0.16105    0.01552  -10.380  < 2e-16 ***
#t:Nonsyn    -0.18616    0.02682   -6.942 3.86e-12 ***
#c:Nonsyn    -0.14228    0.02775   -5.127 2.94e-07 ***
#g:Nonsyn     0.09324    0.02868    3.251  0.00115 ** 
        
        
        
AIC(mod.g1,mod.g2,mod.g3)
#       df       AIC
#mod.g1 21 -75866.15
#mod.g2 19 -75870.06
#mod.g3 18 -75869.38
#mod.g4 17 -75280.82

sum.g1<-summary(mod.g1)
sum.g2<-summary(mod.g2)

modcoef1<-sum.g1$coefficients
coef1<-modcoef1[[1]]
write.csv(coef1,"Output/GLM/BetaReg_mod.g1_Ts_NA.Q35.csv")
modcoef2<-sum.g2$coefficients
coef2<-modcoef2[[1]]
write.csv(coef2,"Output/GLM/BetaReg_mod.g2_Ts_NA.Q35.csv")



#######################################################################################################
#######################################################################################################


#3. zero replaced 
glmData<-glmdata2[[3]]

################

mod.g1 <- betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                          Core +E1 +HVR1++E2 +NS1 +NS2++NS3 +NS4A+NS5A+NS5B, data = glmData[glmData$Stop == 0,])
summary(mod.g1)

#             Estimate Std. Error  z value Pr(>|z|) 
#(Intercept) -4.471949   0.023959 -186.653  < 2e-16 ***
#t            0.115095   0.021692    5.306 1.12e-07 ***
#c           -0.547663   0.021355  -25.646  < 2e-16 ***
#g           -0.741804   0.023931  -30.998  < 2e-16 ***
#CpG         -0.071258   0.013957   -5.105 3.30e-07 ***
#Nonsyn      -0.756970   0.023252  -32.555  < 2e-16 ***
#bigAAChange -0.128630   0.014726   -8.735  < 2e-16 ***
#Core        -0.217036   0.025211   -8.609  < 2e-16 ***
#E1           0.045147   0.023773    1.899   0.0576 .  
#HVR1         0.238119   0.050596    4.706 2.52e-06 ***
#E2           0.091418   0.020626    4.432 9.33e-06 ***
#NS1          0.077258   0.034495    2.240   0.0251 *  
#NS2          0.096171   0.022828    4.213 2.52e-05 ***
#NS3          0.003951   0.018760    0.211   0.8332    
#NS4A        -0.053520   0.038971   -1.373   0.1696    
#NS5A        -0.002359   0.019988   -0.118   0.9060    
#NS5B        -0.169019   0.021716   -7.783 7.07e-15 ***
#t:Nonsyn    -0.185829   0.027857   -6.671 2.55e-11 ***
#c:Nonsyn    -0.207779   0.029118   -7.136 9.63e-13 ***
#g:Nonsyn     0.045865   0.030020    1.528   0.1266


#remove the two ns genes
mod.g2 <- update(mod.g1, ~. -NS3 - NS5A)
summary(mod.g2)
#            Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.47086    0.01910 -234.079  < 2e-16 ***
#t            0.11512    0.02169    5.308 1.11e-07 ***
#c           -0.54765    0.02135  -25.651  < 2e-16 ***
#g           -0.74189    0.02393  -31.001  < 2e-16 ***
#CpG         -0.07114    0.01395   -5.102 3.37e-07 ***
#Nonsyn      -0.75686    0.02325  -32.554  < 2e-16 ***
#bigAAChange -0.12869    0.01472   -8.745  < 2e-16 ***
#Core        -0.21815    0.02070  -10.539  < 2e-16 ***
#E1           0.04402    0.01896    2.322   0.0202 *  
#HVR1         0.23698    0.04854    4.882 1.05e-06 ***
#E2           0.09029    0.01481    6.095 1.09e-09 ***
#NS1          0.07614    0.03138    2.427   0.0152 *  
#NS2          0.09505    0.01776    5.353 8.65e-08 ***
#NS4A        -0.05465    0.03625   -1.508   0.1316    
#NS5B        -0.17015    0.01631  -10.435  < 2e-16 ***
#t:Nonsyn    -0.18594    0.02785   -6.676 2.46e-11 ***
#c:Nonsyn    -0.20788    0.02911   -7.140 9.33e-13 ***
#g:Nonsyn     0.04597    0.03002    1.531   0.1256 


mod.g3 <- update(mod.g2, ~. -NS4A)
summary(mod.g3)        

mod.g4 <- update(mod.g3, ~. -g:Nonsyn)
summary(mod.g4)        


AIC(mod.g1,mod.g2,mod.g3, mod.g4)
#       df       AIC
#mod.g1 21 -75616.03
#mod.g2 19 -75619.86
#mod.g3 18 -75619.52
#mod.g4 17 -75619.15

sum.g1<-summary(mod.g1)
sum.g2<-summary(mod.g2)

modcoef1<-sum.g1$coefficients
coef1<-modcoef1[[1]]
write.csv(coef1,"Output/GLM/BetaReg_mod.g1_Ts_zero.Q35.csv")
modcoef2<-sum.g2$coefficients
coef2<-modcoef2[[1]]
write.csv(coef2,"Output/GLM/BetaReg_mod.g2_Ts_zero.Q35.csv")



