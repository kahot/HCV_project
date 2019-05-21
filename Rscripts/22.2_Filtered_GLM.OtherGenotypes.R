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
glmdata<-list()
fnames<-c("1B","3A")
for (i in 1:length(fnames)){
        filename<-fnames[i]
        df<-read.csv(paste0("Output_all/GlmData_",filename,".csv"),stringsAsFactors = F)
        df<-df[,-1]
        glmdata[[i]]<-df
        names(glmdata)[i]<-filename
}



# Format the data:
### addd gene annotation info for all genes
genes<-read.csv("Data/HCV_annotations_joined.csv")
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
        write.csv(glmData,paste0("Output_all/Glmdata_withGenes.",names(glmdata)[k],".csv"))
}

################
#2. run GLM for each dataset

# 2.1  Genotype 1B
glmData<-glmdata2[[1]]
glmData2<-glmData[complete.cases(glmData),] #7882 for 1B

# check the AIC for using 'gene' as opposed to breaking down the genes
mod0<-betareg(mean ~ t + c + g + CpG + Nonsyn +bigAAChange+ t:Nonsyn + c:Nonsyn + g:Nonsyn +gene,data = glmData2[glmData2$Stop == 0,])
summary(mod0)
AIC(mod0) #-72234.93

mod1<-betareg(mean ~ t + c + g + CpG + Nonsyn +bigAAChange+ t:Nonsyn + c:Nonsyn + g:Nonsyn +
              Core +E1 +HVR1 +E2 +NS1 +NS2++NS3 +NS4A+NS5A+NS5B,data = glmData2[glmData2$Stop == 0,])
summary(mod1)
#Coefficients (mean model with logit link):
#Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.514164   0.028822 -156.620  < 2e-16 ***
#t            0.099787   0.026487    3.767 0.000165 ***
#c           -0.503271   0.025815  -19.495  < 2e-16 ***
#g           -0.681773   0.028354  -24.045  < 2e-16 ***
#CpG         -0.094342   0.016499   -5.718 1.08e-08 ***
#Nonsyn      -0.686016   0.027408  -25.030  < 2e-16 ***
#bigAAChange -0.119298   0.016689   -7.148 8.79e-13 ***
#Core        -0.173673   0.028882   -6.013 1.82e-09 ***
#E1          -0.012440   0.027994   -0.444 0.656762    
#HVR1         0.150378   0.073308    2.051 0.040235 *  
#E2           0.097604   0.023908    4.082 4.46e-05 ***
#NS1          0.120803   0.039817    3.034 0.002414 ** 
#NS2          0.126081   0.026240    4.805 1.55e-06 ***
#NS3          0.014452   0.021669    0.667 0.504809    
#NS4A        -0.033898   0.044191   -0.767 0.443034    
#NS5A         0.001549   0.023108    0.067 0.946545    
#NS5B        -0.134684   0.024863   -5.417 6.06e-08 ***
#t:Nonsyn    -0.171492   0.033242   -5.159 2.48e-07 ***
#c:Nonsyn    -0.152125   0.033513   -4.539 5.65e-06 ***
#g:Nonsyn     0.042478   0.034696    1.224 0.220836 

#remove the two genes:  
mod2 <- update(mod1, ~. -NS3 - NS5A)
summary(mod2)
#Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.50662    0.02305 -195.526  < 2e-16 ***
#t            0.09989    0.02648    3.772 0.000162 ***
#c           -0.50347    0.02580  -19.516  < 2e-16 ***
#g           -0.68213    0.02835  -24.064  < 2e-16 ***
#CpG         -0.09448    0.01649   -5.729 1.01e-08 ***
#Nonsyn      -0.68586    0.02739  -25.037  < 2e-16 ***
#bigAAChange -0.11929    0.01668   -7.152 8.58e-13 ***
#Core        -0.18114    0.02366   -7.657 1.91e-14 ***
#E1          -0.01991    0.02260   -0.881 0.378461    
#HVR1         0.14296    0.07142    2.002 0.045307 *  
#E2           0.09014    0.01726    5.221 1.78e-07 ***
#NS1          0.11336    0.03623    3.129 0.001754 ** 
#NS2          0.11861    0.02036    5.825 5.71e-09 ***
#NS4A        -0.04136    0.04098   -1.009 0.312876    
#NS5B        -0.14215    0.01853   -7.669 1.73e-14 ***
#t:Nonsyn    -0.17174    0.03324   -5.167 2.38e-07 ***
#c:Nonsyn    -0.15228    0.03350   -4.546 5.48e-06 ***
#g:Nonsyn     0.04275    0.03469    1.232 0.217894 


#remove other ns genes
mod3 <- update(mod2, ~. -E1 -NS4A)
summary(mod3)
#Coefficients (mean model with logit link):
#Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.51070    0.02284 -197.478  < 2e-16 ***
#t            0.09983    0.02648    3.770 0.000163 ***
#c           -0.50323    0.02580  -19.504  < 2e-16 ***
#g           -0.68187    0.02835  -24.051  < 2e-16 ***
#CpG         -0.09444    0.01649   -5.726 1.03e-08 ***
#Nonsyn      -0.68555    0.02740  -25.024  < 2e-16 ***
#bigAAChange -0.11882    0.01667   -7.128 1.02e-12 ***
#Core        -0.17723    0.02345   -7.558 4.10e-14 ***
#HVR1         0.14677    0.07136    2.057 0.039706 *  
#E2           0.09401    0.01698    5.537 3.08e-08 ***
#NS1          0.11728    0.03609    3.249 0.001157 ** 
#NS2          0.12252    0.02011    6.092 1.12e-09 ***
#NS5B        -0.13825    0.01827   -7.568 3.78e-14 ***
#t:Nonsyn    -0.17247    0.03324   -5.188 2.12e-07 ***
#c:Nonsyn    -0.15233    0.03351   -4.546 5.46e-06 ***
#g:Nonsyn     0.04179    0.03470    1.204 0.22841


mod4 <- update(mod3, ~. -g:Nonsyn)
summary(mod4)        

AIC(mod1,mod2,mod3,mod4)
#     df       AIC
#mod1 21 -72442.54
#mod2 19 -72444.91
#mod3 17 -72448.12
#mod4 16 -72448.68

sum.g1<-summary(mod1)
sum.g2<-summary(mod2)
sum.g3<-summary(mod3)

modcoef1<-sum.g1$coefficients
coef1<-modcoef1[[1]]
write.csv(coef1,"Output_all/GLM/BetaReg_mod1_1B.csv")
modcoef2<-sum.g2$coefficients
coef2<-modcoef2[[1]]
write.csv(coef2,"Output_all/GLM/BetaReg_mod2_1B.csv")
modcoef3<-sum.g3$coefficients
coef3<-modcoef3[[1]]
write.csv(coef3,"Output_all/GLM/BetaReg_mod3_1B.csv")

#######################################################################################################

#2.2.3A
glmData<-glmdata2[[2]]
glmData2<-glmData[complete.cases(glmData),] #7802 for 3A

mod.g1 <- betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                          Core +E1 +HVR1++E2 +NS1 +NS2++NS3 +NS4A+NS5A+NS5B, data = glmData2[glmData2$Stop == 0,])
summary(mod.g1)


#            Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.86625    0.02639 -184.367  < 2e-16 ***
#t            0.04428    0.02447    1.810 0.070303 .  
#c           -0.36971    0.02472  -14.956  < 2e-16 ***
#g           -0.51226    0.02770  -18.496  < 2e-16 ***
#CpG         -0.02761    0.01538   -1.795 0.072599 .  
#Nonsyn      -0.48756    0.02568  -18.988  < 2e-16 ***
#bigAAChange -0.05449    0.01615   -3.374 0.000740 ***
#Core        -0.06431    0.02803   -2.294 0.021787 *  
#E1           0.04936    0.02712    1.820 0.068810 .  
#HVR1         0.31422    0.06389    4.918 8.73e-07 ***
#E2           0.06925    0.02360    2.934 0.003346 ** 
#NS1          0.08021    0.04059    1.976 0.048125 *  
#NS2          0.06692    0.02617    2.558 0.010542 *  
#NS3          0.00631    0.02128    0.297 0.766845    
#NS4A         0.02086    0.04270    0.488 0.625199    
#NS5A         0.05991    0.02246    2.668 0.007641 ** 
#NS5B        -0.06717    0.02451   -2.741 0.006132 ** 
#t:Nonsyn    -0.10529    0.03105   -3.392 0.000695 ***
#c:Nonsyn    -0.29227    0.03252   -8.987  < 2e-16 ***
#g:Nonsyn    -0.12079    0.03394   -3.559 0.000372 ***
        

#remove the two genes 
mod.g2 <- update(mod.g1, ~. -NS3 -NS5A)
summary(mod.g2)
#            Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.50662    0.02305 -195.526  < 2e-16 ***
#t            0.043483   0.024477    1.777 0.075650 .  
#c           -0.370883   0.024731  -14.997  < 2e-16 ***
#g           -0.511055   0.027707  -18.445  < 2e-16 ***
#CpG         -0.027580   0.015389   -1.792 0.073106 .  
#Nonsyn      -0.489317   0.025686  -19.050  < 2e-16 ***
#bigAAChange -0.053205   0.016144   -3.296 0.000982 ***
#Core        -0.086959   0.022913   -3.795 0.000148 ***
#E1           0.026871   0.021822    1.231 0.218178    
#HVR1         0.291781   0.061876    4.716 2.41e-06 ***
#E2           0.046657   0.017222    2.709 0.006746 ** 
#NS1          0.057623   0.037264    1.546 0.122024    
#NS2          0.044379   0.020616    2.153 0.031342 *  
#NS4A        -0.001687   0.039559   -0.043 0.965982    
#NS5B        -0.089786   0.018446   -4.868 1.13e-06 ***
#t:Nonsyn    -0.104653   0.031058   -3.370 0.000753 ***
#c:Nonsyn    -0.289309   0.032533   -8.893  < 2e-16 ***
#g:Nonsyn    -0.121892   0.033955   -3.590 0.000331 ***   
        

mod.g2.2<- update(mod.g1, ~. -NS3 -NS4A)
summary(mod.g2.2)


mod.g3 <- update(mod.g2, ~. -E1 -NS4A)
summary(mod.g3)        


#            Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.83984    0.02062 -234.702  < 2e-16 ***
#t            0.04403    0.02447    1.799 0.071975 .  
#c           -0.37066    0.02473  -14.991  < 2e-16 ***
#g           -0.51128    0.02771  -18.453  < 2e-16 ***
#CpG         -0.02780    0.01539   -1.806 0.070866 .  
#Nonsyn      -0.48898    0.02568  -19.042  < 2e-16 ***
#bigAAChange -0.05372    0.01614   -3.329 0.000873 ***
#Core        -0.09017    0.02270   -3.972 7.12e-05 ***
#HVR1         0.28846    0.06180    4.668 3.04e-06 ***
#E2           0.04343    0.01693    2.565 0.010322 *  
#NS1          0.05436    0.03713    1.464 0.143184    
#NS2          0.04113    0.02037    2.019 0.043506 *  
#NS5B        -0.09298    0.01818   -5.115 3.14e-07 ***
#t:Nonsyn    -0.10469    0.03105   -3.372 0.000747 ***
#c:Nonsyn    -0.28987    0.03253   -8.912  < 2e-16 ***
#g:Nonsyn    -0.12149    0.03395   -3.578 0.000346 ***        
        
        
AIC(mod.g1,mod.g2,mod.g2.2, mod.g3)
#       df       AIC
#mod.g1   21 -73624.43
#mod.g2   19 -73617.28
#mod.g2.2 19 -73628.17
#mod.g3   17 -73619.76

sum.g1<-summary(mod.g1)
sum.g2<-summary(mod.g2)
sum.g22<-summary(mod.g2.2)

modcoef1<-sum.g1$coefficients
coef1<-modcoef1[[1]]
write.csv(coef1,"Output_all/GLM/BetaReg_mod1_3A.csv")
modcoef2<-sum.g2$coefficients
coef2<-modcoef2[[1]]
write.csv(coef2,"Output_all/GLM/BetaReg_mod2_3A.csv")
modcoef22<-sum.g22$coefficients
coef22<-modcoef22[[1]]
write.csv(coef22,"Output_all/GLM/BetaReg_mod2.2_3A.csv")





###################
## plot the effects of each component for the 3 genotypes

effects1<-read.csv("Output_all/GLM/Mod1_effects_comparison1.csv", stringsAsFactors = F)
effects<-read.csv("Output_all/GLM/Mod1_effects_comparison.csv", stringsAsFactors = F)
effects$factor<-factor(effects$factor, levels=effects$factor[1:19])


cols2<-c("#66CCEE","#EE6677B3" ,"#228833B3")
p1<-ggplot(effects, aes(factor(factor),percent, fill =genotype)) +geom_bar(stat="identity", position="dodge")+
        scale_fill_manual(values=cols2)+
        theme_test() +
        theme(axis.text=element_text(size=13),
               axis.title=element_text(size=14,face="bold"))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        labs(x="", y="Estimated effects (%)")

ggsave("Output_all/GLM/Compare_effects_glm2.pdf", plot=p1, width = 14, height = 8)
plot(effects, col )