#GLM/Beta regression preparation:

library(tidyverse)
library(zoo)
library(purrr)
library(MASS)
library(betareg)

source("Rscripts/baseRscript.R")

## Run 15.Filterout_lowfrew.R first ##
# Read the summary data files
filtered<-list()
fnames<-c("TsMutFreq", "Ts_NA0.01", "Ts_zero0.01")
for (i in 1:3){
        filename<-fnames[i]
        df<-read.csv(paste0("Output/Mut.freq.filtered/Summary_",filename,".csv"),stringsAsFactors = F)
        df<-df[,-1]
        filtered[[i]]<-df
        names(filtered)[i]<-filename
}


#### For preparing GLM data formatting

library(miscTools)
##########
# 1. Transition original (reads>1000)
nucord <- c("a", "t", "c", "g")
glmform<-as.matrix(data.frame(a="",t="",c="",g="", Syn ="",Nonsyn="",Stop=""))
nuc<-as.matrix(data.frame(a="",t="",c="",g=""))
for (k in 2:3){
        dat<-filtered[[k]]
        dat1<-dat[,c("pos","makesCpG","bigAAChange","mean")]
        
        for (i in 1:nrow(dat)){
                atcg <- c(0,0,0,0)
                atcg[which(nucord == dat[i,]$ref)] <- 1
                nuc<-insertRow(nuc,i,atcg)
                nonsyn <- as.numeric(regexpr("nonsyn",dat[i,]$Type) > 0)
                stop <- as.numeric(regexpr("stop",dat[i,]$Type) > 0)
                syn<-as.numeric(regexpr("^syn",dat[i,]$Type) > 0)
                new<-c(atcg,syn,nonsyn,stop)
                glmform<-insertRow(glmform,i,new)
        }
        GlmData<-cbind(dat1$pos,glmform[1:nrow(dat1),])
        GlmData<-cbind(GlmData,dat1[,2:4])
        colnames(GlmData)[9]<-"CpG"
        write.csv(GlmData, paste0("Output/GLM/GlmData_",names(filtered)[k],".csv"))
}
        


############################################

#Read data file
GlmData<-read.csv("Output/GLM/GlmData_TsMutFreq.csv")
GlmData.na<-read.csv("Output/GLM/GlmData_Ts_NA0.01.csv")
GlmData.0<-read.csv("Output/GLM/GlmData_Ts_zero0.01.csv")

GlmData<-GlmData[,-1]
GlmData.na<-GlmData.na[,-1]
GlmData.0<-GlmData.0[,-1]

### select the files you want to run the analysis

#1. Origina transition data (reads>1000)
GlmData1<-GlmData
#2. NA replaced (mf<0.001)
GlmData1<-GlmData.na
#3. zero replaced 
GlmData1<-GlmData.0


### addd gene annotation info
genes<-read.csv("Data/HCV_annotations2.csv")
genenames<-genes$Gene

gene<-c()
for (i in 2:12){
gene<-c(gene, rep(i, times=genes$end[i]-genes$start[i]+1))
}

n<-data.frame(pos=342:(length(gene)+341))
g<-cbind(n,gene)

colnames(GlmData1)[1]<-"pos"
GlmData1<-merge(GlmData1, g, by ="pos")

### Add a column to show whether HVR1 or not 
GlmData1$hvr<-0
GlmData1$hvr[GlmData1$gene==4]<-1
#write.csv(GlmData1, "Output/GLM/GLMData_Ts_Original_with_gene_info.csv")

GlmData1<-GlmData1[!is.na(GlmData1$mean),] #7969 obs

## for GlmData.0
y.transf.betareg <- function(y){
        n.obs <- sum(!is.na(y))
        (y * (n.obs - 1) + 0.5) / n.obs
}

mod1 <- betareg(y.transf.betareg (mean) ~ t + c + g + CpG + CpG:t + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn:CpG + Nonsyn + t:Nonsyn:CpG +bigAAChange +hvr,  data = GlmData1[GlmData1$Stop == 0,])
summary(mod1)

#Coefficients (mean model with logit link):
#              Estimate Std. Error  z value Pr(>|z|)    
#(Intercept)     -4.83868    0.03781 -127.975  < 2e-16 ***
#t             0.21666    0.04568    4.743 2.10e-06 ***
#c            -0.55270    0.04147  -13.327  < 2e-16 ***
#g            -0.87031    0.04521  -19.250  < 2e-16 ***
#CpG          -0.04365    0.05898   -0.740    0.459    
#Nonsyn       -1.57734    0.05083  -31.032  < 2e-16 ***
#bigAAChange  -0.23237    0.02658   -8.744  < 2e-16 ***
#hvr           0.46652    0.08910    5.236 1.64e-07 ***
#t:CpG        -0.01177    0.07248   -0.162    0.871    
#t:Nonsyn     -0.45221    0.06247   -7.239 4.54e-13 ***
#c:Nonsyn     -0.00151    0.05683   -0.027    0.979    
#g:Nonsyn      0.34648    0.05847    5.925 3.11e-09 ***
#CpG:Nonsyn   -0.07070    0.07693   -0.919    0.358    
#t:CpG:Nonsyn  0.01189    0.10251    0.116    0.908          

#remove CpG:t & t:Nonsyn:CpG
mod1.1<-update(mod1, ~.-CpG:t - t:Nonsyn:CpG)
mod1.2 <- update(mod1.1, ~. -Nonsyn:CpG)
mod1.3 <- update(mod1.2, ~. -t:Nonsyn -c:Nonsyn -g:Nonsyn ) 

summary(mod1.1)
#Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.835500   0.032267 -149.861  < 2e-16 ***
#t            0.211994   0.035464    5.978 2.26e-09 ***
#c           -0.555882   0.036508  -15.226  < 2e-16 ***
#g           -0.873486   0.040701  -21.461  < 2e-16 ***
#CpG         -0.051457   0.034277   -1.501    0.133    
#Nonsyn      -1.580546   0.045370  -34.837  < 2e-16 ***
#bigAAChange -0.232372   0.026570   -8.746  < 2e-16 ***
#hvr          0.466080   0.089061    5.233 1.67e-07 ***
#t:Nonsyn    -0.447502   0.049495   -9.041  < 2e-16 ***
#c:Nonsyn     0.001689   0.052082    0.032    0.974    
#g:Nonsyn     0.349683   0.053870    6.491 8.51e-11 ***
#CpG:Nonsyn  -0.062838   0.049810   -1.262    0.207       

summary(mod1.2)->summary1.2
#Coefficients (mean model with logit link):
#             -4.82346    0.03075 -156.856  < 2e-16 ***
#t            0.21101    0.03546    5.950 2.68e-09 ***
#c           -0.56781    0.03519  -16.136  < 2e-16 ***
#g           -0.88538    0.03952  -22.402  < 2e-16 ***
#CpG         -0.08142    0.02489   -3.271  0.00107 ** 
#Nonsyn      -1.60396    0.04158  -38.576  < 2e-16 ***
#bigAAChange -0.23192    0.02657   -8.728  < 2e-16 ***
#hvr          0.46346    0.08910    5.202 1.98e-07 ***
#t:Nonsyn    -0.44692    0.04950   -9.028  < 2e-16 ***
#c:Nonsyn     0.02498    0.04868    0.513  0.60784    
#g:Nonsyn     0.37293    0.05059    7.372 1.69e-13 ***    

#####

summary(mod1.3)

AIC(mod1,mod1.1,mod1.2,mod1.3)        
#       df       AIC
#mod1   15 -84021.10
#mod1.1 13 -84025.07
#mod1.2 12 -84025.48
#mod1.3  9 -83724.78

##setting up to plot the estiamted values ###
modcoef1<-summary1.2$coefficients
modcoef<-modcoef1[[1]]
coef.vals <- modcoef[,1]
coef.pSE.vals <- coef.vals + modcoef[,2]
coef.mSE.vals <- coef.vals - modcoef[,2]




library(lmtest)
lrtest(mod1.2)  #comparing the goodness of fit of two statistical models (against a Null hypo)
#Model 1: mean ~ t + c + g + CpG + Nonsyn + bigAAChange + t:Nonsyn + c:Nonsyn + g:Nonsyn
#Model 2: mean ~ 1
#  #Df LogLik  Df  Chisq Pr(>Chisq)    
#1  12  36134                          
#2   2  32810 -10 6648.4  < 2.2e-16 ***

library(emmeans)
joint_tests(mod1.2)
# model term   df1 df2 F.ratio p.value
#t              1 Inf   19.837  <.0001
#c              1 Inf 1649.397  <.0001
#g              1 Inf 2476.290  <.0001
#Nonsyn         1 Inf 1910.033  <.0001
#t:c            1 Inf   18.785  <.0001
#t:g            1 Inf   24.853  <.0001
#t:Nonsyn       1 Inf   44.208  <.0001
#c:g            1 Inf  573.733  <.0001
#c:Nonsyn       1 Inf  107.219  <.0001
#g:Nonsyn       1 Inf  477.632  <.0001
#t:c:g          1 Inf   23.074  <.0001
#t:c:Nonsyn     1 Inf   59.512  <.0001
#t:g:Nonsyn     1 Inf   52.686  <.0001
#c:g:Nonsyn     1 Inf   48.688  <.0001
#t:c:g:Nonsyn   1 Inf   68.160  <.0001



plot(fitted(mod1.2), residuals(mod1.2)) 
plot(mod1.2)

write.csv(modcoef,"Output/SummaryStats/BetaRegression_BestFit_Model1.2_Ts_zero.csv")




##################
### Add a column to show whether HVR1 or not 
GlmData1<-GlmData
#2. NA replaced (mf<0.001)
GlmData1<-GlmData.na
#3. zero replaced 
GlmData1<-GlmData.0

### addd gene annotation info for all genes
genes<-read.csv("Data/HCV_annotations2.csv")
genenames<-genes$Gene
gene<-c()
for (i in 2:12){
        gene<-c(gene, rep(i, times=genes$end[i]-genes$start[i]+1))
}

n<-data.frame(pos=342:(length(gene)+341))
g<-cbind(n,gene)

colnames(GlmData1)[1]<-"pos"
GlmData1<-merge(GlmData1, g, by ="pos")

glmData<-GlmData1

for (i in 2:12){
        gname<-paste(genes$Gene[i])
        n<-ncol(glmData)
        glmData[,n+1]<-0
        colnames(glmData)[n+1]<-gname
        glmData[glmData$gene==i,n+1]<-1
}

co<-which(colnames(glmData)=="NS1(P7)" )
colnames(glmData)[co]<-"NS1"

mod.g1 <- betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                        Core +E1 +HVR1++E2 +NS1 +NS2++NS3 +NS4A+NS5A+NS5B, data = glmData[glmData$Stop == 0,])
summary(mod.g1)
#             Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -3.620455   0.034514 -104.897  < 2e-16 ***
#t            0.063153   0.034243    1.844 0.065143 .  
#c            0.145872   0.031000    4.706 2.53e-06 ***
#g            0.148146   0.032472    4.562 5.06e-06 ***
#CpG          0.003717   0.018429    0.202 0.840178    
#Nonsyn      -0.253962   0.032635   -7.782 7.15e-15 ***
#bigAAChange -0.147208   0.016057   -9.168  < 2e-16 ***
#Core        -0.017849   0.029211   -0.611 0.541168    
#E1           0.138897   0.027911    4.976 6.48e-07 ***
#HVR1         0.469722   0.058964    7.966 1.64e-15 ***
#E2           0.162176   0.024854    6.525 6.79e-11 ***
#NS1          0.178615   0.039251    4.551 5.35e-06 ***
#NS2          0.137850   0.027529    5.007 5.52e-07 ***
#NS3          0.082232   0.023048    3.568 0.000360 ***
#NS4A         0.010778   0.046888    0.230 0.818202    
#NS5A         0.079391   0.024294    3.268 0.001084 ** 
#NS5B        -0.016841   0.026128   -0.645 0.519223    
#t:Nonsyn    -0.147189   0.040255   -3.656 0.000256 ***
#c:Nonsyn    -0.008708   0.037081   -0.235 0.814329    
#g:Nonsyn     0.041123   0.037561    1.095 0.273595   


mod.g2 <- update(mod.g1, ~. -NS3 - NS5A)
summary(mod.g2)
#            Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -3.553043   0.029299 -121.268  < 2e-16 ***
#t            0.064400   0.034275    1.879 0.060255 .  
#c            0.144905   0.031024    4.671 3.00e-06 ***
#g            0.147422   0.032504    4.535 5.75e-06 ***
#CpG          0.002469   0.018440    0.134 0.893503    
#Nonsyn      -0.255572   0.032664   -7.824 5.11e-15 ***
#bigAAChange -0.145290   0.016064   -9.044  < 2e-16 ***
#Core        -0.084264   0.022974   -3.668 0.000245 ***
#E1           0.072667   0.021338    3.406 0.000660 ***
#HVR1         0.403654   0.056211    7.181 6.92e-13 ***
#E2           0.095859   0.017125    5.598 2.17e-08 ***
#NS1          0.112417   0.034908    3.220 0.001280 ** 
#NS2          0.071590   0.020830    3.437 0.000589 ***
#NS4A        -0.055324   0.043345   -1.276 0.201825    
#NS5B        -0.083074   0.018931   -4.388 1.14e-05 ***
#t:Nonsyn    -0.148935   0.040292   -3.696 0.000219 ***
#c:Nonsyn    -0.008906   0.037110   -0.240 0.810350    
#g:Nonsyn     0.040870   0.037598    1.087 0.277031  


mod.g3 <- update(mod.g1, ~. -c:Nonsyn -g:Nonsyn)
summary(mod.g3)
mod.g4 <- update(mod.g3, ~. -NS4A -NS5B)
summary(mod.g4)


AIC(mod.g1,mod.g2, mod.g3, mod.g4)
#       df       AIC
#mod.g1 21 -39413.44
#mod.g2 19 -39403.13
#mod.g3 19 -39414.41
#mod.g4 17 -39417.80


mod.g5<-update(mod.g4, ~. -Core)
summary(mod.g5)
mod.g6<-update(mod.g5, ~. -CpG)
summary(mod.g6)

AIC(mod.g5,mod.g6)


sum.g1<-summary(mod.g1)
sum.g2<-summary(mod.g2)

modcoef1<-sum.g1$coefficients
coef1<-modcoef1[[1]]
write.csv(coef1,"Output/SummaryStats/BetaReg_Genes_Mod.g1_TsNA0.01.csv")
modcoef2<-sum.g2$coefficients
coef2<-modcoef2[[1]]
write.csv(coef2,"Output/SummaryStats/BetaReg_Genes_Mod.g2_TsNA0.01.csv")



coef.vals <- modcoef[,1]
coef.pSE.vals <- coef.vals + modcoef[,2]
coef.mSE.vals <- coef.vals - modcoef[,2]

