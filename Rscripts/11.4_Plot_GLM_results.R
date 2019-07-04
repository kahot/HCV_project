# run GLM first (11.1 & 11.2)

library(plotrix)
library(sfsmisc)
source("Rscripts/baseRscript.R")
#####
HCVFiles_overview3<-list.files("Output_all/Overview3/",pattern="overview3.csv")


####  Make Plots  #####
source("Rscripts/GLMPlotFunctions.R")
#create alpha
cols4.60<-paste0(cols4,"66")
        
glmData<-read.csv("Output1A/GLM/GlmDataFull.Ts.Q35.csv",row.names = 1)

d<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",row.names = 1)
d<-d[d$pos>=342,]
#GLM results file: modcoeff data -use model 2 (g2)
#modcoef<-read.csv(paste0("Output1A/GLM/BetaReg_mod.g2_Ts.Q35.csv"),stringsAsFactors = F)
modcoef<-read.csv(paste0("Output1A/GLM/BetaReg_BestModel.Q35.csv"),stringsAsFactors = F)

rownames(modcoef)<-modcoef$X
modcoef<-modcoef[,-1]
modcoef<-modcoef[c("(Intercept)","t","c","g","CpG","Nonsyn","bigAAChange","t:Nonsyn","c:Nonsyn","g:Nonsyn"),]
coef.vals <- modcoef[,1]
coef.pSE.vals <- coef.vals + modcoef[,2]
coef.mSE.vals <- coef.vals - modcoef[,2]

pdf(paste0("Output1A/GLM/Mod_best.pdf"),width=10.5,height=7.5)
layout(matrix(1:2, nrow = 1))
par(mar = c(4, 4.5, 1.5, 1))
makePlot.axisbreak(main = "Synonymous Sites")
# syn=0, nonsyn=1, maesCpG=1
plotDat(0, 1, 0,cols4.60[1], .1)
plotDat(0, 0, 0,cols4.60[2], -.1)
plotVals(0,1,0, cols4[1], .1)
plotVals(0,0,0, cols4[2], -.1)
        
abline(v = 1:3 + .5, col = "black")
#legend("topleft", c("CpG-forming", "non-CpG-forming"), col = cols[1:2], pch = 16, bg = "white")
legend("bottomright", c("No drastic AA change (non-CpG-forming)", "No drastic AA change (CpG-forming)", "Drastic AA change (non-CpG-forming)",  "Drastic AA change (CpG-forming)"), 
        col = cols4[c(2,1,4,3)], pch = 16, bg = "white" )
makePlot.axisbreak(main = "Non-synonymous Sites")
plotDat(1, 1, 0, cols4.60[1], -.1)
plotDat(1, 0, 0, cols4.60[2], -.3)
plotDat(1, 1, 1, cols4.60[4], .3)
plotDat(1, 0, 1, cols4.60[3], .1)

plotVals(1,0,0, cols4[2], -.3)
plotVals(1,1,0, cols4[1], -.1)
plotVals(1,1,1, cols4[4], .3)
plotVals(1,0,1, cols4[3], .1)

abline(v = 1:3 + .5, col = "black")
dev.off()








### Try plotting Stop sites (only bigAAchange=1 and C or G only)
Type=3
CpGorNo=0
bigAAChangeOrNo=1
colVal=cols2[1]
offsetval=.1
plotDat2 <- function(Type, CpGorNo, bigAAChangeOrNo, colVal, offsetval){
        
        pchpar <- 16
        cexpar <- .6
        
        nsinds <- which(d$Type == c("syn", "nonsyn","stop")[Type])
        cpginds <- which(d$makesCpG == CpGorNo)
        aachangeinds <- which(d$bigAAChange == bigAAChangeOrNo)
        
        datInds <- intersect(intersect(nsinds, cpginds), aachangeinds)
        
        xinds <- as.character(d$ref[datInds])
        for(i in 1:4){
                xinds[xinds == nucord[i]] <- i
        }
        xinds <- as.numeric(xinds)
        yinds <- d$mean[datInds]
        xjit <- rnorm(length(xinds), offsetval, .03) #jitter
        
        yinds[yinds == 0] <- 10^(-4.4)
        
        points(xinds + xjit, yinds, col = t.col(colVal),  pch = pchpar, cex = cexpar)
}

makePlot.axisbreak(main = "Stop Sites")
plotDat2(3, 1, 1, cols2[1], .1)
plotDat2(3, 0, 1, cols2[2], -.1)
plotDat2(3, 0, 0, cols2[2], -.1)
plotDat2(3, 0, 0, cols2[2], -.1)
plotDat2(3, 0, 0, cols2[2], -.1)
plotVals(0,1,0, cols2[1], .1)
plotVals(0,0,0, cols2[2], -.1)
