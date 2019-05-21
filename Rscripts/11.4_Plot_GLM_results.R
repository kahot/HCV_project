# run GLM first (11.1 & 11.2)

library(plotrix)
library(sfsmisc)

#####
HCVFiles_overview3<-list.files("Output_all/Overview3/",pattern="overview3.csv")


####  Make Plots  #####
source("Rscripts/GLMPlotFunctions.R")


cols2<-c("#66CCEE66", "#22883366" ,"#CCBB4466", "#EE667766", "#AA337766")
filename<-c("Ts.Q35","Ts_NA.Q35","Ts_zero.Q35")
glmData<-read.csv("Output/GLM/GlmData_Ts.csv")
glmData<-glmData[,-1]

for (i in 1:3){
        d<-read.csv(paste0("Output/Mut.freq.filtered/Summary_",filename[i],".csv"))
        d<-d[which(d$pos==342):nrow(d),-1]

       
        #specify the modcoeff dataframe:
        modcoef<-read.csv(paste0("Output/GLM/BetaReg_mod.g2_",filename[i],".csv"),stringsAsFactors = F)
        rownames(modcoef)<-modcoef$X
        modcoef<-modcoef[,-1]
        modcoef<-modcoef[c("(Intercept)","t","c","g","CpG","Nonsyn","bigAAChange","t:Nonsyn","c:Nonsyn","g:Nonsyn"),]

        coef.vals <- modcoef[,1]
        coef.pSE.vals <- coef.vals + modcoef[,2]
        coef.mSE.vals <- coef.vals - modcoef[,2]

        pdf(paste0("Output/GLM/Mod2_",filename[i],".pdf"),width=12,height=7.5)
        layout(matrix(1:2, nrow = 1))
        par(mar = c(4, 4.5, 1.5, 1))
        makePlot.axisbreak(main = "Synonymous Sites")
        # syn=0, nonsyn=1, maesCpG=1
        plotDat(0, 1, 0, cols2[1], .1)
        plotDat(0, 0, 0, cols2[2], -.1)
        plotVals(0,1,0, cols[1], .1)
        plotVals(0,0,0, cols[2], -.1)
        
        abline(v = 1:3 + .5, col = "black")
        #legend("topleft", c("CpG-forming", "non-CpG-forming"), col = cols[1:2], pch = 16, bg = "white")
        legend("bottomright", c("No drastic AA change (non-CpG-forming)", "No drastic AA change (CpG-forming)", "Drastic AA change (non-CpG-forming)",  "Drastic AA change (CpG-forming)"), col = cols[c(2,1,3,4)], pch = 16, bg = "white" )
        makePlot.axisbreak(main = "Non-synonymous Sites")
        plotDat(1, 0, 1, cols2[3], .1)
        plotDat(1, 0, 0, cols2[2], -.3)
        plotDat(1, 1, 0, cols2[1], -.1)
        plotDat(1, 1, 1, cols2[4], .3)
        plotVals(1,0,1, cols[3], .1)
        plotVals(1,0,0, cols[2], -.3)
        plotVals(1,1,0, cols[1], -.1)
        plotVals(1,1,0, cols[4], .3)

        abline(v = 1:3 + .5, col = "black")
        dev.off()
}





