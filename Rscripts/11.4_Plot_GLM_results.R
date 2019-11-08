# run GLM first (11.1 & 11.2)

library(plotrix)
library(sfsmisc)
library(colorspace)
source("Rscripts/baseRscript.R")

colors2<-qualitative_hcl(6, palette="Dark3")

#####
HCVFiles_overview3<-list.files("Output_all/Overview3/",pattern="overview3.csv")


####  Make Plots  #####
<<<<<<< HEAD
source("Rscripts/BetaRefPlotFunctions.R")

=======

source("Rscripts/GLMPlotFunctions.R")
>>>>>>> Updated analysis scrits
#create alpha
cols4.60<-paste0(cols4,"66")
        
glmData<-read.csv("Output1A/GLM/BetaRegFull.Ts.FilteredData.csv",row.names = 1)

d<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",row.names = 1)
d<-d[d$pos>=342,]
#GLM results file: modcoeff data -use model 2 (g2)
modcoef<-read.csv(paste0("Output1A/GLM/BetaReg_mod.g2_Ts.Q35.csv"),stringsAsFactors = F)
#modcoef<-read.csv(paste0("Output1A/GLM/BetaReg_BestModel.Q35.csv"),stringsAsFactors = F)

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


######
effects1<-read.csv("Output1A/GLM/BetaReg.effects.csv", stringsAsFactors = F)
effects1$factor<-factor(effects1$factor, levels=effects1$factor[1:17])
effects1$percent<-effects1$percent*100

ggplot(effects1[effects1$genotype=="1A",], aes(factor,percent)) +
        geom_bar(stat="identity", color=colors2[4], fill=paste0(colors2[4],"CC"))+
        theme_test() +
        theme(axis.text=element_text(size=13), axis.title.y=element_text(size=13))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
        theme(panel.grid.major.y = element_line(color="gray80",linetype=5))+
        labs(x="", y="Estimated effects (%)")

ggsave("Output1A/SummaryFig.Filtered/BetaRef.effects.pdf", width = 9, height = 6)

####
# Set up data 
modcoef<-read.csv(paste0("Output1A/GLM/BetaReg_BestModel.Q35.csv"),stringsAsFactors = F)
rownames(modcoef)<-modcoef$X
modcoef<-modcoef[,-1]
coef.vals <- modcoef[,1]
coef.pSE.vals <- coef.vals + modcoef[,2]
coef.mSE.vals <- coef.vals - modcoef[,2]

d<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",row.names = 1)
d<-d[d$pos>=342,c("pos","mean","ref","Type","makesCpG","bigAAChange")]

genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector

d<-merge(d, genetable, by="pos")

<<<<<<< HEAD
=======

####
# plot syn vs nonsyn estiamted and observed mf
source("Rscripts/BetaRegPlotFunction2.R")
pdf("Output1A/GLM/Betareg.plot.pdf", width = 5.5,height = 5)
col.par <- "gray95"
makePlotAll(main="")
plotMFs(0, paste0(colors2[5],"66"),  -.1)
plotMFs(1, paste0(colors2[1],"33"),  .1)

plotEsp(0,1, -.1)
plotEsp(1,1, .1)

legend("topright", c("syn", "nonsyn"), col = colors2[c(5,1)], pch = 16, bg = "white", cex=.8, bty = "white" )
dev.off()


legend("topright", c("syn", "nonsyn"), col = colors2[c(5,1)], pch = 16, cex=1.5,bg = "white", bty = "white" )



####
pdf("Output1A/GLM/TtoC.2genes.pdf", width = 5.5, height = 6)

col.par <- "gray95"
layout(matrix(1:2, nrow = 1))
par(mfrow=c(1,2), mai = c(.5, 0.7, .5, 0.1))
makePlot1(main=expression("T" %->% "C : Core"))

#plotPoints(NT,NsOrNo, CpGorNo, Core, HVR1,colVal, offset)
plotDots("T",0,0,1,0, paste0(colors2[1],"66"),  -.1)
plotDots("T",0,1,1,0, paste0(colors2[4],"66"),  .1)
plotDots("T",1,0,1,0, paste0(colors2[1],"66"),  -.1)
plotDots("T",1,1,1,0, paste0(colors2[4],"66"),  .1)

plotPoint("T",0,0,1,0, colors2[1], -.1)
plotPoint("T",0,1,1,0, colors2[4], .1)
plotPoint("T",1,0,1,0, colors2[1], -.1)
plotPoint("T",1,1,1,0, colors2[4], .1)


makePlot2(main=expression("T" %->% "C : HVR1"))
legend("bottomright", c("nonCpG", "CpG"), col = colors2[c(1,4)], pch = 16, bg = "white", cex=.7, bty = "n" )

plotDots("T",0,0,0,1, paste0(colors2[1],"66"),  -.1)
plotDots("T",0,1,0,1, paste0(colors2[4],"66"),  .1)
plotDots("T",1,0,0,1, paste0(colors2[1],"66"),  -.1)
plotDots("T",1,1,0,1, paste0(colors2[4],"66"),  .1)

plotPoint("T",0,0,0,1, colors2[1], -.1)
plotPoint("T",0,1,0,1, colors2[4], .1)
plotPoint("T",1,0,0,1, colors2[1], -.1)
plotPoint("T",1,1,0,1, colors2[4], .1)
dev.off()


## same color as HIV
pdf(paste0("Output1A/GLM/Mod_samecolorHIV.pdf"),width=10.5,height=6.5)
layout(matrix(1:2, nrow = 1))
par(mar = c(4, 4.5, 1.5, 1))
makePlot.axisbreak(main = "Synonymous Sites")
# syn=0, nonsyn=1, maesCpG=1
plotDat(0, 1, 0,cols[1], .1)
plotDat(0, 0, 0,cols[2], -.1)
plotVals(0,1,0, cols[1], .1)
plotVals(0,0,0, cols[2], -.1)

abline(v = 1:3 + .5, col = "black")
#legend("topleft", c("CpG-forming", "non-CpG-forming"), col = cols[1:2], pch = 16, bg = "white")
legend("bottomright", c("No drastic AA change (non-CpG-forming)", "No drastic AA change (CpG-forming)", "Drastic AA change (non-CpG-forming)",  "Drastic AA change (CpG-forming)"), 
       col = cols[c(2,1,3,4)], pch = 16, bg = "white" )
makePlot.axisbreak(main = "Non-synonymous Sites")
plotDat(1, 1, 0, cols[1], -.1)
plotDat(1, 0, 0, cols[2], -.3)
plotDat(1, 1, 1, cols[4], .3)
plotDat(1, 0, 1, cols[3], .1)

plotVals(1,0,0, cols[1], -.3)
plotVals(1,1,0, cols[2], -.1)
plotVals(1,1,1, cols[4], .3)
plotVals(1,0,1, cols[3], .1)

abline(v = 1:3 + .5, col = "black")
dev.off()


######
effects1<-read.csv("Output1A/GLM/BetaReg.effects.csv", stringsAsFactors = F)
effects1$factor<-factor(effects1$factor, levels=rev(effects1$factor[1:17]))
effects1$percent<-effects1$percent*100

ggplot(effects1[effects1$genotype=="1A",], aes(factor,percent)) +
        geom_bar(stat="identity",width = .8, color=colors2[4], fill=paste0(colors2[4],"CC"))+
        theme_test() +
        geom_hline(yintercept = 0, color = "gray70", size=.3)+
        theme(axis.text=element_text(size=13, color="black"), axis.title.x = element_text(size=14))+
        theme(panel.grid.major.y = element_line(color="gray80",linetype=5))+
        labs(x="", y="Estimated effects (%)")+
        coord_flip()


ggsave("Output1A/SummaryFig.Filtered/BetaRef.effects.vertical.pdf", width = 5, height = 9)



### for slides
ggplot(effects1[effects1$genotype=="1A",], aes(factor,percent)) +
        geom_bar(stat="identity",width = .8, color=colors2[4], fill=paste0(colors2[4],"CC"))+
        theme_test() +
        geom_hline(yintercept = 0, color = "gray70", size=.3)+
        theme(axis.text.y=element_text(size=13, color="black"),axis.text.x=element_blank(), axis.title.x = element_blank())+
        theme(panel.grid.major.y = element_line(color="gray80",linetype=5))+
        labs(x="", y="Estimated effects (%)")+
        coord_flip()
ggsave("Output1A/SummaryFig.Filtered/BetaRef.effects.vertical2.pdf", width = 5, height = 9)





####
# Set up data 
modcoef<-read.csv(paste0("Output1A/GLM/BetaReg_BestModel.Q35.csv"),stringsAsFactors = F)
rownames(modcoef)<-modcoef$X
modcoef<-modcoef[,-1]
coef.vals <- modcoef[,1]
coef.pSE.vals <- coef.vals + modcoef[,2]
coef.mSE.vals <- coef.vals - modcoef[,2]

d<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",row.names = 1)
d<-d[d$pos>=342,c("pos","mean","ref","Type","makesCpG","bigAAChange")]

genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector

d<-merge(d, genetable, by="pos")

>>>>>>> Updated analysis scrits
#pdf(paste0("Output1A/GLM/Mod_best.pdf"),width=10.5,height=7.5)
#par(mar = c(4, 4.5, 1.5, 1))
#makePlot.axisbreak(main =expression("T" %->% "C : Core"))
pdf("Output1A/GLM/TtoC.2genes.pdf", width = 5.5, height = 6)

col.par <- "gray95"
layout(matrix(1:2, nrow = 1))
par(mfrow=c(1,2), mai = c(.5, 0.7, .5, 0.1))
makePlot1(main=expression("T" %->% "C : Core"))

#plotPoints(NT,NsOrNo, CpGorNo, Core, HVR1,colVal, offset)
plotDots("T",0,0,1,0, paste0(colors2[1],"66"),  -.1)
plotDots("T",0,1,1,0, paste0(colors2[4],"66"),  .1)
plotDots("T",1,0,1,0, paste0(colors2[1],"66"),  -.1)
plotDots("T",1,1,1,0, paste0(colors2[4],"66"),  .1)

plotPoint("T",0,0,1,0, colors2[1], -.1)
plotPoint("T",0,1,1,0, colors2[4], .1)
plotPoint("T",1,0,1,0, colors2[1], -.1)
plotPoint("T",1,1,1,0, colors2[4], .1)


makePlot2(main=expression("T" %->% "C : HVR1"))
legend("bottomright", c("nonCpG", "CpG"), col = colors2[c(1,4)], pch = 16, bg = "white", cex=.7, bty = "n" )

plotDots("T",0,0,0,1, paste0(colors2[1],"66"),  -.1)
plotDots("T",0,1,0,1, paste0(colors2[4],"66"),  .1)
plotDots("T",1,0,0,1, paste0(colors2[1],"66"),  -.1)
plotDots("T",1,1,0,1, paste0(colors2[4],"66"),  .1)

plotPoint("T",0,0,0,1, colors2[1], -.1)
plotPoint("T",0,1,0,1, colors2[4], .1)
plotPoint("T",1,0,0,1, colors2[1], -.1)
plotPoint("T",1,1,0,1, colors2[4], .1)
dev.off()

<<<<<<< HEAD
=======





>>>>>>> Updated analysis scrits
###
pdf("Output1A/GLM/AtoG.2genes.pdf", width = 5.5, height = 6)
layout(matrix(1:2, nrow = 1))
par(mfrow=c(1,2), mai = c(.5, 0.7, .5, 0.1))
makePlot1(main=expression("A" %->% "G : Core"))

#plotPoints(NT,NsOrNo, CpGorNo, Core, HVR1,colVal, offset)
plotDots("A",0,0,1,0, paste0(colors2[1],"66"),  -.1)
plotDots("A",0,1,1,0, paste0(colors2[4],"66"),  .1)
plotDots("A",1,0,1,0, paste0(colors2[1],"66"),  -.1)
plotDots("A",1,1,1,0, paste0(colors2[4],"66"),  .1)

plotPoint("A",0,0,1,0, colors2[1], -.1)
plotPoint("A",0,1,1,0, colors2[4], .1)
plotPoint("A",1,0,1,0, colors2[1], -.1)
plotPoint("A",1,1,1,0, colors2[4], .1)


makePlot2(main=expression("A" %->% "G : HVR1"))
legend("bottomright", c("nonCpG", "CpG"), col = colors2[c(1,4)], pch = 16, bg = "white", cex=.7, bty = "n" )

plotDots("A",0,0,0,1, paste0(colors2[1],"66"),  -.1)
plotDots("A",0,1,0,1, paste0(colors2[4],"66"),  .1)
plotDots("A",1,0,0,1, paste0(colors2[1],"66"),  -.1)
plotDots("A",1,1,0,1, paste0(colors2[4],"66"),  .1)

plotPoint("A",0,0,0,1, colors2[1], -.1)
plotPoint("A",0,1,0,1, colors2[4], .1)
plotPoint("A",1,0,0,1, colors2[1], -.1)
plotPoint("A",1,1,0,1, colors2[4], .1)

dev.off()


#Test
library(MASS)
library(DescTools)

eff1A<-effects2[effects2$genotype=="1A",c("factor","effect")]

eff3A<-effects2[effects2$genotype=="3A",c("factor","effect")]
#chisquare test
eff.table<-eff1A
row.names(eff.table)<-eff.table$factor 
eff.table$Geno1B<-eff3A$effect
colnames(eff.table)[2]<-"Geno1A"
eff.table<-eff.table[,-1]

chisq.test(eff.table)
#X-squared = 0.00069968, df = 16, p-value = 1







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
