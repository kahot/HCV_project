library(zoo)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(ggthemes)
library(sfsmisc)
library(colorspace)
library(cowplot)

source("Rscripts/baseRscript.R")
source("Rscripts/label_scientific.R")


colors2<-qualitative_hcl(6, palette="Dark3")
scaleFUN <- function(x) sprintf("%.2f", x)
col2_light<-qualitative_hcl(6, palette="Set3")
###########

df<-read.csv("Output1A/SelCoeff/SC.csv", stringsAsFactors = F, row.names = 1)
#coding regions only
df<-df[df$pos>=342,]

#add gene ingo
genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
genes$Gene[6]<-"NS1"
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector
end<-df$pos[nrow(df)]
genetable<-genetable[genetable$pos>=342&genetable$pos<=end,]
sc<-merge(df, genetable, by="pos")

#CpG creating mutation freq.
cpg<-list()
cpg.se<-list()
for (i in 1:11){
        cpg[i]<-mean(sc$mean[sc$gene==genes$Gene[i+1]& sc$makesCpG==1], na.rm=T)
        names(cpg)[i]<-genes$Gene[i+1]
        cpg.se[i]<-std.error(sc$mean[sc$gene==genes$Gene[i+1]& sc$makesCpG==1], na.rm=T)
        names(cpg.se)[i]<-genes$Gene[i+1]
}
cpg.mf<-as.data.frame(do.call(rbind,cpg))
cpgse<-as.data.frame(do.call(rbind,cpg.se))
cpg.mf$cpg.se<-cpgse[,1]
cpg.mf$gene<-rownames(cpg.mf)
colnames(cpg.mf)<-c("CpG", "CpG.se", "gene")

#CpG<-cpg.mf[,c(3,1,2)]
#colnames(CpG)[2:3]<-c("MF", "SE")
#CpG$Kind<-"CpG"

ag<-sc[sc$ref=="a" | sc$ref=="t",]
AveMF<-aggregate(ag$mean,by=list(ag$gene), mean, na.rm = TRUE)
colnames(AveMF)<-c("gene","MF")
SEMF<-aggregate(ag$mean,by=list(ag$gene), std.error, na.rm = TRUE)
colnames(SEMF)<-c("gene","mf.se")
mf<-merge(AveMF, SEMF, by="gene")
Cpg.sum<-merge(mf,cpg.mf, by="gene")


mfA<-sc[sc$ref=="a",]
mfT<-sc[sc$ref=="t",]
aveAs<-aggregate(mfA$mean,by=list(mfA$gene), mean, na.rm = TRUE)
colnames(aveAs)<-c("gene","MF.A")
aveTs<-aggregate(mfT$mean,by=list(mfT$gene), mean, na.rm = TRUE)
colnames(aveTs)<-c("gene","MF.T")
ATs<-merge(aveAs, aveTs, by="gene")

Cpg.sum<-merge(Cpg.sum, ATs, by="gene")
Cpg.sum$Ratio<-Cpg.sum$CpG/rowMeans()
Cpg.sum$Ratio2<-Cpg.sum$CpG/Cpg.sum$MF
Cpg.sum$mean.AT<-rowMeans(Cpg.sum[,7:8]) 


Cpg.sum$gene<-factor(Cpg.sum$gene, levels=c("Core", "E1", "HVR1", "E2","NS1","NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

cpgPlot1<-ggplot(Cpg.sum,aes(x=gene,y=Ratio))+
        scale_y_continuous(labels=scaleFUN)+
        geom_point(color=colors2[3], size=3)+
        theme_bw()+
        theme(axis.title.x=element_blank())+ylab("CpG MF ratio")+
        geom_hline(yintercept = 1, col="gray50")+
        theme(panel.grid.major.x = element_blank())

cpgPlot2<-ggplot(Cpg.sum,aes(x=gene,y=CpG))+
        scale_y_continuous(labels=scaleFUN)+
        geom_point(color=colors2[3], size=3)+
        theme_bw()+
        theme(axis.title.x=element_blank())+ylab("CpG creating mut freq")+
        theme(panel.grid.major.x = element_blank())

###
AveSC<-aggregate(sc$EstSC,by=list(sc$gene), mean, na.rm = TRUE)
colnames(AveSC)<-c("gene","SC")
SESC<-aggregate(sc$EstSC,by=list(sc$gene), std.error, na.rm = TRUE)
colnames(SESC)<-c("gene","se")

#
#AT content of each genes:
ATcontents<-list()
for (i in 1:11){
        genename<-genes$Gene[i+1]
        seq1<-as.character(sc$ref[sc$gene==genename])
        ATcontents[i]<-1-GC(seq1)
        names(ATcontents)[i]<-genename
        
}

ATs<-as.data.frame(do.call(rbind,ATcontents))
ATs$gene<-rownames(ATs)
colnames(ATs)[1]<-"AT"


SC_summary2<-merge(AveSC,SESC,by="gene")
SC_summary2<-merge(SC_summary2,ATs,by="gene")
SC_summary2<-merge(SC_summary2,Cpg.sum, by="gene")
write.csv(SC_summary2, "Output1A/SummaryStats/SC-AT-CpG-oddsRatio.csv")
SC_summary2$gene<-factor(SC_summary2$gene, levels=c("Core", "E1", "HVR1", "E2","NS1","NS2","NS3","NS4A","NS4B","NS5A","NS5B"))




scPlot<-ggplot(SC_summary2,aes(gene,y=SC))+
        geom_point(color="royalblue", size=3)+
        geom_errorbar(aes(ymin=SC-se, ymax=SC+se), width=.3,color="royalblue" )+
        theme_bw()+
        theme(axis.title.x=element_blank())+ylab("Average selection coefficient")+
        theme(panel.grid.major.x = element_blank())

ATplot<-ggplot(SC_summary2,aes(x=gene,y=AT))+
        geom_point(color="#FF7A16", size=3.5)+
        theme_bw()+
        theme(axis.title.x=element_blank())+ylab("AT content")+
        theme(panel.grid.major.x = element_blank())



title <- ggdraw() + draw_label("Selection coefficients,AT contents, CpG", fontface='bold')
plot_grid(title,scPlot, ATplot, cpgPlot1, nrow = 4, labels = c("", "","", ""),
          rel_heights = c(0.2, 1, 1,1))
ggsave(filename="Output1A/SelCoeff/SC-AT-CpG.ratio.pdf",width =7, height =6)



title <- ggdraw() + draw_label("Selection coefficients,AT contents, CpG creating MF", fontface='bold')
plot_grid(title,scPlot, ATplot, cpgPlot2, nrow = 4, labels = c("", "","", ""),
          rel_heights = c(0.2, 1, 1,1))
ggsave(filename="Output1A/SelCoeff/SC-AT-CpG.pdf",width =7, height =6)


## Correlation:

results1<-cor.test(SC_summary2$AT,SC_summary2$Ratio, method = "spearman")
print(results1)
# rho =-0.3727273
# p-value = 0.2606

results2<-cor.test(SC_summary2$AT,SC_summary2$SC, method = "spearman")
print(results2)
#      rho 
#0.2181818 
#p-value = 0.5209

results3<-cor.test(SC_summary2$AT,SC_summary2$CpG, method = "spearman")
print(results3)
results3<-cor.test(SC_summary2$GC,SC_summary2$CpG, method = "spearman")
print(results3)
#S = 372, p-value = 0.02306
#rho -0.6909091 


## Plot the relationship
#GC content of each genes:
GCcontents<-list()
for (i in 1:11){
        genename<-genes$Gene[i+1]
        seq1<-as.character(sc$ref[sc$gene==genename])
        GCcontents[i]<-GC(seq1)
        names(GCcontents)[i]<-genename
        
}

GCs<-as.data.frame(do.call(rbind,GCcontents))
GCs$gene<-rownames(GCs)
colnames(GCs)[1]<-"GC"

SC_summary2<-merge(SC_summary2, GCs, by="gene")
results3<-cor.test(SC_summary2$GC,SC_summary2$CpG, method = "spearman")
print(results3)

sum<-merge(SC_summary2, Cpg.sum[,c(1,6:10)], by="gene")
cor.test(sum$Ratio.y,sum$SC, method = "spearman")


pdf("Output1A/SummaryFig.Filtered/GC-Cpgcreting.mf.pdf", height = 5, width = 5)
plot(SC_summary2$GC,SC_summary2$CpG, ylab="", xlab="",pch=16,col=colors2[5],cex=3)
mtext("GC Content",1, 2.5, cex=1.2)
mtext("CpG-creating mutation frequency",2,2.2, cex=1.2)
abline(lm(SC_summary2$CpG~SC_summary2$GC), col = "gray70")
text(x=max(SC_summary2$GC),y=min(SC_summary2$CpG)+0.0002, labels=expression(paste(rho, " = 0.691***")),cex=1.1, adj=1)
dev.off()


# Mut freq and GC contents are not related, (but CpG creating mutations are)
cor.test(SC_summary2$GC,SC_summary2$MF, method = "spearman")
#p-value = 0.7138
#rho  0.1272727 


AveMF<-aggregate(sc$mean,by=list(sc$gene), mean, na.rm = TRUE)
colnames(AveMF)<-c("gene","Ts.MF")
SEMF<-aggregate(sc$mean,by=list(sc$gene), std.error, na.rm = TRUE)
colnames(SEMF)<-c("gene","mf.se")
mf<-merge(AveMF, SEMF, by="gene")

SC_summary2<-merge(SC_summary2, mf, by="gene")
cor.test(SC_summary2$GC,SC_summary2$Ts.MF, method = "spearman")
#p-value = 0.9244
#rho =-0.03636364 

#### 
