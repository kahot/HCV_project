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
mean(df$EstSC) # 0.002681807  

#coding regions only
df<-df[df$pos>=342,]
mean(df$EstSC) #0.002602148
mean(df$EstSC[df$Type=="syn"])  #0.001128887
mean(df$EstSC[df$Type=="nonsyn"])  #0.003327282
std.error(df$EstSC) #1.961157e-05
std.error(df$EstSC[df$Type=="syn"])  #1.40941e-05
std.error(df$EstSC[df$Type=="nonsyn"])  #2.330815e-05


#### Calcualte the mean SC and GC_contents for each gene
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
AveSC<-aggregate(sc$EstSC,by=list(sc$gene), mean, na.rm = TRUE)
colnames(AveSC)<-c("gene","SC")
SESC<-aggregate(sc$EstSC,by=list(sc$gene), std.error, na.rm = TRUE)
colnames(SESC)<-c("gene","se")

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

#CpG creating mutation freq.
#mut.freq & se
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

# Total A/T mut freq
ag<-sc[sc$ref=="a" | sc$ref=="t",]
AveMF<-aggregate(ag$mean,by=list(ag$gene), mean, na.rm = TRUE)
colnames(AveMF)<-c("gene","MF")
SEMF<-aggregate(ag$mean,by=list(ag$gene), std.error, na.rm = TRUE)
colnames(SEMF)<-c("gene","mf.se")
mf<-merge(AveMF, SEMF, by="gene")


# create data.frame for plotting cpg mf and mf together
CpG<-cpg.mf[,c(3,1,2)]
colnames(CpG)[2:3]<-c("MF", "SE")
CpG$Kind<-"CpG"
mutf<-mf
colnames(mutf)[2:3]<-c("MF", "SE")
mutf$Kind<-"All"
CpG<-rbind(CpG,mutf)        
CpG$gene<-factor(CpG$gene, levels=c("Core", "E1", "HVR1", "E2","NS1","NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

#proportion of CpG creating mutations
#cpgPercent<-list()
#for (i in 1:12){
#        dt<-sc[sc$gene==genes$Gene[i],]
#        cpgPercent[i]<-table(dt$makesCpG)[2]/nrow(dt)*100
#        names(cpgPercent)[i]<-genes$Gene[i]
#}
#CpGs<-as.data.frame(do.call(rbind,cpgPercent))
#CpGs$gene<-rownames(CpGs)
#colnames(CpGs)[1]<-"CpG"


SC_summary<-merge(AveSC,SESC,by="gene")
SC_summary<-merge(SC_summary,GCs,by="gene")
SC_summary<-merge(SC_summary,ATs,by="gene")
SC_summary<-merge(SC_summary,cpg.mf, by="gene")
SC_summary<-merge(SC_summary,mf, by="gene")

write.csv(SC_summary,"Output1A/SelCoeff/SC.GCcontents_summary_by_genes.csv")

###########  Plot the average Sel Coeff by genes
SC_summary$gene<-factor(SC_summary$gene, levels=c("Core", "E1", "HVR1", "E2","NS1","NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

#ggboxplot(sc,x="gene",y="EstSC", xlab="Gene", ylab="Estimated selective coefficient",color="gene")+
#        theme(legend.position="none")
#ggboxplot(sc,x="gene",y="EstSC", xlab="Gene", ylab="Estimated selective coefficient",color="gray40")+
#        theme(legend.position="none")
#ggsave(filename=paste0("Output1A/SelCoeff/SC_by.gene.Boxplot.pdf"), width = 10, height = 7)

ggplot(SC_summary,aes(gene,y=SC))+
        geom_point(color="royalblue", size=3)+
        geom_errorbar(aes(ymin=SC-se, ymax=SC+se), width=.2,color="royalblue" )+
        theme_bw()+
        theme(axis.title.x=element_blank())+ylab("Average selection coefficient")+
        theme(panel.grid.major.x = element_blank())

ggsave(filename="Output1A/SelCoeff/Ave.SC_by.gene.pdf",width =5, height = 4)



##
GCplot<-ggplot(SC_summary,aes(x=gene,y=GC))+
        geom_point(color="#FF7A16", size=3)+
        theme_bw()+
        theme(axis.title.x=element_blank())+ylab("GC content")+
        theme(panel.grid.major.x = element_blank())

scPlot<-ggplot(SC_summary,aes(gene,y=SC))+
        geom_point(color="royalblue", size=3)+
        geom_errorbar(aes(ymin=SC-se, ymax=SC+se), width=.3,color="royalblue" )+
        theme_bw()+
        theme(axis.title.x=element_blank())+ylab("Average selection coefficient")+
        theme(panel.grid.major.x = element_blank())

ATplot<-ggplot(SC_summary,aes(x=gene,y=AT))+
        geom_point(color="#FF7A16", size=3.5)+
        theme_bw()+
        theme(axis.title.x=element_blank())+ylab("AT content")+
        theme(panel.grid.major.x = element_blank())


cpgPlot<-ggplot(CpG,aes(x=gene,y=MF, color=Kind))+
        scale_y_continuous(labels=scaleFUN)+scale_color_manual(values=colors2[c(3,6)])+
        geom_point(position=position_dodge(width=0.8), size=3)+
        geom_errorbar(aes(ymin=MF-SE, ymax=MF+SE), width=.2, size=.2,position=position_dodge(width=0.8))+
        theme_bw()+
        theme(axis.title.x=element_blank())+ylab("Mutation frequency")+
        theme(panel.grid.major.x = element_blank())+
        theme(legend.position = "none") 


title <- ggdraw() + draw_label("Selection coefficients,AT contents, CpG creating MF", fontface='bold')
plot_grid(title,scPlot, ATplot, cpgPlot, nrow = 4, labels = c("", "","", ""),
          rel_heights = c(0.2, 1, 1,1))
ggsave(filename="Output1A/SelCoeff/SC-AT-CpG.pdf",width =7, height =6)





### PLOT SCs across the genome
startnuc<-sc$pos[1]
endnuc<- sc$pos[nrow(sc)]  # from 11.Mut.Freq calculation. If using the last row -> Trans$pos[nrow(Trnas)]
n<-data.frame("pos"=c(startnuc:endnuc))

SCs<-merge(n,sc,by="pos",all.x=T) 
range(SCs$EstSC,na.rm = T) #0.0001359244 0.0793168573
ylow=0.00008
yhigh=0.09

title<-"SC_acroosGenoem1"
pdf(paste0("Output1A/SelCoeff/",title,".pdf"),width=15,height=4.5)
plot(EstSC~pos, data=SCs,t="n",log='y',yaxt='n',xlab='Genome position',ylab="Estimated selection coefficient",
     ylim=c(ylow,yhigh),xlim=c(265,8618))
eaxis(side = 2, at = 10^((0):(-(5))), cex=2)
for(k in 1:5){abline(h = 1:10 * 10^(-k), col = "gray80")}

#add rolling average
#roll100<-rollmean(SCs$mean, k=100, na.rm=T)
#SCs$roll100<-c(rep(NA, times=50),roll100,c(rep(NA, times=49)))
#lines(roll100~pos,data=SCs, col="#4477AA",lwd=1.5)

# roling average of 50
roll50<-rollmean(SCs$EstSC, k=50, na.rm=T)
SCs$roll50<-c(rep(NA, times=25),roll50,c(rep(NA, times=24)))

#scColors<-c("#AEE5F6","#4477AA","#8CD0C8","#009988","#228833")

scColors<-c("#AEE5F6","#4477AA","#EEBAB9","#EE6677","#228833")

for (i in 1:(nrow(genes)-1)){
        gname<-genes$Gene[i]
        if (i%%2==0) {
                points(SCs$EstSC[SCs$gene==gname]~SCs$pos[SCs$gene==gname], pch=20,col=scColors[1],cex=0.5)
                lines(SCs$roll50[SCs$gene==gname]~SCs$pos[SCs$gene==gname], col=scColors[2],lwd=1.5)}
        if (i%%2==1) {
                points(SCs$EstSC[SCs$gene==gname]~SCs$pos[SCs$gene==gname], pch=20,col=scColors[3],cex=0.5)
                lines(SCs$roll50[SCs$gene==gname]~SCs$pos[SCs$gene==gname], col=scColors[4],lwd=1.5)}
        #if (i%%3==2) lines(SCs$roll50[SCs$gene==gname]~SCs$pos[SCs$gene==gname], col=scColors[3],lwd=1.5)
}

for (j in 1:(nrow(genes)-1)){
        xleft<-genes$start[j]
        xright<-genes$start[j+1]
        if (j%%2==1){
                rect(xleft,ylow,xright,ylow*1.4,density = NULL, angle = 45,col=scColors[3],border ="gray20",lwd=1.5)
                if (j==1) text(200,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                if (j==9) text(xleft+80,1.8*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                else {text(xright-(xright-xleft)/2,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)}
        }
        if (j%%2==0){
                rect(xleft,ylow,xright,1.4*ylow,density = NULL, angle = 45,col=scColors[1],border ="gray20",lwd=1.5)
                if (j==4|j==6) text(xleft+50,1.8*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                if (j==12) text(xleft+700,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                else text(xright-(xright-xleft)/2,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        }
        
        
}


       
box()
dev.off()

#no5'UTR

title<-"SC_acroosGenoem2"
yhigh=0.03
pdf(paste0("Output1A/SelCoeff/",title,".pdf"),width=15,height=4.5)
plot(EstSC~pos, data=SCs,t="n",log='y',yaxt='n',xlab='Genome position',ylab="Estimated selection coefficient",
     ylim=c(ylow,yhigh),xlim=c(342,8618))
eaxis(side = 2, at = 10^((0):(-(5))), cex=2)
for(k in 1:5){abline(h = 1:10 * 10^(-k), col = "gray80")}

# roling average of 50
roll50<-rollmean(SCs$EstSC, k=50, na.rm=T)
SCs$roll50<-c(rep(NA, times=25),roll50,c(rep(NA, times=24)))

#scColors<-c("#AEE5F6","#4477AA","#8CD0C8","#009988","#228833")

scColors<-c("#AEE5F6","#4477AA","#EEBAB9","#EE6677","#228833")

for (i in 2:(nrow(genes)-1)){
        gname<-genes$Gene[i]
        if (i%%2==0) {
                points(SCs$EstSC[SCs$gene==gname]~SCs$pos[SCs$gene==gname], pch=20,col=scColors[1],cex=0.5)
                lines(SCs$roll50[SCs$gene==gname]~SCs$pos[SCs$gene==gname], col=scColors[2],lwd=1.5)}
        if (i%%2==1) {
                points(SCs$EstSC[SCs$gene==gname]~SCs$pos[SCs$gene==gname], pch=20,col=scColors[3],cex=0.5)
                lines(SCs$roll50[SCs$gene==gname]~SCs$pos[SCs$gene==gname], col=scColors[4],lwd=1.5)}
        #if (i%%3==2) lines(SCs$roll50[SCs$gene==gname]~SCs$pos[SCs$gene==gname], col=scColors[3],lwd=1.5)
}

for (j in 1:(nrow(genes)-1)){
        xleft<-genes$start[j]
        xright<-genes$start[j+1]
        if (j%%2==1){
                rect(xleft,ylow,xright,ylow*1.4,density = NULL, angle = 45,col=scColors[3],border ="gray20",lwd=1.5)
                if (j==1) text(200,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                if (j==9) text(xleft+80,1.8*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                if (j!=1&j!=9) text(xright-(xright-xleft)/2,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        }
        if (j%%2==0){
                rect(xleft,ylow,xright,1.4*ylow,density = NULL, angle = 45,col=scColors[1],border ="gray20",lwd=1.5)
                if (j==4|j==6) text(xright-(xright-xleft)/2,1.8*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                if (j==12) text(xleft+700,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                if (j!=4&j!=6&j!=12) text(xright-(xright-xleft)/2,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        }
        
        
}

box()
dev.off()






##Change color for all genes
Colors<-qualitative_hcl(11, palette = "Pastel 1")
add.alpha <- function(col, alpha=1){
        if(missing(col))
                stop("Please provide a vector of colours.")
        apply(sapply(col, col2rgb)/255, 2, 
              function(x) 
                      rgb(x[1], x[2], x[3], alpha=alpha))  
}


scColorsAlpha<-add.alpha(Colors, alpha=0.5)
scColors<-qualitative_hcl(11, palette = "Dynamic")
#scColors<-c("#77AADD","#99DDFF","#44BB99","#BBCC33","#AAAA00","#EEDD88","#EE8866","#FFAABB","DDDDDD")
ylow=0.00009
yhigh=0.03
pdf(paste0("Output1A/SelCoeff/SC_acrossGenome.multicolors.pdf"),width=15,height=4.5)
plot(EstSC~pos, data=SCs,t="n",log='y',yaxt='n',xlab='',ylab="",
     ylim=c(ylow,yhigh),xlim=c(340,8618))

mtext("Genome position", side=1, line=2.5, cex=1.4)
mtext("Estimated selection coefficient", side=2, line=2.5, cex=1.4)
eaxis(side = 2, at = 10^((0):(-(5))), cex=2)
for(k in 1:5){abline(h = 1:10 * 10^(-k), col = "gray80")}

for (i in 2:(nrow(genes)-1)){
        gname<-genes$Gene[i]
        points(SCs$EstSC[SCs$gene==gname]~SCs$pos[SCs$gene==gname], pch=20,col=scColorsAlpha[i-1],cex=0.5)
        lines(SCs$roll50[SCs$gene==gname]~SCs$pos[SCs$gene==gname], col=scColors[i-1],lwd=1.5)
        }


for (j in 1:(nrow(genes)-1)){
        xleft<-genes$start[j]
        xright<-genes$start[j+1]
        rect(xleft,ylow,xright,ylow*1.4,density = NULL, angle = 45,col=scColors[j-1],border ="gray20",lwd=1.5)
        if (j==1) text(200,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        if (j==4|j==6|j==9) text(xleft+50,.88*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        else text(xleft+200,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
      }
        

box()
dev.off()







###################
##Histogram of Selection Coefficients with Median 
#source("Rscripts/SelectionCoefficientSummary.R")
#1.Transition Mutation
colors<-rep(paste0(colors2[c(1,3,5,2)], "99"), times=2)


sc2<-sc[sc$pos>=342,]
sc2<-sc2[,c("pos", "ref", "Type","EstSC")]

plots<-list()
meanSC<-list()
k=1
for (type in c("syn","nonsyn")){
for (i in c("a","t","c","g")) {
        
                dt<-sc2[sc2$Type==type & sc2$ref==i,]
                #nt<-toupper(i)
                if (i=="a") {linecol=colors2[1];fillcol=paste0(colors2[1],"99")}
                if (i=="t") {linecol=colors2[3];fillcol=paste0(colors2[3],"99")}
                if (i=="c") {linecol=colors2[5];fillcol=paste0(colors2[5],"99")}
                if (i=="g") {linecol=colors2[2];fillcol=paste0(colors2[2],"99")}
                
               
               p1<-eval(substitute(
                       ggplot(dt, aes(x=EstSC))+geom_histogram(color=linecol, fill=fillcol)+
                        scale_x_continuous(trans = 'log10',label=label_scientific, limits = c(0.00001,1))+
                        theme_bw()+
                        ylab("Count")+
                        xlab("")+
                        geom_vline(aes(xintercept=mean(EstSC)),
                                   color="blue", linetype="dashed", size=.5)
                        , list(k=k)))
                plots[[k]]<-p1
                meanSC[[k]]<-mean(dt$EstSC, na.rm=T)
                k=k+1
        }
}             
par(mfrow=c(2,4))
multiplot(plotlist=plots, layout=matrix(c(1,2,3,4,5,6,7,8), nrow=2, byrow=T))
ggsave("Output1A/SelCoeff/hist.pdf", width=8, height=6)


#for (type in c("syn","nonsyn")){
#        for (i in c("a","t","c","g")) {

type="syn"
dt<-sc2[sc2$Type=="syn" & sc2$ref=="a",]
A.syn<-ggplot(dt, aes(x=EstSC))+geom_histogram(color="gray60", fill=paste0(colors2[1],"CC"))+
                scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00001,1))+
                theme_bw()+ylab("Count")+xlab("")+ggtitle("A:Syn")+theme(plot.title = element_text(hjust = 0.5, size=12))+
                geom_vline(aes(xintercept=mean(EstSC, na.rm=T)),color="blue", linetype="dashed", size=.5)
dt<-sc2[sc2$Type=="nonsyn" & sc2$ref=="a",]
A.nonsyn<-ggplot(dt, aes(x=EstSC))+geom_histogram(color="gray60", fill=paste0(colors2[1],"CC"))+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00001,1))+
        theme_bw()+ylab("Count")+xlab("")+ggtitle("A:Nonsyn")+theme(plot.title = element_text(hjust = 0.5, size=12))+
        geom_vline(aes(xintercept=mean(EstSC, na.rm=T)),color="blue", linetype="dashed", size=.5)


dt<-sc2[sc2$Type=="syn" & sc2$ref=="t",]
T.syn<-ggplot(dt, aes(x=EstSC))+geom_histogram(color="gray60", fill=paste0(colors2[3],"CC"))+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00001,1))+
        theme_bw()+ylab("Count")+xlab("")+ggtitle("T:Syn")+theme(plot.title = element_text(hjust = 0.5, size=12))+
        geom_vline(aes(xintercept=mean(EstSC, na.rm=T)),color="blue", linetype="dashed", size=.5)
dt<-sc2[sc2$Type=="nonsyn" & sc2$ref=="t",]
T.nonsyn<-ggplot(dt, aes(x=EstSC))+geom_histogram(color="gray60", fill=paste0(colors2[3],"CC"))+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00001,1))+
        theme_bw()+ylab("Count")+xlab("")+ggtitle("T:Nonsyn")+theme(plot.title = element_text(hjust = 0.5, size=12))+
        geom_vline(aes(xintercept=mean(EstSC, na.rm=T)),color="blue", linetype="dashed", size=.5)

dt<-sc2[sc2$Type=="syn" & sc2$ref=="c",]
C.syn<-ggplot(dt, aes(x=EstSC))+geom_histogram(color="gray60", fill=paste0(colors2[5],"CC"))+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00001,1))+
        theme_bw()+ylab("Count")+xlab("")+ggtitle("C:Syn")+theme(plot.title = element_text(hjust = 0.5, size=12))+
        geom_vline(aes(xintercept=mean(EstSC, na.rm=T)),color="blue", linetype="dashed", size=.5)
dt<-sc2[sc2$Type=="nonsyn" & sc2$ref=="c",]
C.nonsyn<-ggplot(dt, aes(x=EstSC))+geom_histogram(color="gray60", fill=paste0(colors2[5],"CC"))+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00001,1))+
        theme_bw()+ylab("Count")+xlab("")+ggtitle("C:Nonsyn")+theme(plot.title = element_text(hjust = 0.5, size=12))+
        geom_vline(aes(xintercept=mean(EstSC, na.rm=T)),color="blue", linetype="dashed", size=.5)

dt<-sc2[sc2$Type=="syn" & sc2$ref=="g",]
G.syn<-ggplot(dt, aes(x=EstSC))+geom_histogram(color="gray60", fill=paste0(colors2[2],"CC"))+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00001,1))+
        theme_bw()+ylab("Count")+xlab("")+ggtitle("G:Syn")+theme(plot.title = element_text(hjust = 0.5, size=12))+
        geom_vline(aes(xintercept=mean(EstSC, na.rm=T)),color="blue", linetype="dashed", size=.5)
dt<-sc2[sc2$Type=="nonsyn" & sc2$ref=="g",]
G.nonsyn<-ggplot(dt, aes(x=EstSC))+geom_histogram(color="gray60", fill=paste0(colors2[2],"CC"))+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00001,1))+
        theme_bw()+ylab("Count")+xlab("")+ggtitle("G:Nonsyn")+theme(plot.title = element_text(hjust = 0.5, size=12))+
        geom_vline(aes(xintercept=mean(EstSC, na.rm=T)),color="blue", linetype="dashed", size=.5)



quartz()
plot_grid(A.syn,T.syn,C.syn,G.syn,A.nonsyn,T.nonsyn, C.nonsyn,G.nonsyn,nrow = 2, ncol=4,
          rel_heights = c(1,1))
ggsave("Output1A/SelCoeff/SC_hist.pdf", width=8, height=5)



######
filename<-paste0("Output1A/SelCoeff/SC_Histograms_median.pdf")
pdf(filename, width = 12, height = 7)
par(mfrow=c(2,4))
par(mar = c(4,4.3,2.5,1))
breaks=seq(-6,0,by=.2)

k=1
for (type in c("syn","nonsyn")){
        for (i in c("a","t","c","g")) {
                datavector<-sc$EstSC[sc$Type==type & sc$ref==i]
                nt<-toupper(i)
                hist(log10(datavector), breaks = breaks, xlim=c(-5,-0),col=colors[k], 
                     main = paste0(nt, " : ", type), ylab="Count", xlab="",
                     xaxt="n", cex.main=1.7, las=1, cex.axis =1,cex.lab=1.4)
                abline(v=median(log10(datavector), na.rm = TRUE),col="white",lwd=2,lty=1)
                abline(v=median(log10(datavector), na.rm = TRUE),col="red",lwd=2,lty=2)
                x=seq(-5,0,by=1)
                labels <- sapply(x,function(y)as.expression(bquote(10^ .(y))))
                axis(1, at=x,labels=labels, col.axis="black", line=-.7, cex.axis=1)
                mtext(text = "Sselection coefficient", side = 1, line = 1.5, cex=.9)
                k=k+1
        }
}

dev.off()

#Same plot with Mean
filename<-paste0("Output1A/SelCoeff/SC_Histograms_mean.pdf")
pdf(filename, width = 12, height = 7)
par(mfrow=c(2,4))
par(mar = c(4,4.3,2.5,1))
breaks=seq(-6,0,by=.2)
k=1
for (type in c("syn","nonsyn")){
        for (i in c("a","t","c","g")) {
                datavector<-sc$EstSC[sc$Type==type & sc$ref==i]
                nt<-toupper(i)
                hist(log10(datavector), breaks = breaks, xlim=c(-5,-0),col=colors[k], 
                     main = paste0(nt, " : ", type), ylab="Count", xlab="",
                     xaxt="n", cex.main=1.7, las=1, cex.axis =1,cex.lab=1.4)
                abline(v=mean(log10(datavector), na.rm = TRUE),col="white",lwd=2,lty=1)
                abline(v=mean(log10(datavector), na.rm = TRUE),col="red",lwd=2,lty=2)
                x=seq(-5,0,by=1)
                labels <- sapply(x,function(y)as.expression(bquote(10^ .(y))))
                axis(1, at=x,labels=labels, col.axis="black", line=-.7, cex.axis=1)
                mtext(text = "Sselection coefficient", side = 1, line = 1.5, cex=.9)
                k=k+1
        }
}
dev.off()

#######
##########################################
#Plot summary of 1) sel coef by type by nucleotide       
#1. 
k=1
transSC<-list()
for (i in c("a","t","c","g")) {
        for (type in c("syn","nonsyn")){
                datavector<-sc2$EstSC[sc2$Type==type & sc2$ref==i]
                nt<-toupper(i)
                vname<-paste0(nt,".",type)
                dat<-data.frame(base=rep(nt,times=length(datavector)),
                                type=rep(type, times=length(datavector)), S.C.=datavector)
                transSC[[k]]<-dat
                names(transSC)[k]<-vname
                k=k+1
        }
}

scdata<-do.call(rbind, transSC)

ggboxplot(scdata, x="base",y="S.C.", color="type", xlab="",
          ylab="Estimated SC", palette = colors2[c(5,1)])+scale_y_continuous(trans = 'log10', labels=label_scientific)+
        geom_vline(xintercept = c(1:3)+0.5, color="gray60")+
        theme(panel.border = element_rect(color="black"))
        

z=rep(c(0.7,0.3),times=4)
x<-1:4
ybreaks<- c(1:10 * 10^c(-4),1:10 * 10^c(-3),1:10 * 10^c(-2),1:10 * 10^c(-1)) 
        
c(10^-4,10^-3,10^-1)
ggplot(scdata,aes(x=base, y=S.C., fill=factor(type)))+geom_boxplot(outlier.alpha =.4, outlier.color = "gray60")+
        scale_y_continuous(trans = 'log10',breaks=c(0.0001,0.001,0.01), minor_breaks=ybreaks, labels=label_scientific2)+
        labs(x="Nucleotide",y="Estimated selection coefficient")+
        scale_fill_manual(values=col2_light[c(5,1)], ) + theme_bw()+
        theme(legend.title = element_blank()) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10), panel.grid.major.x=element_blank())+
        geom_vline(xintercept = c(1:3)+0.5, color="gray60")
ggsave("Output1A/SelCoeff/SC.byNT2.pdf", width = 5,height = 4)

## wilcox test by NT
NT<-c("a","c","t","g")
Ncomb<-t(combn(NT,2))
WilcoxTest.nt<-data.frame(matrix(ncol=4,nrow=nrow(Ncomb)))
colnames(WilcoxTest.nt)<-c("NT1","NT2","test","P.value")

for (i in 1:nrow(Ncomb)) {
        vec1<-sc2$EstSC[sc2$ref==Ncomb[i,1]]
        vec2<-sc2$EstSC[sc2$ref==Ncomb[i,2]]
        result<-wilcox.test(vec1, vec2, alternative = "less", paired = FALSE) 
        
        WilcoxTest.nt$NT1[i]<-Ncomb[i,1]
        WilcoxTest.nt$NT2[i]<-Ncomb[i,2]
        WilcoxTest.nt$test[i]<-"less"
        WilcoxTest.nt$P.value[i]<-result[[3]]
}   

WilcoxTest.nt2<-data.frame(matrix(ncol=4,nrow=nrow(Ncomb)))
colnames(WilcoxTest.nt2)<-c("NT1","NT2","test","P.value")

for (i in 1:nrow(Ncomb)) {
        vec1<-sc2$EstSC[sc2$ref==Ncomb[i,1]]
        vec2<-sc2$EstSC[sc2$ref==Ncomb[i,2]]
        result<-wilcox.test(vec1, vec2, alternative = "greater", paired = FALSE) 
        
        WilcoxTest.nt2$NT1[i]<-Ncomb[i,1]
        WilcoxTest.nt2$NT2[i]<-Ncomb[i,2]
        WilcoxTest.nt2$test[i]<-"greater"
        WilcoxTest.nt2$P.value[i]<-result[[3]]
}   

WilcoxTest.nt<-rbind(WilcoxTest.nt,WilcoxTest.nt2)
write.csv(WilcoxTest.nt,"Output1A/SummaryStats/SC_WilcoxTestResults_byNT.csv")




##########################################
#Plot summary of 1) sel coef _compare CpG creating vs. non-CpGcreating        
#1. Transition SYN
k=1
transSC<-list()
for (i in c("a","t")) {
        for (cp in c(0,1)){
                datavector<-df$EstSC[df$Type=="syn" & df$ref==i &df$makesCpG==cp]
                nt<-toupper(i)
                if (cp==0) cpg<-''
                if (cp==1) cpg<-"(CpG)"
                vname<-paste0(nt,cpg)
                dat<-data.frame(base=rep(vname,times=length(datavector)), sc=datavector)
                transSC[[k]]<-dat
                names(transSC)[k]<-vname
                k=k+1
        }
}

scdata<-do.call(rbind, transSC)

z=c(0.7,0.3,0.7,0.3)
ggplot(scdata,aes(x=base,y=sc,fill=base))+geom_boxplot(alpha=z)+
        scale_y_continuous(trans = 'log10', labels=label_scientific)+
        labs(x="Nucleotide",y="Estimated selection coefficient")+
        scale_fill_manual(values=c("#66CCEEB3","#66CCEE4D","#228833B3","#2288334D")) + theme_classic()+
        guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10))+
        ggtitle("CpG vs. nonCpG syn transition") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
        theme(plot.title = element_text(hjust = 0.5))+
        geom_vline(xintercept = 2.5)

ggsave(filename="Output1A/SelCoeff/SC_CpGvsNonCpG_syn_Ts.pdf",width=4, height=4, units='in',device='pdf')


### Transition Nonsyn

k=1
transSC2<-list()
for (i in c("a","t")) {
                for (cp in c(0,1)){
                        datavector<-df$EstSC[df$Type=="nonsyn" & df$ref==i &df$makesCpG==cp]
                        nt<-toupper(i)
                        if (cp==0) cpg<-''
                        if (cp==1) cpg<-"(CpG)"
                        vname<-paste0(nt,cpg)
                        dat<-data.frame(base=rep(vname,times=length(datavector)), sc=datavector)
                        transSC2[[k]]<-dat
                        names(transSC2)[k]<-vname
                        k=k+1
                        }
        }

scdata2<-do.call(rbind, transSC2)

z=c(0.7,0.3,0.7,0.3)
ggplot(scdata2,aes(x=base,y=sc,fill=base))+geom_boxplot(alpha=z)+labs(x="Nucleotide",y="Estimated selection coefficient")+
        scale_y_continuous(trans = 'log10', labels=label_scientific)+
        scale_fill_manual(values=c("#66CCEEB3","#66CCEE4D","#228833B3","#2288334D")) + theme_classic()+
        guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10))+
        ggtitle("CpG vs. nonCpG nonsyn transition") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
        theme(plot.title = element_text(hjust = 0.5))+
        geom_vline(xintercept = 2.5)

ggsave(filename="Output1A/SelCoeff/SC_CpGvsNonCpG_Nonsyn_Ts.pdf",width=4, height=4, units='in',device='pdf', plot=scplot2)



#####
#dinucleotides
rho(, wordsize = 2, alphabet = s2c("acgt"))
zscore(sequence, simulations = NULL, modele, exact = FALSE, alphabet = s2c("acgt"), ... )
