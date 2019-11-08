library(zoo)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(ggthemes)
library(sfsmisc)
library(colorspace)

source("Rscripts/baseRscript.R")

###########
SC<-list()
fnames<-c("Ts", "Ts_NA", "Ts_zero")
for (i in 1:length(fnames)){
        df<-read.csv(paste0("Output1A/SelCoeff/SC_",fnames[i],"_summary.csv"))
        print(fnames[i])
        print(mean(df$EstSC,na.rm=T))
        SC[[i]]<-df
        names(SC)[i]<-fnames[i]
}
#[1] "Ts"
#[1] 0.002602148
#[1] "Ts_NA"
#[1] 0.002336401
#[1] "Ts_zero"
#[1] 0.002796673      

#### Calcualte the mean SC and GC_contents for each gene
genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)

gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector

Trans<-SC[[1]]

end<-Trans$pos[nrow(Trans)]
genetable<-genetable[genetable$pos>=342&genetable$pos<=end,]

sc<-merge(Trans, genetable, by="pos")
AveSC<-aggregate(sc$EstSC,by=list(sc$gene), mean, na.rm = TRUE)
colnames(AveSC)<-c("gene","SC")
SESC<-aggregate(sc$EstSC,by=list(sc$gene), std.error, na.rm = TRUE)
colnames(SESC)<-c("gene","se")

#   Group.1           x
#1     Core 0.002741075
#2       E1 0.002440335
#3       E2 0.002453301
#4     HVR1 0.001800567
#5  NS1(P7) 0.002223507
#6      NS2 0.002364056
#7      NS3 0.002665575
#8     NS4A 0.002701439
#9     NS4B 0.002638799
#10    NS5A 0.002580037
#11    NS5B 0.002915698


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

SC_summary<-merge(AveSC,GCs,by="gene")
SC_summary<-merge(SC_summary,SESC,by="gene")
write.csv(SC_summary,"Output1A/SummaryStats/SelCoeff_summary_by_genes.csv")

###########  Plot the average Sel Coeff by genes

title<-"Average selection coefficients by genes"

ggboxplot(sc,x="gene",y="EstSC", xlab="Gene", ylab="Estimated selective coefficient",color="gene")+
        theme(legend.position="none")

ggboxplot(sc,x="gene",y="EstSC", xlab="Gene", ylab="Estimated selective coefficient",color="gray40")+
        theme(legend.position="none")

ggsave(filename=paste0("Output1A/SelCoeff/SelCoef_by.gene.Boxplot.pdf"), width = 10, height = 7)



ggplot(SC_summary,aes(gene,y=SC))+
        geom_point(color="royalblue")+
        geom_errorbar(aes(ymin=SC-se, ymax=SC+se), width=.2,color="royalblue" )+
        theme_bw()+
        theme(axis.title.x=element_blank())+ylab("Average selection coefficient")+
        theme(panel.grid.major.x = element_blank())

ggsave(filename="Output1A/SelCoeff/Ave.SC_by.gene.pdf",width =5, height = 4)





### PLOT SCs across the genome
startnuc<-sc$pos[1]
endnuc<- sc$pos[nrow(sc)]  # from 11.Mut.Freq calculation. If using the last row -> Trans$pos[nrow(Trnas)]
n<-data.frame("pos"=c(startnuc:endnuc))

SCs<-merge(n,sc,by="pos",all.x=T) 
range(SCs$EstSC,na.rm = T) #0.0001359244 0.0136154784
ylow=0.00008
yhigh=0.03
cols <- c("#66CCEE","#228833","#CCBB44","#EE6677","#AA3377","#4477AA","#BBBBBB")
title<-"SC_acroosGenoem_twocolors"
pdf(paste0("Output1A/SelCoeff/",title,".pdf"),width=15,height=7.5)
#pdf(paste0("Output1A/SummaryFigures/SC_acrossGenome_twocolors.pdf"),width=15,height=4.5)
plot(EstSC~pos, data=SCs,t="n",log='y',yaxt='n',xlab='Genome position (H77)',ylab="Estimated selection coefficient",
     main=paste0(title),ylim=c(ylow,yhigh),xlim=c(340,8618))
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
                if (j==9) text(xleft+50,.88*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                else text(xleft+200,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                }
        if (j%%2==0){
                rect(xleft,ylow,xright,1.4*ylow,density = NULL, angle = 45,col=scColors[1],border ="gray20",lwd=1.5)
                if (j==4|j==6) text(xleft+50,.88*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                else text(xleft+200,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        }
                
}
       
box()
dev.off()

##Change color for all genes
hcl_palettes(plot = TRUE)

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
pdf(paste0("Output1A/SelCoeff/SC_acrossGenome.pdf"),width=15,height=4.5)
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
colors<-c("#66CCEECC", "#228833CC" ,"#CCBB44CC", "#EE6677CC","#66CCEECC", "#228833CC", "#CCBB44CC", "#EE6677CC")
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















##############################
#### Transversion #####

Tv1<-read.csv("Output/Mut.freq.filtered/Summary_Tv1.MutFreq_filtered.csv")
Tv2<-read.csv("Output/Mut.freq.filtered/Summary_Tv2.MutFreq_filtered.csv")

Tv1<-Tv1[,-1]
Tv2<-Tv2[,-1]


#calculate the Sel Coeff for each site using the mean mutation frequencies (Transveresion mutations)
mutrates<-read.csv("Data/Geller.mutation.rates.csv")
Tv1$TVSmutrate.tv1[Tv1$ref=="a"]<-mean(mutrates$mut.rate[mutrates$mutations=="AC"])
Tv1$TVSmutrate.tv1[Tv1$ref=="c"]<-mean(mutrates$mut.rate[mutrates$mutations=="CA"])
Tv1$TVSmutrate.tv1[Tv1$ref=="g"]<-mean(mutrates$mut.rate[mutrates$mutations=="GC"])
Tv1$TVSmutrate.tv1[Tv1$ref=="t"]<-mean(mutrates$mut.rate[mutrates$mutations=="UA"])

Tv2$TVSmutrate.tv2[Tv2$ref=="a"]<-mean(mutrates$mut.rate[mutrates$mutations=="AU"])
Tv2$TVSmutrate.tv2[Tv2$ref=="c"]<-mean(mutrates$mut.rate[mutrates$mutations=="CG"])
Tv2$TVSmutrate.tv2[Tv2$ref=="g"]<-mean(mutrates$mut.rate[mutrates$mutations=="GU"])
Tv2$TVSmutrate.tv2[Tv2$ref=="t"]<-mean(mutrates$mut.rate[mutrates$mutations=="UG"])



Tv1$EstSC<-""
Tv2$EstSC<-""
for (i in 1:nrow(Tv1)){
        Tv1$EstSC[i] <- EstimatedS(Tv1$TVSmutrate.tv1[i],Tv1[i,'mean'])
        Tv2$EstSC[i] <- EstimatedS(Tv2$TVSmutrate.tv2[i],Tv2[i,'mean'])
}
Tv1$EstSC<-as.numeric(Tv1$EstSC)
Tv2$EstSC<-as.numeric(Tv2$EstSC)
write.csv(Tv1,"Output/SelCoeff/SC_Transv1_summary.csv")
write.csv(Tv2,"Output/SelCoeff/SC_Transv2_summary.csv")

#only the coding regions
TV1<-Tv1[Tv1$pos>=342&Tv1$pos<=8609,]  #remove the last end of the genomes (the same end as Transition mutations)
TV2<-Tv2[Tv2$pos>=342&Tv2$pos<=8609,]  

#sanity check! 
mean(TV1$EstSC,na.rm=T) #[1] 0.004585691 (Ave. sel coefficient of TV1 mutations)
mean(TV2$EstSC,na.rm=T) #[1] 0.004414212 (Ave. sel coefficient of TV2 mutations)


#### Calcualte the mean SC and GC_contents for each gene
if(FLASE){
        genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
        genes[6,1]<-"NS1(P7)"
        gene.vector<-c()
        for (i in 1:(nrow(genes)-1)){
                gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
        }
        genetable<-data.frame("pos"=c(1:length(gene.vector)))
        genetable$gene<-gene.vector
        
        end<-scts$pos[nrow(TV1)]
        genetable<-genetable[genetable$pos>=342&genetable$pos<=end,]
}

sc.tv1<-merge(TV1, genetable, by="pos")
sc.tv2<-merge(TV2, genetable, by="pos")

AveSC.tv1<-aggregate(sc.tv1$EstSC,by=list(sc.t$gene), mean, na.rm = TRUE)
AveSC.tv2<-aggregate(sc.tv2$EstSC,by=list(sc.t$gene), mean, na.rm = TRUE)


SC_TV_summary<-cbind(AveSC.tv1,AveSC.tv2$x)
colnames(SC_TV_summary)<-c("gene","SC.tv1","SC.tv2")
SC_TV_summary<-merge(SC_TV_summary,GCs,by="gene")
colnames(SC_TV_summary)[4]<-c("GCcontent")
write.csv(SC_TV_summary,"Output/SummaryStats/SelCoeff_TV_summary_by_genes.csv")

###########  Plot the average Sel Coeff by genes
title<-"Average selection coefficients by genes (transversions)"

rownames(SC_TV_summary)<-SC_TV_summary[,1]
tv<-SC_TV_summary[,2:3]
tvt<-t(tv)
pdf(paste0("Output/SummaryFigures/SC_tvs_ave_by_gene.pdf"),width=7,height=5.5)
barplot(tvt,beside=T,cex.names=.6)
dev.off()



######## PLOT SCs across the genome
startnuc<-sc.t$pos[1]
endnuc<- 8609  # from 11.Mut.Freq calculation. If using the last row -> Trans$pos[nrow(Trnas)]
n<-data.frame("pos"=c(startnuc:endnuc))

sc1<-merge(n,sc.tv1,by="pos",all.x=T) 
sc2<-merge(n,sc.tv2,by="pos",all.x=T) 

sc1<-sc1[,c(1,2,4,5,11,12,14,15,18,19)]
sc2<-sc2[,c(1,2,4,5,11,12,14,15,18,19)]


scTVs<-cbind(sc1,sc2$EstSC)
colnames(scTVs)[11]<-"EstSC2"
apply(SNPFreq[2:(s+1)], 1, mean, na.rm=T)
scTVs$ave.sc<-apply(scTVs[c("EstSC","EstSC2")], 1,mean,na.rm=T)

range(scTVs$EstSC,na.rm = T) #1.893386e-05 1.353137e-01
ylow=0.00001
yhigh=0.2
cols <- c("#66CCEE","#228833","#CCBB44","#EE6677","#AA3377","#4477AA","#BBBBBB")

title<-"Selection coefficients over the genome (transversion)"
pdf(paste0("Output/SummaryFigures/",title,".pdf"),width=15,height=7.5)
plot(EstSC~pos, data=scTVs,t="n",log='y',yaxt='n',xlab='Genome position (H77)',ylab=paste0(yax),
     main=paste0(title),ylim=c(ylow,yhigh),xlim=c(340,8618))
eaxis(side = 2, at = 10^((0):(-(5))), cex=2)
for(k in 1:5){abline(h = 1:10 * 10^(-k), col = "gray80")}

#points(sc1$EstSC~sc1$pos, pch=20,col=scColors[1],cex=0.3)
#points(sc2$EstSC~sc2$pos, pch=20,col=scColors[1],cex=0.3)

#add rolling average
# roling average of 50
roll50<-rollmean(scTVs$ave.sc, k=50, na.rm=T)
scTVs$roll50<-c(rep(NA, times=25),roll50,c(rep(NA, times=24)))

scColors<-c("#AEE5F6","#4477AA","#8CD0C8","#009988","#228833")
for (i in 2:(nrow(genes)-1)){
        gname<-genes$Gene[i]
        if (i%%2==0) {
                points(sc1$EstSC[sc1$gene==gname]~sc1$pos[sc1$gene==gname], pch=20,col=scColors[1],cex=0.3)
                points(sc2$EstSC[sc2$gene==gname]~sc2$pos[sc2$gene==gname], pch=20,col=scColors[1],cex=0.3)
                lines(scTVs$roll50[sc1$gene==gname]~sc1$pos[sc1$gene==gname], col=scColors[2],lwd=1.5)
                }
        if (i%%2==1) {
                points(sc1$EstSC[sc1$gene==gname]~sc1$pos[sc1$gene==gname], pch=20,col=scColors[3],cex=0.5)
                points(sc2$EstSC[sc2$gene==gname]~sc2$pos[sc2$gene==gname], pch=20,col=scColors[3],cex=0.5)
                lines(scTVs$roll50[sc2$gene==gname]~sc2$pos[sc2$gene==gname], col=scColors[4],lwd=1.5)
                }
}

for (j in 1:(nrow(genes)-1)){
        xleft<-genes$start[j]
        xright<-genes$start[j+1]
        if (j%%2==1){
                rect(xleft,ylow,xright,ylow*1.4,density = NULL, angle = 45,col=scColors[3],border ="gray20",lwd=1.5)
                if (j==1) text(200,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                if (j==9) text(xleft+50,.9*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                else text(xleft+200,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        }
        if (j%%2==0){
                rect(xleft,ylow,xright,1.4*ylow,density = NULL, angle = 45,col=scColors[1],border ="gray20",lwd=1.5)
                if (j==4|j==6) text(xleft+50,.9*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                else text(xleft+200,1.2*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        }
        
}

box()
dev.off()






##### 2) Histogram of tranversion SCs


filename<-paste0("Output/SelCoeff/filtered/Tranversion_SelCoef_Histograms_mean.pdf")
pdf(filename, width = 12, height = 7)
par(mfrow=c(2,4))
par(mar = c(4,4.3,2.5,1))
breaks=seq(-6,0,by=.2)

k=1
for (type in c("syn","nonsyn")){
        for (i in c("a","t","c","g")) {
                datavector1<-sc1$EstSC[sc1$Type.tv1==type&sc1$ref==i]
                datavector2<-sc2$EstSC[sc2$Type.tv2==type&sc2$ref==i]
                datavector1<-c(datavector1,datavector2)
                nt<-toupper(i)
                hist(log10(datavector1), breaks = breaks, xlim=c(-5,-0),col=colors[k], 
                     main = paste0(nt, " : ", type), ylab="Count", xlab="",
                     xaxt="n", cex.main=1.7, las=1, cex.axis =1,cex.lab=1.4)
                abline(v=mean(log10(datavector1), na.rm = TRUE),col="white",lwd=2,lty=1)
                abline(v=mean(log10(datavector1), na.rm = TRUE),col="red",lwd=2,lty=2)
                x=seq(-5,0,by=1)
                labels <- sapply(x,function(y)as.expression(bquote(10^ .(y))))
                axis(1, at=x,labels=labels, col.axis="black", line=-.7, cex.axis=1)
                mtext(text = "Sselection coefficient", side = 1, line = 1.5, cex=.9)
                k=k+1
        }
}

dev.off()


filename<-paste0("Output/SelCoeff/filtered/Tranversion_SelCoeff_Histograms_median.pdf")
pdf(filename, width = 12, height = 7)
par(mfrow=c(2,4))
par(mar = c(4,4.3,2.5,1))
breaks=seq(-6,0,by=.2)
k=1
for (type in c("syn","nonsyn")){
        for (i in c("a","t","c","g")) {
                datavector1<-sc1$EstSC[sc1$Type.tv1==type&sc1$ref==i]
                datavector2<-sc2$EstSC[sc2$Type.tv2==type&sc2$ref==i]
                datavector1<-c(datavector1,datavector2)
                nt<-toupper(i)
                hist(log10(datavector1), breaks = breaks, xlim=c(-5,-0),col=colors[k], 
                     main = paste0(nt, " : ", type), ylab="Count", xlab="",
                     xaxt="n", cex.main=1.7, las=1, cex.axis =1,cex.lab=1.4)
                abline(v=median(log10(datavector1), na.rm = TRUE),col="white",lwd=2,lty=1)
                abline(v=median(log10(datavector1), na.rm = TRUE),col="red",lwd=2,lty=2)
                x=seq(-5,0,by=1)
                labels <- sapply(x,function(y)as.expression(bquote(10^ .(y))))
                axis(1, at=x,labels=labels, col.axis="black", line=-.7, cex.axis=1)
                mtext(text = "Sselection coefficient", side = 1, line = 1.5, cex=.9)
                k=k+1
}}
dev.off()



##########################################
#Plot summary of 1) sel coef _compare CpG creating vs. non-CpGcreating        
#1. Transition SYN
k=1
transSC<-list()
for (i in c("a","t")) {
        for (cp in c(0,1)){
                datavector<-sc.t$EstSC[sc.t$Type=="syn" & sc.t$ref==i &sc.t$makesCpG==cp]
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
scplot1<-ggplot(scdata,aes(x=base,y=sc,fill=base))+geom_boxplot(alpha=z)+labs(x="Nucleotide",y="Estimated selection coefficient")+
        scale_fill_manual(values=c("#66CCEEB3","#66CCEE4D","#228833B3","#2288334D")) + theme_classic()+
        guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10))+
        ggtitle("CpG vs. nonCpG syn transition") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
        theme(plot.title = element_text(hjust = 0.5))+
        geom_vline(xintercept = 2.5)

ggsave(filename="Output/SelCoeff/filtered/CpGvsNonCpG_syn_Transition.pdf",width=4, height=4, units='in',device='pdf', plot=scplot1)


### Transition Nonsyn

k=1
transSC2<-list()
for (i in c("a","t")) {
                for (cp in c(0,1)){
                        datavector<-sc.t$EstSC[sc.t$Type=="nonsyn" & sc.t$ref==i &sc.t$makesCpG==cp]
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
scplot2<-ggplot(scdata2,aes(x=base,y=sc,fill=base))+geom_boxplot(alpha=z)+labs(x="Nucleotide",y="Estimated selection coefficient")+
        scale_fill_manual(values=c("#66CCEEB3","#66CCEE4D","#228833B3","#2288334D")) + theme_classic()+
        guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10))+
        ggtitle("CpG vs. nonCpG nonsyn transition") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
        theme(plot.title = element_text(hjust = 0.5))+
        geom_vline(xintercept = 2.5)

ggsave(filename="Output/SelCoeff/filtered/CpGvsNonCpG_Nonsyn_Transition.pdf",width=4, height=4, units='in',device='pdf', plot=scplot2)


###
# 2.Transversion syn     

k=1
tvSC<-list()
for (i in c("a","t","c","g")) {
        for (cp in c(0,1)){
                datavector1<-sc.tv1$EstSC[sc.tv1$Type.tv1=="syn" & sc.tv1$ref==i &sc.tv1$makesCpG.tv1==cp]
                datavector2<-sc.tv2$EstSC[sc.tv2$Type.tv2=="syn" & sc.tv2$ref==i &sc.tv2$makesCpG.tv2==cp]
                datavec<-c(datavector1,datavector2)
                
                nt<-toupper(i)
                if (cp==0) cpg<-''
                if (cp==1) cpg<-"(CpG)"
                vname<-paste0(nt,cpg)
                dat<-data.frame(base=rep(vname,times=length(datavec)), sc=datavec)
                tvSC[[k]]<-dat
                names(tvSC)[k]<-vname
                k=k+1
        }
}
svdata<-do.call(rbind, tvSC)

z=c(0.7,0.3,0.7,0.3,0.7,0.3,0.7,0.3)
scplot3<-ggplot(svdata,aes(x=base,y=sc,fill=base))+geom_boxplot(alpha=z)+labs(x="Nucleotide",y="Estimated selection coefficient")+
        scale_fill_manual(values=c("#66CCEE","#66CCEE","#228833","#228833","#CCBB44","#CCBB44","#EE6677","#EE6677")) + theme_classic()+
        guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10))+
        ggtitle("CpG vs. nonCpG syn transversion") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
        theme(plot.title = element_text(hjust = 0.5))+
        geom_vline(xintercept = c(2.5,4.5,6.5), color="gray60")

ggsave(filename="Output/SelCoeff/filtered/CpGvsNonCpG_syn_Transv.pdf",width=6, height=4, units='in',device='pdf', plot=scplot3)

#nonsyn transversion
k=1
tvSC2<-list()
for (i in c("a","t","c","g")) {
        for (cp in c(0,1)){
                datavector1<-sc.tv1$EstSC[sc.tv1$Type.tv1=="nonsyn" & sc.tv1$ref==i &sc.tv1$makesCpG.tv1==cp]
                datavector2<-sc.tv2$EstSC[sc.tv2$Type.tv2=="nonsyn" & sc.tv2$ref==i &sc.tv2$makesCpG.tv2==cp]
                datavec<-c(datavector1,datavector2)
                
                nt<-toupper(i)
                if (cp==0) cpg<-''
                if (cp==1) cpg<-"(CpG)"
                vname<-paste0(nt,cpg)
                dat<-data.frame(base=rep(vname,times=length(datavec)), sc=datavec)
                tvSC2[[k]]<-dat
                names(tvSC2)[k]<-vname
                k=k+1
        }
}
svdata2<-do.call(rbind, tvSC2)

z=c(0.7,0.3,0.7,0.3,0.7,0.3,0.7,0.3)
scplot4<-ggplot(svdata2,aes(x=base,y=sc,fill=base))+geom_boxplot(alpha=z)+labs(x="Nucleotide",y="Estimated selection coefficient")+
        scale_fill_manual(values=c("#66CCEE","#66CCEE","#228833","#228833","#CCBB44","#CCBB44","#EE6677","#EE6677")) + theme_classic()+
        guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10))+
        ggtitle("CpG vs. nonCpG nonsyn transversion") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
        theme(plot.title = element_text(hjust = 0.5))+
        geom_vline(xintercept = c(2.5,4.5,6.5), color="gray60")

ggsave(filename="Output/SelCoeff/filtered/CpGvsNonCpG_Nonsyn_Transv.pdf",width=6, height=4, units='in',device='pdf', plot=scplot4)



