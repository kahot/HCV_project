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
scColors<-c("#AEE5F6","#4477AA","#EEBAB9","#EE6677","#228833")

colors2<-qualitative_hcl(6, palette="Dark3")
col2_light<-qualitative_hcl(6, palette="Set3")
scColors2<-c("#EC4A4D","#0055EC")



###########

df<-read.csv("Output1A/SelCoeff/SC.csv", stringsAsFactors = F, row.names = 1)
mean(df$EstSC) # 0.002681807  
#coding regions only
df<-df[df$pos>=342,]

#add the gene info
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

endnuc<- sc$pos[nrow(sc)]  # from 11.Mut.Freq calculation. If using the last row -> Trans$pos[nrow(Trnas)]
n<-data.frame("pos"=c(342:endnuc))

SCs<-merge(n,sc,by="pos",all.x=T) 
range(SCs$EstSC,na.rm = T) #0.0001359244 0.0136154784

##
# create a table for syn only and nonsyn only to calculate separate rolling average
ns<-SCs[SCs$Type=="nonsyn",]
syn<-SCs[SCs$Type=="syn",]
ns<-merge(n, ns, by="pos",all.x=T) 
syn<-merge(n, syn, by="pos",all.x=T) 

ns.roll50<-rollmean(ns$EstSC, k=50, na.rm=T, align="center")
#ns.roll100<-rollmean(ns$EstSC, k=100, na.rm=T,align="center")
SCs$ns.roll50<-c(rep(NA, times=25),ns.roll50,c(rep(NA, times=24)))

syn.roll50 <-rollmean(syn$EstSC, k=50, na.rm=T,align="center")
#syn.roll100<-rollmean(syn$EstSC, k=100, na.rm=T,align="center")
SCs$syn.roll50<-c(rep(NA, times=25),syn.roll50,c(rep(NA, times=24)))
#SCs$syn.roll100<-c(rep(NA, times=),syn.roll100,c(rep(NA, times=49)))




### PLOT SCs across the genome
ylow=0.00008
yhigh=0.03
title<-"SC_acroosGenoem2"

pdf(paste0("Output1A/SelCoeff/",title,".pdf"),width=15,height=4.5)
plot(EstSC~pos, data=SCs,t="n",log='y',yaxt='n',xlab='Genome position',ylab="Estimated selection coefficient",
     ylim=c(ylow,yhigh),xlim=c(342,8618))
eaxis(side = 2, at = 10^((0):(-(5))), cex=2)
for(k in 1:5){abline(h = 1:10 * 10^(-k), col = "gray80")}

# roling average of 50
roll50<-rollmean(SCs$EstSC, k=50, na.rm=T)
SCs$roll50<-c(rep(NA, times=25),roll50,c(rep(NA, times=24)))

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

##########
# plot syn and nonsyn with different colors
ylow1=0.00009
yhigh=0.03
title<-"SC_syn.nonsyn"
pdf(paste0("Output1A/SelCoeff/",title,".pdf"),width=14,height=5)
plot(EstSC~pos, data=SCs,t="n",log='y',yaxt='n',xlab='Genome position',ylab="Estimated selection coefficient",
     ylim=c(ylow1,yhigh),xlim=c(265,8618))
eaxis(side = 2, at = 10^((0):(-(5))), cex=2)
for(k in 1:5){abline(h = 1:10 * 10^(-k), col = "gray80")}

for (i in 1:endnuc){
        if (is.na(SCs$Type[i])==T) next
        if (SCs$Type[i]=="stop") next
        if (SCs$Type[i]=="syn") {c=colors2[2]}
        if (SCs$Type[i]=="nonsyn"&SCs$ref[i]%in%c("c","g")) {c=colors2[1]}
        if (SCs$Type[i]=="nonsyn"&SCs$ref[i]%in%c("a","t")) {c=colors2[5]}
        points(SCs$pos[i],SCs$EstSC[i],pch=21,col='gray30',lwd=0.3, bg=paste0(c,"B3"),cex=.4)
}

ylow<-0.00007
for (j in 2:(nrow(genes)-1)){
        xleft<-genes$start[j]
        xright<-genes$start[j+1]
        
        if ((j==4|j==6|j==9)){
                rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="gray40")
                text(xleft+80, 1.44*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                #mtext(paste0(genes$Gene[j]),side= 1, line=-0.1, at= xleft+80, col="black", cex=0.8)
        }
        else if (j==12){
                rect(xleft,ylow,genes$end[j],1.3*ylow,density = NULL, angle = 45,col="white",border ="gray40")
                text(xleft+600,1.15*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        }
        else{rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="gray40")
                text(xright-(xright-xleft)/2,1.15*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)}
}

# roling average of 50
roll50<-rollmean(SCs$EstSC, k=50, na.rm=T)
SCs$roll50<-c(rep(NA, times=25),roll50,c(rep(NA, times=24)))
lines(roll50~pos,data=SCs, col="#4477AA",lwd=1.5)

abline(v=genes$end, col="gray80", lwd=.5)
#Add legend
legpos=300; legposV=0.04
rect(legpos, 0.42*legposV, (legpos+1000), 1.05*legposV, density = NULL, angle = 45,col=alpha("white",1))
points((legpos+100),legposV*0.8,pch=21,bg=colors2[2],col=1,cex=1)
text((legpos+150),legposV*0.8,"Syn",adj=0, cex=1)
points((legpos+100),legposV*0.64,pch=21,bg=colors2[5],col=1,cex=1)
text((legpos+150),legposV*0.64,"Nonsyn, A/T",adj=0, cex=1)
points((legpos+100),legposV*0.51,pch=21,bg=colors2[1],col=1,cex=1)
text((legpos+150),legposV*0.51,"Nonsyn, C/G",adj=0, cex=1)

box()
dev.off()



### plot syn and nonsyn with different colors & separate rollmean lines
ylow1=0.00009
yhigh=0.03
title<-"SC_syn.nonsyn2"

pdf(paste0("Output1A/SelCoeff/",title,".pdf"),width=14,height=5)
plot(EstSC~pos, data=SCs,t="n",log='y',yaxt='n',xlab='Genome position',ylab="Estimated selection coefficient",
     ylim=c(ylow1,yhigh),xlim=c(265,8618))
eaxis(side = 2, at = 10^((0):(-(5))), cex=2)
for(k in 1:5){abline(h = 1:10 * 10^(-k), col = "gray80")}

for (i in 1:endnuc){
        if (is.na(SCs$Type[i])==T) next
        if (SCs$Type[i]=="stop") next
        if (SCs$Type[i]=="syn") {c=colors2[2]}
        if (SCs$Type[i]=="nonsyn"&SCs$ref[i]%in%c("c","g")) {c=colors2[1]}
        if (SCs$Type[i]=="nonsyn"&SCs$ref[i]%in%c("a","t")) {c=colors2[5]}
        points(SCs$pos[i],SCs$EstSC[i],pch=21,col='gray30',lwd=0.3, bg=paste0(c,"B3"),cex=.4)
}

ylow<-0.00007
for (j in 2:(nrow(genes)-1)){
        xleft<-genes$start[j]
        xright<-genes$start[j+1]
        
        if ((j==4|j==6|j==9)){
                rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="gray40")
                text(xleft+80, 1.44*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                #mtext(paste0(genes$Gene[j]),side= 1, line=-0.1, at= xleft+80, col="black", cex=0.8)
        }
        else if (j==12){
                rect(xleft,ylow,genes$end[j],1.3*ylow,density = NULL, angle = 45,col="white",border ="gray40")
                text(xleft+600,1.15*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        }
        else{rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="gray40")
                text(xright-(xright-xleft)/2,1.15*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)}
}

# roling average of 50
lines(ns.roll50~pos,data=SCs, col=scColors2[1],lwd=1)
lines(syn.roll50~pos,data=SCs, col=scColors2[2],lwd=1)

abline(v=genes$end, col="gray80", lwd=.5)

#Add legend
legpos=300; legposV=0.04
rect(legpos, 0.42*legposV, (legpos+900), 1.05*legposV, density = NULL, angle = 45,col=alpha("white",1))
points((legpos+100),legposV*0.8,pch=21,bg=colors2[2],col=1,cex=1)
text((legpos+150),legposV*0.8,"Syn",adj=0, cex=0.8)
points((legpos+100),legposV*0.64,pch=21,bg=colors2[5],col=1,cex=1)
text((legpos+150),legposV*0.64,"Nonsyn, A/T",adj=0, cex=.8)
points((legpos+100),legposV*0.51,pch=21,bg=colors2[1],col=1,cex=1)
text((legpos+150),legposV*0.51,"Nonsyn, C/G",adj=0, cex=.8)



box()
dev.off()


