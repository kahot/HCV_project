library(purrr)
library(tidyverse)
library(zoo)
library(colorspace)
source("Rscripts/label_sceintific.R")

colors2<-qualitative_hcl(6, palette="Dark3")


source("Rscripts/baseRscript.R")

mfs<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F, row.names = 1)
genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector

genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
genenames<-genes$Gene[1:12]

### Plot mutation freq. across the genome based on the mutation types 

n<-data.frame("pos"=c(1:maxpos))
mfs<-merge(n,mfs,by="pos",all.x=T)

pdf("Output1A/SummaryFig.Filtered/MutFreq_by_Types2.pdf",width=15,height=6)
maxnuc=mfs$pos[nrow(mfs)]
par(mar = c(3,5,1,2))
#selcoeffcolumn <-SC3$mean 

plot(mfs$pos[1:maxnuc],mfs$mean[1:maxnuc],
        log="y", ylab="Ave. mutation freq.",cex.lab=1.4,
        yaxt="n", xlab="",xaxt='n',
        col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-4,0.1),xlim=c(342,maxnuc))
axis(1,at=c(seq(500,8500,by=1000)), labels=c(seq(500,8500,by=1000)))
eaxis(side = 2, at = 10^((-1):(-(4))), cex=2)

for(i in 2:4){abline(h = 1:10 * 10^(-i), col = "gray80")}

for (i in 1:maxnuc){
        if (is.na(mfs$Type[i])==T) next
        if (mfs$Type[i]=="stop") {points(mfs$pos[i],mfs$mean[i],pch=21,col='gray30',lwd=0.3, bg="black",cex=.4)
                next}
        if (mfs$Type[i]=="syn") {c=colors2[2]}
        if (mfs$Type[i]=="nonsyn"&mfs$ref[i]%in%c("c","g")) {c=colors2[1]}
        if (mfs$Type[i]=="nonsyn"&mfs$ref[i]%in%c("a","t")) {c=colors2[5]}
        points(mfs$pos[i],mfs$mean[i],pch=21,col='gray30',lwd=0.3, bg=paste0(c,"B3"),cex=.4)
}
        #Add legend
legpos=300; legposV=0.1
rect(legpos, 0.42*legposV, (legpos+1000), 1.05*legposV, density = NULL, angle = 45,col=alpha("white",1))
points((legpos+100),legposV*0.9,pch=21,bg=colors2[2],col=1,cex=1)
text((legpos+150),legposV*0.9,"Syn",adj=0, cex=1)
points((legpos+100),legposV*0.74,pch=21,bg=colors2[5],col=1,cex=1)
text((legpos+150),legposV*0.74,"Non-syn, A/T",adj=0, cex=1)
points((legpos+100),legposV*0.6,pch=21,bg=colors2[1],col=1,cex=1)
text((legpos+150),legposV*0.6,"Non-syn, C/G",adj=0, cex=1)
points((legpos+100),legposV*0.49,pch=21,bg=1,col=1,cex=1)
text((legpos+150),legposV*0.49,"Nonsense",adj=0, cex=1)

ylow<-0.00025
for (j in 2:(nrow(genes)-1)){
        xleft<-genes$start[j]
        xright<-genes$start[j+1]
        
        if ((j==4|j==6|j==9)){
                rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="gray40")
                text(xleft+80, 1.4*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                #mtext(paste0(genes$Gene[j]),side= 1, line=-0.1, at= xleft+80, col="black", cex=0.8)
        }
        else if (j==12){
                rect(xleft,ylow,genes$end[j],1.3*ylow,density = NULL, angle = 45,col="white",border ="gray40")
                text(xleft+600,1.15*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        }
        else{rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="gray40")
                text(xright-(xright-xleft)/2,1.15*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)}
}

roll100<-rollmean(mfs$mean, k=100, na.rm=T)
mfs$roll100<-c(rep(NA, times=50),roll100,c(rep(NA, times=49)))
lines(roll100~pos,data=mfs, col="#001ade", lwd=0.8)
abline(v=genes$end, col="gray80", lwd=.5)

        
box()
dev.off()


####

endnuc<-mfs$pos[nrow(mfs)]
SNPFreq<-mfs

n<-data.frame("pos"=c(342:endnuc))
SNPFreqs<-merge(n,SNPFreq,by="pos",all.x=T)
pdf("Output1A/SummaryFig.Filtered//Ts.MF.acrossGenome.pdf",width=15,height=6.5)
plot(mean~pos, data=SNPFreqs,t="n",log='y',yaxt='n',xlab='Genome position',ylab="Ave.transition mut.freq.",
     ylim=c(0.0002,0.1),xlim=c(340,8500))
eaxis(side = 2, at = 10^((-1):(-(4))), cex=2)
for(i in 2:4){abline(h = 1:10 * 10^(-i), col = "gray80")}
points(mean~pos, data=SNPFreqs,pch=20,col=paste0(colors2[5],"66"),cex=0.5)

#add rolling average
roll100<-rollmean(SNPFreqs$mean, k=100, na.rm=T)
SNPFreqs$roll100<-c(rep(NA, times=50),roll100,c(rep(NA, times=49)))
lines(roll100~pos,data=SNPFreqs, col="blue")


abline(v=genes$end, col="gray80", lwd=.5)

ylow<-0.0002;yhigh<-.3
for (j in 1:nrow(genes)){
        xleft<- genes$start[j]
        xright<-genes$start[j+1]
        if (j==1){
                
                rect(-200,ylow,xright,ylow*1.6,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(150,1.27*ylow,paste0(genes$Gene[j]),col="black", cex=0.7)
        }
        else if (j==6){
                rect(xleft,ylow,xright,ylow*1.6,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xleft+100,1.27*ylow,paste0(genes$Gene[j]),col="black", cex=0.7)
                
        }
        else if (j==4|j==9){
                rect(xleft,ylow,xright,1.6*ylow,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xleft+50,2.1*ylow,paste0(genes$Gene[j]),col="black", cex=0.7)
        }
        else if (j==12){
                rect(xleft,ylow,genes$end[j],1.6*ylow,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xleft+600,1.27*ylow,paste0(genes$Gene[j]),col="black", cex=0.7)
        }
        else{rect(xleft,ylow,xright,1.6*ylow,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xright-(xright-xleft)/2,1.27*ylow,paste0(genes$Gene[j]),col="black", cex=0.7)}
}

box()
dev.off()


## By genes
mf1<-mfs[!is.na(mfs$mean),]
mf1<-merge(mf1, genetable, by="pos", all.x=T )
SumMFGenes<-aggregate(mf1$mean,by=list(mf1$gene),FUN=mean)
SumMF.se<-aggregate(mf1$mean,by=list(mf1$gene),FUN=std.error)

sumG<-cbind(SumMFGenes, SumMF.se$x)
colnames(sumG)<-c("Gene","Mean","SE")
sumG$Gene<-factor(sumG$Gene, levels=c("Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

write.csv(sumG, "Output1A/MutFreq.filtered/MF_Summary_Table.by.gene.csv")

mf2<-mf1[,c("pos","mean","gene")]

ggplot(sumG, aes(x=Gene, y=Mean))+scale_y_continuous(trans="log10", 
          breaks = c(0.001,0.005,0.007,0.009,0.01),labels=c(10^-3,0.005,0.007,'',10^-2))+
        geom_point(data=mf2, aes(x=gene, y=mean), position=position_jitter(width=0.2),stat = "identity", col=paste0(colors2[5],"4D"), size=0.2)+
        geom_point()+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, position=position_dodge(width=0.3))+
        theme_bw()+theme(axis.text=element_text(size=11), axis.title=element_text(size=13))+
        ylab("Average mutation frequency")+
        theme(panel.grid.major.x=element_blank(),axis.title.y = element_text(size=13), panel.grid.minor.y = element_blank())+
        geom_vline(xintercept = c(1:11)+0.5,  
                   color = "gray70", size=.5)+
        theme(axis.title.x=element_blank())
ggsave(filename="Output1A/MutFreq.filtered/Ave.mf.byGene.pdf", width = 8.5, height = 5)


