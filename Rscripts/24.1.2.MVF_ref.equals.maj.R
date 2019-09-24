library(stringr)
library(tidyverse)
library(zoo)
library(reshape2)
library(colorspace)
library(RColorBrewer)
library(sfsmisc)

source("Rscripts/baseRscript.R")

#Geonotype colors 
colors2<-qualitative_hcl(6, palette="Dark3")
# 1A=1, 1B=3, 3A=5, c("#E16A86","#50A315","#009ADE")

geno<-c("1A","1B","3A")

sameDF<-read.csv("Output_all/Unfiltered/Ts.Same_Mean_3genotypes.csv", stringsAsFactors = F, row.names = 1)
sameDF$gene[sameDF$gene=="NS1(P7)"]<-"NS1"
sameDF$gene<-factor(sameDF$gene, levels=c("5' UTR","Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))


genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
genenames<-genes$Gene[1:12]
#############################################################
sameDF$roll1001a<-c(rep(NA, times=99), rollmean(sameDF$mean.1A, k=100, na.rm=T))
sameDF$roll1001b<-c(rep(NA, times=99), rollmean(sameDF$mean.1B, k=100, na.rm=T))
sameDF$roll1003a<-c(rep(NA, times=99), rollmean(sameDF$mean.3A, k=100, na.rm=T))

pdf(paste0("Output_all/REFequalsMAJ/MVF_AcrossGenome.pdf"),width=12,height=4.5)

plot(mean.1A~merged.pos, data=sameDF,t="n",log='y',yaxt='n',xlab='Genome position', ylab="Minor variant frequency",
     ylim=c(10^-4,0.4),xlim=c(264,8500))
eaxis(side = 2, at = 10^((-0):(-(4))), cex=2)
for(k in 1:4){abline(h = 1:10 * 10^(-k), col = "gray90",lty=1,lwd=.5)}

points(mean.1B~merged.pos,pch=16, data=sameDF, col=paste0(colors2[3],"B3"), cex=0.2)
points(mean.3A~merged.pos,pch=16, data=sameDF, col=paste0(colors2[5],"B3"), cex=0.2)
points(mean.1A~merged.pos,pch=16, data=sameDF, col=paste0(colors2[1],"B3"), cex=0.2)

lines(roll1001b~merged.pos,data=sameDF, col=colors2[3],lwd=1)
lines(roll1003a~merged.pos,data=sameDF, col=colors2[5],lwd=1)
lines(roll1001a~merged.pos,data=sameDF, col=colors2[1],lwd=1)

legend("topright",legend=c("1A","1B","3A"), col=colors2[c(1,3,5)], pch=16, bty = "n", cex=1)

genes2<-genes[1:12,]
#genes2$Gene[6]<-"NS1"

abline(v=genes2$end, col="gray80", lwd=.5)

ylow<-0.0001;yhigh<-.4
for (j in 1:nrow(genes2)){
        xleft<- genes2$start[j]
        xright<-genes2$start[j+1]
        if (j==1){
                
                rect(-200,ylow,xright,ylow*1.6,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(150,1.27*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else if (j==6){
                rect(xleft,ylow,xright,ylow*1.6,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xleft+100,1.27*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
                
        }
        else if (j==4|j==9){
                rect(xleft,ylow,xright,1.6*ylow,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xleft+50,2.1*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else if (j==12){
                rect(xleft,ylow,genes2$end[j],1.6*ylow,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xleft+600,1.27*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else{rect(xleft,ylow,xright,1.6*ylow,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xright-(xright-xleft)/2,1.27*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)}
}

box()
dev.off()


####


#Tvs
Sum2<-read.csv(paste0("Output_all/Unfiltered/Tvs.same_Mean_metadata.1A.csv"),stringsAsFactors = F, row.names = 1)
colnames(Sum2)[ncol(Sum2)]<-paste0("mean.tvs.1A")
for (f in 2:3){
        d<-read.csv(paste0("Output_all/Unfiltered/Tvs.same_Mean_metadata.",geno[f],".csv"),stringsAsFactors = F)
        d<-d[,c(2,ncol(d))]
        colnames(d)[2]<-paste0("mean.tvs.",geno[f])
        Sum2<-merge(Sum2,d,by="merged.pos")
}
write.csv(Sum2,"Output_all/Unfiltered/Tvs.same_Mean_3genotypes.csv")

#all
Sum3<-read.csv(paste0("Output_all/Unfiltered/All.same_Mean_metadata.1A.csv"),stringsAsFactors = F, row.names = 1)
colnames(Sum3)[ncol(Sum3)]<-paste0("mean.all.1A")
for (f in 2:3){
        d<-read.csv(paste0("Output_all/Unfiltered/All.same_Mean_metadata.",geno[f],".csv"),stringsAsFactors = F)
        d<-d[,c(2,ncol(d))]
        colnames(d)[2]<-paste0("mean.all.",geno[f])
        Sum3<-merge(Sum3,d,by="merged.pos")
}
write.csv(Sum3,"Output_all/Unfiltered/All.same_Mean_3genotypes.csv")

sum1<-sameDF[,c("merged.pos", "mean.1A","mean.1B", "mean.3A")]
sum2<-Sum2[,c("merged.pos","mean.tvs.1A", "mean.tvs.1B", "mean.tvs.3A")]
sum3<-Sum3[,c("merged.pos","mean.all.1A", "mean.all.1B", "mean.all.3A")]

SumT<-merge(sum1, sum2, by="merged.pos")
SumT<-merge(SumT, sum3, by="merged.pos")



tb2<-data.frame(type=colnames(SumT)[2:10])
tb2$Genotype<-rep(c("1A","1B","3A"),times=3)
tb2$Type<-c(rep("Transition",times=3), rep("Transversion",times=3), rep("Total MV",times=3))
for (i in 1:nrow(tb2)){
        n<-i+1
        tb2$Mean[i]<-mean(SumT[,n], na.rm=T)
        tb2$SE[i]<-std.error(SumT[,n], na.rm=T)
}
write.csv(tb2, "Output_all/REFequalsMAJ/MVFreq_Summary_3genotypes.csv")


ggplot(tb2,aes(x=Type,y=Mean,group=Genotype, color=Genotype))+
        geom_point(position=position_dodge(width=0.3),size =1.5)+scale_color_manual(values=colors2[c(1,3,5)])+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, position=position_dodge(width=0.3))+
        theme(axis.title.x=element_blank())+ylab("Mean minor variant frequency")+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=13),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:2)+0.5,  
                   color = "gray60", size=.5)
        
ggsave(filename="Output_all/REFequalsMAJ/MVF.by.type_by.genotypes.pdf",width = 5, height = 5.8)



#barplot
ggplot(tb2,aes(x=Type,y=Mean, ymin=Mean-SE, ymax=Mean+SE, fill=Genotype))+
        geom_bar(position=position_dodge(.9), stat="identity",width=0.8)+
        scale_fill_manual(values=paste0(colors2[c(1,3,5)],"E6"), labels=c("1A","1B","3A"))+
        geom_errorbar(position=position_dodge(.9), width=.2, color="gray30")+
        theme(axis.title.x=element_blank())+ylab("Mean minor variant frequency")+
        theme_linedraw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        geom_vline(xintercept = c(1:2)+0.5,  
                   color = "gray60", size=.4)+
        theme(panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(linetype=2, colour="gray60"),
              panel.grid.minor.y = element_blank())
ggsave(filename="Output_all/REFequalsMAJ/MVF_by.type_by.genotypes_barplot2.pdf",width = 6, height = 4.5)



####
## plot mean and SE for each gene
mf<-data.frame()
#geno<-c("1A","1B","3A")
for(g in 1:3){
        dat<-sameDF[,c("merged.pos","gene",paste0("mean.",geno[g]))]
        colnames(dat)[3]<-"mean"
        dat$genotype<-geno[g]
        mf<-rbind(mf,dat)
        }


mf1<-mf[!is.na(mf$mean),]
SumMFGenes<-aggregate(mf1$mean,by=list(mf1$genotype,mf1$gene),FUN=mean)
SumMF.se<-aggregate(mf1$mean,by=list(mf1$genotype,mf1$gene),FUN=std.error)

sumG<-cbind(SumMFGenes, SumMF.se$x)
colnames(sumG)<-c("Genotype","Gene","Mean","SE")
sumG$Gene<-factor(sumG$Gene, levels=c("5' UTR","Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))


ggplot(sumG, aes(x=Gene, y=Mean, group=Genotype, color=Genotype))+
        geom_point(position=position_dodge(width=0.3))+scale_fill_manual(values=colors2[c(1,3,5)])+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, position=position_dodge(width=0.3))+
        theme_bw()+theme(axis.text=element_text(size=11), axis.title=element_text(size=13))+ylab("Mean mutation frequency")+
        theme(panel.grid.major.x=element_blank(),axis.title.y = element_text(size=13))+
        geom_vline(xintercept = c(1:11)+0.5,  
                   color = "gray70", size=.5)
ggsave(filename="Output_all/REFequalsMAJ/MVF_by.gene_by.genotype.pdf", width = 8.5, height = 5)



ggplot(sumG,aes(x=Gene,y=Mean, ymin=Mean-SE, ymax=Mean+SE, fill=Genotype))+
        geom_bar(position=position_dodge(.9), stat="identity",width=0.8)+
        scale_fill_manual(values=paste0(colors2[c(1,3,5)],"E6"), labels=c("1A","1B","3A"))+
        geom_errorbar(position=position_dodge(.9), width=.2, color="gray40")+
        theme(axis.title.x=element_blank())+ylab("Mean minor variant frequency")+
        theme_linedraw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        geom_vline(xintercept = c(1:11)+0.5,  
                   color = "gray60", size=.4)+
        theme(panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(linetype=2, colour="gray60"),
              panel.grid.minor.y = element_blank())

ggsave(filename="Output_all/REFequalsMAJ/MVF_by.gene_by.genotype_Barplot.pdf", width = 8.5, height = 5)



### Plot MVF separately for Nonsyn and syn

mf2<-data.frame()
for(g in 1:3){
        dat<-sameDF[,c("merged.pos","gene",paste0("mean.",geno[g]),paste0("Type.",geno[g]))]
        colnames(dat)[3:4]<-c("mean", "Type")
        dat$genotype<-geno[g]
        mf2<-rbind(mf2,dat)
}


mf2<-mf2[!is.na(mf2$mean),]
#remove the stop sites
mfstop<-mf2[mf2$Type=="stop",]
mf3<-mf2
#plot(mfstop$mean[mfstop$genotype=="1A"])
mf2<-mf2[mf2$Type!="stop",]
#remove 5'UTR
mf2<-mf2[mf2$gene!="5' UTR",]
SumMFGenes2<-aggregate(mf2$mean,by=list(mf2$genotype,mf2$gene, mf2$Type),FUN=mean)
SumMF.se2<-aggregate(mf2$mean,by=list(mf2$genotype,mf2$gene, mf2$Type),FUN=std.error)

sumG2<-cbind(SumMFGenes2, SumMF.se2$x)
colnames(sumG2)<-c("Genotype","Gene","Type","Mean","SE")
sumG2$Gene<-factor(sumG2$Gene, levels=c("Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))


ggplot(sumG2, aes(x=Gene, y=Mean, shape=Type,color=Genotype))+
        geom_point(position=position_dodge(width=0.8),size =2)+scale_color_manual(values=colors2[c(1,3,5)])+
        scale_shape_manual(values=c(4,19))+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, size=.2, position=position_dodge(width=0.8))+
        theme(axis.title.x=element_blank())+ylab("Minor variant frequency")+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray60", size=.4)+
        guides(
                shape = guide_legend(order = 2),
                color=guide_legend(order=1)
        )

ggsave(filename="Output_all/REFequalsMAJ/MVF_by.gene_by.genotype_S.NS.pdf", width = 8, height = 5)


#### include nonsense (stop) mutations
#remove 5' UTR
mf3.2<-mf3[mf3$gene!="5' UTR",]
SumMFGenes3<-aggregate(mf3.2$mean,by=list(mf3.2$genotype,mf3.2$gene, mf3.2$Type),FUN=mean)
SumMF.se3<-aggregate(mf3.2$mean,by=list(mf3.2$genotype,mf3.2$gene, mf3.2$Type),FUN=std.error)

sumG3<-cbind(SumMFGenes3, SumMF.se3$x)
colnames(sumG3)<-c("Genotype","Gene","Type","Mean","SE")
sumG3$Gene<-factor(sumG3$Gene, levels=c("5' UTR","Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))


ggplot(sumG3, aes(x=Gene, y=Mean, shape=Type,color=Genotype))+
        geom_point(position=position_dodge(width=0.8),size =2)+scale_color_manual(values=colors2[c(1,3,5)])+
        scale_shape_manual(values=c(4,17,19))+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, size=.2, position=position_dodge(width=0.8))+
        theme(axis.title.x=element_blank())+ylab("Minor variant frequency")+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray60", size=.4)+
        guides(
                shape = guide_legend(order = 2),
                color=guide_legend(order=1)
        )

ggsave(filename="Output_all/REFequalsMAJ/MVF_by.gene_by.genotype_S.NS.stop.pdf", width = 8.5, height = 5)





### 
# stat test
Comb<-t(combn(1:3, 2))
combnames<-c("1A-1B","1A-3A","1B-3A")

results.Ts<-list()
results.Tvs<-list()
results.all<-list()

for (i in 1:3){
        k=1
        n1<-Comb[i,1]+k
        n2<-Comb[i,2]+k
        r1<-wilcox.test(SumT[,n1], SumT[,n2], alternative = "greater", paired = FALSE) 
        results.Ts[[i]]<-r1
        names(results.Ts)[i]<-combnames[i]
        k=4
        n1<-Comb[i,1]+k
        n2<-Comb[i,2]+k
        r2<-wilcox.test(SumT[,n1], SumT[,n2], alternative = "greater", paired = FALSE) 
        results.Tvs[[i]]<-r2
        names(results.Tvs)[i]<-combnames[i]
        k=7
        n1<-Comb[i,1]+k
        n2<-Comb[i,2]+k
        r3<-wilcox.test(SumT[,n1], SumT[,n2], alternative = "greater", paired = FALSE) 
        results.all[[i]]<-r3
        names(results.all)[i]<-combnames[i]
 
}

sink("Output_all/REFequalsMAJ/MutFreq_mean_comparison_genotypes.txt")
for (i in 1:3){
        print(names(results.Ts[i]))
        print("Transition")
        print(results.Ts[[i]][3])
        print("Transversion")
        print(results.Tvs[[i]][3])
        print("Total")
        print(results.all[[i]][3])
}
sink(NULL)

#[1] "1A-1B,  
#Transition,   p.value = 5.534574e-12
#Transversion $p.value = 9.023772e-19
# "Total"     $p.value = 4.7754e-15
#
#[1] "1A-3A"
#"Transition   $p.value = 7.458138e-28
#"Transversion"$p.value = 1.374172e-101
#"Total"       $p.value = 2.327866e-44

#[1] "1B-3A"
#"Transition"  $p.value = 9.593004e-06
#"Transversion"$p.value = 2.768475e-42
#"Total"       $p.value = 2.498357e-11


sumTs<-sameDF[,c("merged.pos", "gene", "mean.1A","mean.1B", "mean.3A")]
genenames<-unique(genes$Gene)
Comb<-t(combn(1:3, 2))
combnames<-c("1A-1B","1A-3A","1B-3A")

wilcox.res<-data.frame(matrix(nrow=12, ncol=3))
rownames(wilcox.res)<-genenames[1:12]
colnames(wilcox.res)<-combnames                       
for (i in 1:3){
        k=2
        n1<-Comb[i,1]+k
        n2<-Comb[i,2]+k
        for (g in 1:12){
                r1<-wilcox.test(sumTs[sumTs$gene==genenames[g],n1], sumTs[sumTs$gene==genenames[g],n2], alternative = "greater", paired = FALSE) 
                wilcox.res[g,i]<-r1[[3]][1]
        }
}

write.csv(wilcox.res,"Output_all/REFequalsMAJ/MVF_wilcoxTest_byGene.csv",row.names = T)


### 
# stat test for comparing nonsyn vs syn MVF 
genenames<-unique(genes$Gene)

sumTs<-sameDF[,c("merged.pos", "gene", "mean.1A","mean.1B", "mean.3A","Type.1A","Type.1B","Type.3A")]

wilcox.resS<-data.frame(matrix(nrow=12, ncol=3))
wilcox.resNS<-data.frame(matrix(nrow=12, ncol=3))
wilcox.res2<-data.frame(matrix(nrow=12, ncol=3))
rownames(wilcox.resS)<-genenames
colnames(wilcox.resS)<- combnames                 
rownames(wilcox.resNS)<-genenames
colnames(wilcox.resNS)<- combnames                   
rownames(wilcox.res2)<-genenames
colnames(wilcox.res2)<- geno                    

for (i in 1:3){
        d<-sumTs[,c("merged.pos", "gene", paste0("mean.",geno[i]),paste0("Type.",geno[i]))]
        colnames(d)[3:4]<-c("mean", "Type")
        
        k=2
        n1<-Comb[i,1]+k
        n2<-Comb[i,2]+k
        
        if (i==1) d1<-sumTs[sumTs$Type.1A==sumTs$Type.1B,]
        if (i==2) d1<-sumTs[sumTs$Type.1A==sumTs$Type.3A,]
        if (i==3) d1<-sumTs[sumTs$Type.1B==sumTs$Type.3A,]
        
        
        for (g in 1:12){
                r1s<-wilcox.test(d1[d1$gene==genenames[g]& d1$Type.1A=="syn",n1], d1[d1$gene==genenames[g] & d1$Type.1A=="syn",n2], alternative = "greater", paired = FALSE) 
                wilcox.resS[g,i]<-r1s[[3]][1]
                r1n<-wilcox.test(d1[d1$gene==genenames[g]& d1$Type.1A=="nonsyn",n1], d1[d1$gene==genenames[g] & d1$Type.1A=="nonsyn",n2], alternative = "greater", paired = FALSE) 
                wilcox.resNS[g,i]<-r1n[[3]][1]
                
                r2<-wilcox.test(d$mean[d$gene==genenames[g]& d$Type=="syn"], d$mean[d$gene==genenames[g] & d$Type=="nonsyn"], alternative = "greater", paired = FALSE) 
                wilcox.res2[g,i]<-r2[[3]][1]
                
                
                
        }
}

write.csv(wilcox.resS, "Output_all/REFequalsMAJ/MVF_wilcoxTest_byGene.syn.csv",row.names = T)
write.csv(wilcox.resNS,"Output_all/REFequalsMAJ/MVF_wilcoxTest_byGene.ns.csv",row.names = T)
write.csv(wilcox.res2, "Output_all/REFequalsMAJ/MVF_wilcoxTest_nonsyn.vs.syn.csv",row.names = T)


