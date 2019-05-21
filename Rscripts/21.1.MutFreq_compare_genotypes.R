library(ggpubr)
library(zoo)
library(sfsmisc)
library(colorspace)
library(plotrix)
library(gplots)


geno<-c("1A","1B","3A")
cols2<-c("#66CCEE","#EE6677" ,"#228833")
colors<-c("#66CCEE99","#EE667799" ,"#22883399")

#read mut freq summary data
for (i in 1:3){
        dat<-read.csv(paste0("Output_all/Ts_summary_metadata.",geno[i],".csv"),stringsAsFactors = F)
        dat<-dat[c(1:8650),-1]
        assign(paste0("geno",geno[i]),dat)
}

#create a mut. freq by gene plot by genotype
for (g in 1:3){
        dat<-get(paste0("geno",geno[g]))
        plot1<-ggboxplot(dat,x="gene",y="mean", main=paste0(geno[g]), ylim=c(0,0.04),xlab="Gene", 
                         ylab="Average mutation frequency",color="gene")+
                         theme(legend.position="none")
        ggsave(filename=paste0("Output_all/Mutfreq_byGene.",geno[g],".pdf"),plot=plot1, width = 10, height = 7)
}

#create a plot with all genotypes
mf<-data.frame()
for (g in 1:3){
        dat<-get(paste0("geno",geno[g]))
        dt<-dat[,c("merged.pos","gene","mean")]
        dt$genotype<-geno[g]
        mf<-rbind(mf,dt)
}

plot2<-ggboxplot(mf,x="gene",y="mean", color="genotype", ylim=c(0,0.04),xlab="Gene",
                  ylab="Average mutation frequency",palette = colors)
ggsave(filename=paste0("Output_all/Mutfreq_byGene_byGenotypes.pdf"),plot=plot2, width = 14, height = 8)


######## Plot Mut Freq across the genome by genes to compare 3 genotypes
colors2<-paste0(cols2,"66")

genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
genenames<-genes$Gene[1:12]

pdf(paste0("Output_all/Mutfreq.50.by_gene.pdf"),width=12,height=18)
par(mfrow = c(6,2))
for (i in 1:length(genenames)){
        plotname<-genenames[i]
        plot(geno1A$mean[geno1A$gene==plotname], pch=16, cex=.7,col=colors2[1], main=plotname, ylab="Mutation frequency", xlab="Genome position")
        legend("topright",legend=c("1A","1B","3A"), col=colors[c(1,2,3)], pch=16, bty = "n")
        points(geno1B$mean[geno1B$gene==plotname], pch=16, cex=.7,col=colors2[2])
        points(geno3A$mean[geno3A$gene==plotname], pch=16, cex=.7,col=colors2[3])
        
        roll501a<-rollmean(geno1A$mean[geno1A$gene==plotname], k=50, na.rm=T)
        roll501b<-rollmean(geno1B$mean[geno1B$gene==plotname], k=50, na.rm=T)
        roll503a<-rollmean(geno3A$mean[geno3A$gene==plotname], k=50, na.rm=T)
        
        roll501a<-c(rep(NA, times=50),roll501a)
        roll501b<-c(rep(NA, times=50),roll501b)
        roll503a<-c(rep(NA, times=50),roll503a)
        
        lines(roll501a,  col=cols2[1],lwd=2)
        lines(roll501b,  col=cols2[2],lwd=2)
        lines(roll503a,  col=cols2[3],lwd=2)
}
dev.off()


#####
# plot rolling average of 3 genotypes together across genome
pdf(paste0("Output_all/RollAve.100_all.pdf"),width=12,height=6)
plot(mean~merged.pos, data=geno1A,t="n",log='y',yaxt='n',xlab='Genome position', ylab="Mutation frequency",
    ylim=c(10^-4,10^0),xlim=c(340,8500))
eaxis(side = 2, at = 10^((-0):(-(4))), cex=2)
for(k in 1:4){abline(h = 1:10 * 10^(-k), col = "gray80")}
legend("topright",legend=c("1A","1B","3A"), col=cols2[c(1,2,3)], pch=16, bty = "n")
roll1001a<-rollmean(geno1A$mean, k=100, na.rm=T)
roll1001b<-rollmean(geno1B$mean, k=100, na.rm=T)
roll1003a<-rollmean(geno3A$mean, k=100, na.rm=T)

roll1001a<-c(rep(NA, times=100),roll1001a)
roll1001b<-c(rep(NA, times=100),roll1001b)
roll1003a<-c(rep(NA, times=100),roll1003a)

lines(roll1001a,  col=cols2[1],lwd=2)
lines(roll1001b,  col=cols2[2],lwd=2)
lines(roll1003a,  col=cols2[3],lwd=2)
dev.off()


##############
# subsamples mut freq from 1B & 3A, and plot across the genome by Genes
colors2<-paste0(cols2,"66")
colors10<-qualitative_hcl(10, palette="Dark3")
genenames<-genes$Gene[1:12]

#random resampling (1B and 3A)

for (g in 2:3){
        dt<-read.csv(paste0("Output_all/Ts_summary_",geno[g],".csv"),stringsAsFactors = F)
        dt<-dt[,-1]
        s<-ncol(dt)-2
        
        if (g==2) n<-12
        if (g==3) n<-20
        
        subset<-list()
        for (i in 1:10){
                sub<-sample.int(s, n, replace = F)
                subset.dt<-dt[,c(1,(sub+1))]
                subset.dt$mean<-rowMeans(subset.dt[,2:(n+1)],na.rm=T)
                subset[[i]]<-subset.dt
        }
        
        pdf(paste0("Output_all/",geno[g],"_subset.MF.byGene.pdf"),width=12,height=18)
        par(mfrow = c(6,2))
        for (j in 1:length(genenames)){
                plotname<-genenames[j]
                data<-get(paste0("geno",geno[g]))
                plot(data$mean[data$gene==plotname], pch=16, cex=.7,col=colors2[g], main=plotname, ylab="Mutation frequency", xlab="Genome position")

                roll503A<-list()
                r50<-data.frame(merged.pos=c(1:9410))
                for ( i in 1:10){
                        df<-subset[[i]]
                        df<-df[which(data$gene==plotname),]
                        roll50<-rollmean(df$mean, k=50, na.rm=T)
                        roll50<-c(rep(NA, times=49),roll50)
                        roll503A[[i]]<-roll50
                        lines(roll50, col=colors10[i], lwd=1)
                }
        }
        dev.off()
}

# Plot resampling resutls of all genotypes to better visualize the trend
pdf(paste0("Output_all/All_Subset.MF.50.byGene.pdf"),width=12,height=18)
par(mfrow = c(6,2))
for (i in 1:length(genenames)){
        plotname<-genenames[i]
        plot(geno1A$mean[geno1A$gene==plotname], pch=16, cex=.7,col=colors2[1], main=plotname, ylab="Mutation frequency", xlab="Genome position")
        legend("topright",legend=c("1A","1B","3A"), col=colors[c(1,2,3)], pch=16, bty = "n")
        points(geno1B$mean[geno1B$gene==plotname], pch=16, cex=.7,col=colors2[2])
        points(geno3A$mean[geno3A$gene==plotname], pch=16, cex=.7,col=colors2[3])
        
        for ( i in 1:10){
                df<-subset1B[[i]]
                df<-df[which(geno1B$gene==plotname),]
                roll50<-rollmean(df$mean, k=50, na.rm=T)
                roll50<-c(rep(NA, times=49),roll50)
                lines(roll50, col=cols2[2], lwd=1)
        }
        
        for ( i in 1:10){
                df<-subset3A[[i]]
                df<-df[which(geno3A$gene==plotname),]
                roll50<-rollmean(df$mean, k=50, na.rm=T)
                roll50<-c(rep(NA, times=49),roll50)
                lines(roll50, col=cols2[3], lwd=1)
        }
        
        roll501a<-rollmean(geno1A$mean[geno1A$gene==plotname], k=50, na.rm=T)
        roll501a<-c(rep(NA, times=49),roll501a)
        lines(roll501a,  col=cols2[1],lwd=1)
}
dev.off()



#####################  
# plot the mut freq means & se

summaryT<-read.csv("Output_all/MutFreqSummary_NonFiltered.Table.csv",stringsAsFactors = F)

p1<-ggplot(summaryT,aes(x=X,y=Mean,group=Genotype, color=Genotype))+
        geom_point(position=position_dodge(width=0.3))+scale_color_manual(values=cols2)+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, position=position_dodge(width=0.3))+
        theme_light()+
        theme(axis.title.x=element_blank())+ylab("Mean mutation frequency")

ggsave(filename="Output_all/Ave.Mutfreq_by.type_by.genotypes.pdf",plot=p1, width = 9, height = 8)


## with SD 
#ggplot(summaryT,aes(x=X,y=Mean,group=Genotype, color=Genotype))+geom_point(position=position_dodge(width=0.3))+scale_color_manual(values=cols2)+
#        geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.1, position=position_dodge(width=0.3))+ theme_linedraw()


## plot mean and SE for each gene
mf2<-mf[!is.na(mf$mean),]
SumMFGenes<-aggregate(mf2$mean,by=list(mf2$genotype,mf2$gene),FUN=mean)
SumMF.se<-aggregate(mf2$mean,by=list(mf2$genotype,mf2$gene),FUN=std.error)

sumG<-cbind(SumMFGenes, SumMF.se$x)
colnames(sumG)<-c("Genotype","Gene","Mean","SE")

colors12<-qualitative_hcl(12, palette="Dark3")
p3<-ggplot(sumG, aes(x=Gene, y=Mean, group=Genotype, color=Genotype))+
        geom_point(position=position_dodge(width=0.3))+scale_color_manual(values=cols2)+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, position=position_dodge(width=0.3))+
        theme_bw()+theme(axis.text=element_text(size=11), axis.title=element_text(size=13))+ylab("Mean mutation frequency")

ggsave(filename="Output_all/Ave.Mutfreq_by.gene_by.genotype.pdf",plot=p3, width = 10, height = 8)
#ggsave(filename="Output_all/Ave.Mutfreq_by.gene_by.genotype2.pdf",plot=p3, width = 8.5, height = 5)

