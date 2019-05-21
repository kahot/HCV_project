library(ggplot2)
library(reshape)
library(ggpubr)
library(ggthemes)
library(tidyverse)
library(zoo)
library(seqinr)
library(pegas)

cols2<-c("#66CCEE","#EE667799" ,"#22883399")
geno<-c("1A","1B","3A")
for (g in 1:3){
       entropy<-read.table(paste0("Data/HCV",geno[g],"_logo_data.txt"))
       colnames(entoropy)<-c("pos","A","C","G","T","Entropy","Low","High","Weight")
       
       mutfreq<-read.csv(paste0("Output",geno[g],"/MutFreq/Ave.MF.Total_Mutations_",geno[g],".csv"), stringsAsFactors = F)
       mutfreq<-mutfreq[,-1]
       beginpos<-mutfreq$pos[1]
       entropy<-entropy[beginpos:nrow(entropy),]
       colnames(entropy)<-c("pos","A","C","G","T","Entropy","Low","High","Weight")
       comparison<-merge(entropy,mutfreq, by="pos")
       results<-cor.test(comparison$Entropy,comparison$mean, method = "spearman")
       print(results)
       pdf(paste0("Output_all/Diversity/Correlation_SeqDiveristy-MutFreq_inverse.",geno[g],".pdf"))
       plot(comparison$mean,-(comparison$Entropy), ylab="", xlab="",
            pch=16,col=cols2[g],cex=0.6, cex.axis=1.2)
       mtext("Intrahost diversity (mean mutation frequency)",1,2.2, cex=1.2)
       mtext("Interhost diversity (-Shannon's entropy)",2,2.2, cex=1.2)
       abline(lm(-(comparison$Entropy)~comparison$mean), col = "gray70")
       rho<-as.numeric(results[[4]])
       rho<-format(round(rho,3), nsmall=3)
       if (results[[3]]>=0.05) star<-""
       if (results[[3]]<0.05&results[[3]]>=0.01) star<-"*"
       if (results[[3]]<0.01&results[[3]]>=0.001) star<-"**"
       if (results[[3]]<0.001) star<-"***"
       
       max(comparison$mean)
       text(x=max(comparison$mean)-0.015,y=-0.2, labels=paste0("rho = ",rho,star),cex=1.1)
       dev.off()
}



## nucleotide diversity correlation between populations vs. within population

# Use genbank 1B dataset
g=2
entropy<-read.table(paste0("Data/HCV1B_ncbi_logo_data.txt"))
colnames(entropy)<-c("row","A","C","G","T","Entropy","Low","High","Weight")
entropy$pos<-c(342:(nrow(entropy)+341))
entropy1B<-entropy
write.csv(entropy1B, "Data/HCV1B_logo_data.txt")

mutfreq<-read.csv(paste0("Output",geno[g],"/MutFreq/Ave.MF.Total_Mutations_",geno[g],".csv"), stringsAsFactors = F)
mutfreq<-mutfreq[342:nrow(mutfreq),-1]

comparison<-merge(entropy,mutfreq, by="pos")

results<-cor.test(comparison$Entropy,comparison$mean, method = "spearman")
print(results)

pdf(paste0("Output_all/Diversity/Cor_SeqDiveristy(NCBI)-MutFreq_inverse.",geno[g],".pdf"))
plot(comparison$mean,-(comparison$Entropy), ylab="", xlab="",
     pch=16,col=cols2[g],cex=0.6, cex.axis=1.2)
mtext("Intrahost diversity (mean mutation frequency)",1,2.2, cex=1.2)
mtext("Interhost diversity (-Shannon's entropy)",2,2.2, cex=1.2)
abline(lm(-(comparison$Entropy)~comparison$mean), col = "gray70")
rho<-as.numeric(results[[4]])
rho<-format(round(rho,3), nsmall=3)
if (results[[3]]>=0.05) star<-""
if (results[[3]]<0.05&results[[3]]>=0.01) star<-"*"
if (results[[3]]<0.01&results[[3]]>=0.001) star<-"**"
if (results[[3]]<0.001) star<-"***"

max(comparison$mean)
text(x=max(comparison$mean)-0.015,y=-0.2, labels=paste0("rho = ",rho,star),cex=1.1)

dev.off()


##############
# compare 1A vs. 3A

#need to adjust for indels 

Positions<-read.csv("Data/MergedPositionInfo.csv")
Positions<-Positions[,-1]

pos<-data.frame(merged.pos=Positions$merged.pos)
geno<-c("1A","1B","3A")
comb<-t(combn(geno,2))
for (g in 1:3){
        g1<-comb[g,1]
        g2<-comb[g,2]
        
        entropy1<-read.table(paste0("Data/HCV",g1,"_logo_data.txt"))
        entropy2<-read.table(paste0("Data/HCV",g2,"_logo_data.txt"))

        colnames(entropy1)<-c("pos","A","C","G","T","Entropy1","Low","High","Weight")
        colnames(entropy2)<-c("pos","A","C","G","T","Entropy2","Low","High","Weight")
        
        orgpos1<-paste0("org.pos.",g1)
        orgpos2<-paste0("org.pos.",g2)
        
        ent1<-merge(Positions,entropy1, by.x=orgpos1,  by.y="pos",all.x=T )
        ent2<-merge(Positions,entropy2, by.x=orgpos2,  by.y="pos",all.x=T )
        
        ent12<-merge(pos, ent1, by="merged.pos", all.x=T)
        ent22<-merge(pos, ent2, by="merged.pos", all.x=T)
        
        comparison<-merge(ent12,ent22, by="merged.pos")
        write.csv(comparison, paste0("Output_all/Diversity/Entropy_",g1,"-",g2,".csv"))
}
        comparison<-comparison[342:8665,]

        results<-cor.test(comparison$Entropy1,comparison$Entropy2, method = "spearman")
        print(results)
        pdf(paste0("Output_all/Diversity/Correlation_entropy.",g1,"-",g2,".pdf"))
        plot(comparison$Entropy1,comparison$Entropy2, ylab="", xlab="",
             pch=16,col=cols2[g],cex=0.6, cex.axis=1.2)
        mtext(paste0("Entropy (",g1,")"),1,2.2, cex=1.2)
        mtext(paste0("Entropy (",g2,")"),2,2.2, cex=1.2)
        abline(lm(comparison$Entropy1~comparison$Entropy2), col = "gray70")
        rho<-as.numeric(results[[4]])
        rho<-format(round(rho,3), nsmall=3)
        if (results[[3]]>=0.05) star<-""
        if (results[[3]]<0.05&results[[3]]>=0.01) star<-"*"
        if (results[[3]]<0.01&results[[3]]>=0.001) star<-"**"
        if (results[[3]]<0.001) star<-"***"
        
        text(x=min(comparison$Entropy1, na.rm=T)+0.2,y=min(comparison$Entropy1,na.rm=T)+0.1, labels=paste0("rho = ",rho,star),cex=1.1)
        dev.off()
}


# compare correlation of mutation frequency at each site
for (g in 1:3){
        g1<-comb[g,1]
        g2<-comb[g,2]
        
        beginpos<-mutfreq$pos[1]
        
        mutf1<-read.csv(paste0("Output",g1,"/MutFreq/Ave.MF.Total_Mutations_",g1,".csv"), stringsAsFactors = F)
        mutf2<-read.csv(paste0("Output",g2,"/MutFreq/Ave.MF.Total_Mutations_",g2,".csv"), stringsAsFactors = F)
        mutf1<-mutf1[,-1]
        mutf2<-mutf2[,-1]
        
        
        orgpos1<-paste0("org.pos.",g1)
        orgpos2<-paste0("org.pos.",g2)
        
        mf1<-merge(Positions,mutf1, by.x=orgpos1,  by.y="pos",all.x=T )
        mf2<-merge(Positions,mutf2, by.x=orgpos2,  by.y="pos",all.x=T )
        
        mf12<-merge(pos, mf1, by="merged.pos",all.x=T)
        mf22<-merge(pos, mf2, by="merged.pos",all.x=T)
        
    
        comparison<-merge(mf12,mf22, by="merged.pos")
        comparison<-comparison[260:8600,]
        results<-cor.test(comparison$mean.x,comparison$mean.y, method = "spearman")
        print(results)
        pdf(paste0("Output_all/Diversity/Correlation_MutFreq.",g1,"-",g2,".pdf"))
        plot(comparison$mean.x,comparison$mean.y, ylab="", xlab="",
             pch=16,col=cols2[g],cex=0.6, cex.axis=1.2)
        mtext(paste0("Mutation frequency (",g1,")"),1,2.2, cex=1.2)
        mtext(paste0("Mutation frequency (",g2,")"),2,2.2, cex=1.2)
        abline(lm(comparison$mean.x~comparison$mean.y), col = "gray70")
        rho<-as.numeric(results[[4]])
        rho<-format(round(rho,3), nsmall=3)
        if (results[[3]]>=0.05) star<-""
        if (results[[3]]<0.05&results[[3]]>=0.01) star<-"*"
        if (results[[3]]<0.01&results[[3]]>=0.001) star<-"**"
        if (results[[3]]<0.001) star<-"***"
        
        text(x=max(comparison$mean.x, na.rm=T)-0.02,y=max(comparison$mean.y,na.rm=T), labels=paste0("rho = ",rho,star),cex=1.1)
        dev.off()
}


###

dfm<-read.csv("Output_all/Diversity/TsMF_merged_3genotypes.csv")
dfm<-dfm[342:8591,]
dfe<-read.csv("Output_all/Diversity/Entropy_merged_3genotypes.csv")
dfe<-dfe[1:(8591-341),]

div<-data.frame(Genotype=geno)
div$MutFreq<-colMeans(dfm[3:5],na.rm=T)
div$mf.SE<-std.error(dfm[3:5],na.rm=T)
div$Entropy<-colMeans(dfe[5:7],na.rm=T)
div$Ent.SE<-std.error(dfe[5:7],na.rm=T)



#calculate nucleotide diversity

plot(div$MutFreq)
plot(div$Entropy)


####
library(scatterplot3d)
library(colorspace)

df<-read.csv("Output_all/TsMF_merged_3genotypes.csv")
df<-df[286:8591,]

scatterplot3d(df[,3:5],pch = 16, color="#66CCEE")


colors12<-qualitative_hcl(12, palette="Dark3")
cls<-colors12[as.numeric(df$gene)]
scatterplot3d(df[,3:5],pch = 16, color=cls)

#assumption test 
shapiro.test(comparison$mean[1:5000]) 
#W = 0.73271, p-value < 2.2e-16
shapiro.test(comparison$mean[5001:nrow(comparison)]) 
#W = 0.78515, p-value < 2.2e-16
qqnorm(comparison$mean)
shapiro.test(sqrt(comparison$mean[1:5000])) 
qqnorm(sqrt(comparison$mean))

results<-cor.test(comparison$Entropy,comparison$mean, method = "spearman") #-0.8098487
#	Spearman's rank correlation rho

#data:  comparison$Entropy and comparison$mean
#S = 1.727e+11, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#        rho 
#-0.8098487



pdf("Output/SummaryFigures/Correlation_SeqDiveristy-MutFreq.pdf")
plot(comparison$mean,comparison$Entropy, ylab="Sequence conservation", xlab="Average mutation frequency",
     pch=21,col="steelblue",bg="lightblue",cex=0.6)
abline(lm(comparison$Entropy~comparison$mean), col = "gray70")
text(x=0.09,y=1.2, labels="rho=-0.810***")
dev.off()



pdf("Output/SummaryFigures/Correlation_SeqDiveristy-MutFreq_inverse.pdf")
plot(comparison$mean,-(comparison$Entropy), ylab="", xlab="",
     pch=16,col="#66CCEEB3",cex=0.6, cex.axis=1.2)
mtext("Average mutation frequency ",1,2.2, cex=1.5)
mtext("-Entropy",2,2.2, cex=1.5)
abline(lm(-(comparison$Entropy)~comparison$mean), col = "gray70")
text(x=0.085,y=-0.2, labels="rho=0.810***",cex=1.2)
dev.off()



