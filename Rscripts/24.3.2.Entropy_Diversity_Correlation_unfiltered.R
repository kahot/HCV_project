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

#calculate the values using unfiltered data:

Summary<-read.csv("Output_all/Unfiltered/Ts.MinorVariant_Mean_3genotypes.csv",stringsAsFactors = F,row.names = 1)

for (g in 1:3){
       entropy<-read.table(paste0("Data/HCV",geno[g],"_logo_data.txt"))
       colnames(entropy)<-c("pos","A","C","G","T","Entropy","Low","High","Weight")
       
       mutfreq<-Summary[,c(paste0("org.pos.",geno[g]),paste0("mean.",geno[g]))]
       colnames(mutfreq)<-c("pos","mean")
       mutfreq<-mutfreq[!is.na(mutfreq$pos),]
       
       beginpos<-mutfreq$pos[!is.na(mutfreq$mean)][1]
       entropy<-entropy[beginpos:nrow(entropy),]
       #colnames(entropy)<-c("pos","A","C","G","T","Entropy","Low","High","Weight")
       comparison<-merge(entropy,mutfreq, by="pos")
       results<-cor.test(comparison$Entropy,comparison$mean, method = "spearman")
       print(results)
       pdf(paste0("Output_all/Diversity/Correlation_SeqDiveristy-MutFreq_inverse.",geno[g],".pdf"))
       plot(comparison$mean,-(comparison$Entropy), ylab="", xlab="",
            pch=16,col=cols2[g],cex=0.5, cex.axis=1.2)
       mtext("Intrahost diversity (mean mutation frequency)",1,2.2, cex=1.2)
       mtext("Interhost diversity (-Shannon's entropy)",2,2.2, cex=1.2)
       abline(lm(-(comparison$Entropy)~comparison$mean), col = "gray70")
       rho<-as.numeric(results[[4]])*-1
       rho<-format(round(rho,3), nsmall=3)
       if (results[[3]]>=0.05) star<-""
       if (results[[3]]<0.05&results[[3]]>=0.01) star<-"*"
       if (results[[3]]<0.01&results[[3]]>=0.001) star<-"**"
       if (results[[3]]<0.001) star<-"***"
       
       text(x=max(comparison$mean, na.rm=T)-0.015,y=-0.2, labels=paste0("rho = ",rho,star),cex=1.1)
       dev.off()
}



## nucleotide diversity correlation between populations vs. within population

# Use genbank 1B dataset for 1B entropy calculation (not enough samples)
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
# compare pairwise
#1. Entropy

#adjust for indels 
Positions<-read.csv("Data/MergedPositionInfo.csv", row.names = 1)

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
        #write.csv(comparison, paste0("Output_all/Diversity/Entropy_",g1,"-",g2,".csv"))

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


#2. Mutartion frequency
Summary<-read.csv("Output_all/Unfiltered/Ts.MinorVariant_Mean_3genotypes.csv",stringsAsFactors = F,row.names = 1)


#1. unfiltered Mut. freq (including the sites Ref != Maj)
# compare correlation of mutation frequency at each site
for (g in 1:3){
        g1<-comb[g,1]
        g2<-comb[g,2]
        
        #beginpos<-mutfreq$pos[1]
        
        mutf<-Summary[,c("merged.pos",paste0("mean.",g1),paste0("mean.",g2))]
        #mutf<-mutf[260:8600,]
        results<-cor.test(mutf[,2],mutf[,3], method = "spearman")
        print(results)
        pdf(paste0("Output_all/Diversity/Unfiltered_Cor_MutFreq.",g1,"-",g2,".pdf"), width=6, height=6)
        plot(mutf[,2],mutf[,3], ylab="", xlab="",
             pch=16,col=cols2[g],cex=0.6, cex.axis=1.2, xlim=c(0,max(mutf[,2], na.rm=T)+0.01),ylim=c(0,max(mutf[,3], na.rm=T)+0.01))
        mtext(paste0("Mutation frequency (",g1,")"),1,2.2, cex=1.2)
        mtext(paste0("Mutation frequency (",g2,")"),2,2.2, cex=1.2)
        abline(lm(mutf[,2]~mutf[,3]), col = "gray70")
        rho<-as.numeric(results[[4]])
        rho<-format(round(rho,3), nsmall=3)
        if (results[[3]]>=0.05) star<-""
        if (results[[3]]<0.05&results[[3]]>=0.01) star<-"*"
        if (results[[3]]<0.01&results[[3]]>=0.001) star<-"**"
        if (results[[3]]<0.001) star<-"***"
        
        text(x=max(mutf[,2], na.rm=T)-0.01,y=max(mutf[,3],na.rm=T), labels=paste0("rho = ",rho,star),cex=1.1)
        dev.off()
        
        #Plot with axis limits to 0.1
        pdf(paste0("Output_all/Diversity/Unfiltered_Cor_MutFreq.",g1,"-",g2,".green.pdf"), width=5.7, height=6)
        plot(mutf[,2],mutf[,3], ylab="", xlab="",
             pch=16,col=cols2[3],cex=0.6, cex.axis=1.2, xlim=c(0,0.1),ylim=c(0,0.1))
        mtext(paste0("Mutation frequency (",g1,")"),1,2.2, cex=1.2)
        mtext(paste0("Mutation frequency (",g2,")"),2,2.2, cex=1.2)
        abline(lm(mutf[,2]~mutf[,3]), col = "gray70")
        rho<-as.numeric(results[[4]])
        rho<-format(round(rho,3), nsmall=3)
        if (results[[3]]>=0.05) star<-""
        if (results[[3]]<0.05&results[[3]]>=0.01) star<-"*"
        if (results[[3]]<0.01&results[[3]]>=0.001) star<-"**"
        if (results[[3]]<0.001) star<-"***"
        
        text(x=0.088,y=max(0.1,na.rm=T), labels=paste0("rho = ",rho,star),cex=1.1)
        dev.off()
}

#2. unfiltered Mut. freq (excluding the sites Ref != Maj)
Summary2<-read.csv("Output_all/Unfiltered/Ts.Same_Mean_3genotypes.csv",stringsAsFactors = F,row.names = 1)
for (g in 1:3){
        g1<-comb[g,1]
        g2<-comb[g,2]
        
        mutf<-Summary2[,c("merged.pos",paste0("mean.",g1),paste0("mean.",g2))]
        #mutf<-mutf[260:8600,]
        results<-cor.test(mutf[,2],mutf[,3], method = "spearman")
        print(results)
        pdf(paste0("Output_all/Diversity/Unfiltered_Same_Cor_MutFreq.",g1,"-",g2,".pdf"), width=6, height=6)
        plot(mutf[,2],mutf[,3], ylab="", xlab="",
             pch=16,col=cols2[g],cex=0.6, cex.axis=1.2, xlim=c(0,max(mutf[,2:3], na.rm=T)+0.01),ylim=c(0,max(mutf[,2:3], na.rm=T)+0.01))
        mtext(paste0("Mutation frequency (",g1,")"),1,2.2, cex=1.2)
        mtext(paste0("Mutation frequency (",g2,")"),2,2.2, cex=1.2)
        abline(lm(mutf[,2]~mutf[,3]), col = "gray70")
        rho<-as.numeric(results[[4]])
        rho<-format(round(rho,3), nsmall=3)
        if (results[[3]]>=0.05) star<-""
        if (results[[3]]<0.05&results[[3]]>=0.01) star<-"*"
        if (results[[3]]<0.01&results[[3]]>=0.001) star<-"**"
        if (results[[3]]<0.001) star<-"***"
        
        text(x=max(mutf[,2:3], na.rm=T)-0.008,y=max(mutf[,2:3],na.rm=T), labels=paste0("rho = ",rho,star),cex=1.1)
        dev.off()
}

#2. unfiltered Mut. freq (excluding the sites Ref != Maj, and excluding the sites with different nucleotide between the genotypes)
Summary2<-read.csv("Output_all/Unfiltered/Ts.Same_Mean_3genotypes.csv",stringsAsFactors = F,row.names = 1)
for (g in 1:3){
        g1<-comb[g,1]
        g2<-comb[g,2]
        
        mutf<-Summary2[,c("merged.pos",paste0("mean.",g1),paste0("ref.",g1),paste0("mean.",g2),paste0("ref.",g2))]
        #mutf<-mutf[260:8600,]
        
        mutf2<-mutf[mutf[,3]==mutf[5],]
        mutf2<-mutf2[complete.cases(mutf2),]
        results<-cor.test(mutf2[,2],mutf2[,4], method = "spearman")
        print(results)
        pdf(paste0("Output_all/Diversity/Unfiltered_SameAncestral_Cor_MutFreq.",g1,"-",g2,".pdf"), width=6, height=6)
        plot(mutf2[,2],mutf2[,4], ylab="", xlab="",
             pch=16,col=cols2[g],cex=0.6, cex.axis=1.2, xlim=c(0,max(mutf2[,c(2)], na.rm=T)+0.01),ylim=c(0,max(mutf2[,c(4)], na.rm=T)+0.01))
        mtext(paste0("Mutation frequency (",g1,")"),1,2.2, cex=1.2)
        mtext(paste0("Mutation frequency (",g2,")"),2,2.2, cex=1.2)
        abline(lm(mutf2[,2]~mutf2[,4]), col = "gray70")
        rho<-as.numeric(results[[4]])
        rho<-format(round(rho,3), nsmall=3)
        if (results[[3]]>=0.05) star<-""
        if (results[[3]]<0.05&results[[3]]>=0.01) star<-"*"
        if (results[[3]]<0.01&results[[3]]>=0.001) star<-"**"
        if (results[[3]]<0.001) star<-"***"
        
        if (g==1) x1<-0.025
        if (g==2) x1<-0.02
        if (g==3) x1<-0.015
        
        text(x=max(mutf[,2], na.rm=T)-x1,y=max(mutf[,4],na.rm=T), labels=paste0("rho = ",rho,star),cex=1.1)
        dev.off()
}


###

dfm<-read.csv("Output_all/Unfiltered/Ts.Same_Mean_3genotypes.csv", row.names = 1)
dfm<-dfm[264:8400,]
#dfe<-read.csv("Output_all/Diversity/Entropy_merged_3genotypes.csv")
#dfe<-dfe[1:(8591-341),]

div<-data.frame(Genotype=geno)
div$MutFreq<-colMeans(dfm[c("mean.1A","mean.1B","mean.3A")],na.rm=T)
div$mf.SE<-std.error(dfm[c("mean.1A","mean.1B","mean.3A")],na.rm=T)
#div$Entropy<-colMeans(dfe[5:7],na.rm=T)
#div$Ent.SE<-std.error(dfe[5:7],na.rm=T)


#1. mutation frequency
df<-dfm
colnames(df)[28:30]<-c("MV.Freq.1A","MV.Freq.1B","MV.Freq.3A")
#2. Entropty
EN<-read.csv("Output_all/Diversity/Entropy_merged_3genotypes.csv", stringsAsFactors = F)
EN<- EN[1:8620,]
library(psych)
#pdf("Output_all/Diversity/Cor_MVF.pdf", width = 7, height = 5)
pairs.panels(df[28:30],method="spearman", pch=".", hist.col = "#66CCEECC",,ellipses=F, lm=T,rug=F, col="#66CCEE")
rect(.3,.8,.6,1.1, col="white",border=NA)
text(.45,.95,labels="0.766",cex=2)
text(.56,1, labels="***", cex=1.3)

rect(.7,.8,.9,1.1, col="white",border=NA)
text(.82,.95,labels="0.692",cex=2)
text(.93,1, labels="***", cex=1.3)

rect(.7,.35,.9,.6, col="white",border=NA)
text(.82,.45,labels="0.682",cex=2)
text(.93,.5, labels="***", cex=1.3)

#dev.off()

#pdf("Output_all/Diversity/Cor_entropy.pdf", width = 8, height = 5.7)
pairs.panels(EN[5:7],method="spearman",pch=".", hist.col = "#66CCEECC",,ellipses=F, lm=T,rug=F, col="#66CCEE")
text(.96,1,labels="***",cex=1.5)
text(.6,1,labels="***",cex=1.5)
text(.96,.5,labels="***",cex=1.5)
#dev.off()



#similar plot using a different package
library(PerformanceAnalytics)
chart.Correlation(df[28:30],method="spearman", pch=16, cex=1, cex.cor=.5, histogram=T, col="gray")

library(corrplot)
x <- cor(df[3:5],use="na.or.complete",method = "spearman")
corrplot(x, type="upper", order="hclust")
