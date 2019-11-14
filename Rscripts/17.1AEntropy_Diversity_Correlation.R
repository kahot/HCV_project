library(ggplot2)
library(colorspace)
library(purrr)
#library(ggpubr)
#library(ggthemes)
#library(tidyverse)
#library(zoo)
#library(seqinr)
#library(pegas)
#library(reshape)


cols2<-c("#66CCEE","#EE667799" ,"#22883399")
colors2<-qualitative_hcl(6, palette="Dark3")

# Using 195 1a sequences
entropy<-read.table("Data/HCV1A_logo_data.txt")
colnames(entropy)<-c("pos","A","C","G","T","Entropy","Low","High","Weight")

# 1. With filtered mut freq
mutfreq<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv", stringsAsFactors = F, row.names = 1)
beginpos<-mutfreq$pos[1]
endpos<-mutfreq$pos[nrow(mutfreq)]
entropy1<-entropy[beginpos:nrow(entropy),]

compare<-merge(entropy,mutfreq, by="pos")
results<-cor.test(compare$Entropy,compare$mean, method = "spearman")
print(results)
#rho=-0.6768761 

rho1<-as.numeric(results[[4]])*-1
rho1<-format(round(rho1,3), nsmall=3)
if (results[[3]]>=0.05) star<-""
if (results[[3]]<0.05&results[[3]]>=0.01) star<-"*"
if (results[[3]]<0.01&results[[3]]>=0.001) star<-"**"
if (results[[3]]<0.001) star<-"***"


pdf("Output1A/SummaryFig.Filtered/Entropy-MF-Filtered_195.pdf")
plot(compare$mean,-(compare$Entropy), ylab="", xlab="",
     pch=16,col=colors2[5],cex=0.6, xaxt="n")
mtext(expression(paste("", italic("In vivo"), "diversity (average mutation frequency)")),1, 2.5, cex=1.2)
mtext("Among-host diversity (-Shannon's entropy)",2,2.2, cex=1.2)
axis(1, at=c(0,0.01,0.02,0.03), labels=c(0,0.01,0.02,0.03))

abline(lm(-(compare$Entropy)~compare$mean), col = "gray70")
text(x=max(compare$mean)-0.005,y=-0.2, labels=expression(paste(rho, " = 0.677***")),cex=1.1)
dev.off()

##GGplot 
#rh <- expression(""~rho == 0.677)
rh1 <- expression(""~rho == 0.677~"***")
ggplot(compare, aes(x=mean, y=Entropy*-1))+
        geom_point(color=colors2[5], size=0.7)+
        theme(axis.title.x=element_blank())+ylab("Among-host diversity (-Shannon's entropy)")+
        xlab(expression(paste("", italic("In vivo "), "diversity (average mutation frequency)")))+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
        geom_text(x=0.03, y=-0.2,label=rh1,  parse = T)
ggsave("Output1A/SummaryFig.Filtered/Entropy-MF.195.ggplot3.pdf", width = 4,height = 4)

        
## Compare unfiltered mut freq
HCVFiles_overview2<-list.files("Output1A/Overview2/",pattern="overview2.csv")
Overview_sum_ref<-list()
for (i in 1:length(HCVFiles_overview2)){ 
        overviews<-read.csv(paste0("Output1A/Overview2/",HCVFiles_overview2[i]),stringsAsFactors=FALSE)
        overviews<-overviews[,-1]
        Overview_sum_ref[[i]]<-overviews
        names(Overview_sum_ref)[i]<-substr(paste(HCVFiles_overview2[i]),start=1,stop=7)
}

Freq<-data.frame(pos=Overview_sum_ref[[1]][1])
for (i in 1:length(Overview_sum_ref)){
        dat<-Overview_sum_ref[[i]]
        dat<-dat[dat$ref==dat$MajNt,]
        fname<-names(Overview_sum_ref)[i]
        low_reads<-which(dat$TotalReads<1000) 
        dat[low_reads,c(8:19)]<-NA
        Freq<-merge(Freq, dat[,c("pos","freq.Ts.ref")],by="pos", all.x=T)
        colnames(Freq)[(i+1)]<-fname
}

write.csv(Freq, "Output1A/MutFreq/Ts.summary.ref.same.unfilter.csv")
Freq$mean<-rowMeans(Freq[,2:196], na.rm=T)
#remove the last 180 rows
Freq<-Freq[1:8357,]
beginpos<-Freq$pos[1]
endpos<-Freq$pos[nrow(Freq)]

entropy2<-entropy[beginpos:nrow(entropy),]

compare<-merge(entropy2,Freq, by="pos")
results<-cor.test(compare$Entropy,compare$mean, method = "spearman")
print(results)
#rho = -0.7432347 
# p-value < 2.2e-16
rho2<-as.numeric(results[[4]])*-1
rho2<-format(round(rho2,3), nsmall=3)
if (results[[3]]>=0.05) star<-""
if (results[[3]]<0.05&results[[3]]>=0.01) star<-"*"
if (results[[3]]<0.01&results[[3]]>=0.001) star<-"**"
if (results[[3]]<0.001) star<-"***"

#Uniltered mut freq vs. Entropy
pdf("Output1A/SummaryFig.Filtered/Entropy-MF-Unfiltered_195.pdf")
plot(compare$mean,-(compare$Entropy), ylab="", xlab="", pch=16,col=colors2[5],cex=0.6)
mtext(expression(paste("", italic("In vivo"), "diversity (average mutation frequency)")),1, 2.5, cex=1.2)
mtext("Among-host diversity (-Shannon's entropy)",2,2.2, cex=1.2)
abline(lm(-(compare$Entropy)~compare$mean), col = "gray70")
text(x=max(compare$mean, na.rm=T),y=-0.2, labels=expression(paste(rho, " = 0.743***")),cex=1.1,adj=1)
dev.off()



#unfiltered minor variant freq vs. Entropy

mvf<-read.csv("Output_all/Unfiltered/Ts.MinorVariant_Mean_metadata.1A.csv", row.names = 1, stringsAsFactors = F)
mvf<-mvf[,c("org.pos.1A","mean.1A")]
colnames(mvf)<-c("pos","mean")
comparison1<-merge(entropy,mvf, by="pos")
results<-cor.test(comparison1$Entropy,comparison1$mean, method = "spearman")
print(results)
#S = 1.6616e+11, p-value < 2.2e-16
#        rho 
#-0.7858632 


pdf(paste0("Output1A/SummaryFig.Filtered/Entropy-MVF.1A.pdf"))
plot(comparison1$mean,-(comparison1$Entropy), ylab="", xlab="",
     pch=16,col=colors2[5],cex=0.6, cex.axis=1.2)
mtext("In vivo diversity (mean mutation frequency)",1,2.2, cex=1.2)
mtext("Among host diversity (-Shannon's entropy)",2,2.2, cex=1.2)
abline(lm(-(comparison1$Entropy)~comparison1$mean), col = "gray70")
rho<-as.numeric(results[[4]])
rho<-format(round(rho,3), nsmall=3)
if (results[[3]]>=0.05) star<-""
if (results[[3]]<0.05&results[[3]]>=0.01) star<-"*"
if (results[[3]]<0.01&results[[3]]>=0.001) star<-"**"
if (results[[3]]<0.001) star<-"***"
paste0(rho,star)
text(x=max(comparison1$mean, na.rm=T)-0.015,y=-1.2, labels=expression(paste(rho, " = 0.786***")),cex=1.1)
dev.off()







####
#NCBI HCV1A sequences + 195 consensus (423 +195 = 618 sequences, coding sequnces only)
g=1
entropy<-read.table(paste0("Data/1A_combind_Logo.txt"))
colnames(entropy)<-c("pos","A","C","G","T","Entropy","Low","High","Weight")
entropy$pos<-342:(nrow(entropy)+341)
mutfreq<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv", stringsAsFactors = F, row.names = 1)
compare<-merge(entropy,mutfreq, by="pos")
results<-cor.test(compare$Entropy,compare$mean, method = "spearman")
print(results)
#      rho 
#-0.677469 
pdf(paste0("Output1A/SummaryFig.Filtered/Entropy-MF_618.pdf"))
plot(comparison1$mean,-(comparison1$Entropy), ylab="", xlab="",
     pch=16,col=colors2[5],cex=0.6, cex.axis=1.2)
mtext("In vivo diversity (mean mutation frequency)",1,2.2, cex=1.2)
mtext("Among host diversity (-Shannon's entropy)",2,2.2, cex=1.2)
abline(lm(-(comparison1$Entropy)~comparison1$mean), col = "gray70")
rho<-as.numeric(results[[4]])
rho<-format(round(rho,3), nsmall=3)
if (results[[3]]>=0.05) star<-""
if (results[[3]]<0.05&results[[3]]>=0.01) star<-"*"
if (results[[3]]<0.01&results[[3]]>=0.001) star<-"**"
if (results[[3]]<0.001) star<-"***"
paste0(rho,star)
text(x=max(comparison1$mean, na.rm=T)-0.015,y=-1.2, labels=expression(paste(rho, " = 0.677***")),cex=1.1)
dev.off()



comparison2<-merge(entropy,mvf, by="pos")
results<-cor.test(comparison2$Entropy,comparison2$mean, method = "spearman")
print(results)
#rho 
#-0.7875385 

pdf(paste0("Output1A/SummaryFig.Filtered/Entropy-MVF.618.pdf"))
plot(comparison1$mean,-(comparison1$Entropy), ylab="", xlab="",
     pch=16,col=colors2[5],cex=0.6, cex.axis=1.2)
mtext("In vivo diversity (mean mutation frequency)",1,2.2, cex=1.2)
mtext("Among host diversity (-Shannon's entropy)",2,2.2, cex=1.2)
abline(lm(-(comparison1$Entropy)~comparison1$mean), col = "gray70")
rho<-as.numeric(results[[4]])
rho<-format(round(rho,3), nsmall=3)
if (results[[3]]>=0.05) star<-""
if (results[[3]]<0.05&results[[3]]>=0.01) star<-"*"
if (results[[3]]<0.01&results[[3]]>=0.001) star<-"**"
if (results[[3]]<0.001) star<-"***"
paste0(rho,star)
text(x=max(comparison1$mean, na.rm=T)-0.015,y=-1.2, labels=expression(paste(rho, " = 0.788***")),cex=1.1)
dev.off()
