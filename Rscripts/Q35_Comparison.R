library(purrr)
library(tidyverse)
library(zoo)
library(plotrix)
library(sfsmisc)



HCVFiles_overview2<-list.files("Output/Overview_RefQ35/",pattern="overview2.csv")
Overview_sum_ref<-list()
Fnames<-c("noRF","noRemap","Q30","Q35")
for (i in 1:length(HCVFiles_overview2)){ 
        overviews<-read.csv(paste0("Output/Overview_RefQ35/",HCVFiles_overview2[i]),stringsAsFactors=FALSE)
        overviews<-overviews[,-1]
        Overview_sum_ref[[i]]<-overviews
        names(Overview_sum_ref)[i]<-Fnames[i]
}

################################################

## put NA when MajNT != ref in all frequencies/sel coeff.
# put NA for the sites with total-reads <100 (3/5/2019)
FilteredOverview1<-list()
for (i in 1:length(Overview_sum_ref)){
        dat<-Overview_sum_ref[[i]]
        filename<-names(Overview_sum_ref)[i]
        low_reads<-which(dat$TotalReads<1000) 
        dat[low_reads,c("freq.Ts.ref","freq.transv1.ref","freq.transv2.ref","freq.transv.ref","freq.mutations.ref","EstSelCoeff","EstSelCoeff_trans1","EstSelCoeff_trans2","EstSelCoeff_transv")]<-NA
        
        mf_high<-which(dat$freq.Ts.ref>=0.2|dat$freq.transv.ref>=0.2) # row numbers of high mut freq
        dat[mf_high,c("freq.Ts.ref","freq.transv1.ref","freq.transv2.ref","freq.transv.ref","freq.mutations.ref")]<-NA
        dat[mf_high,c("EstSelCoeff","EstSelCoeff_trans1","EstSelCoeff_trans2","EstSelCoeff_transv")]<-NA
        FilteredOverview1[[i]]<-dat
        names(FilteredOverview1)[i]<-filename
        #write.csv(dat,paste0("Output/Overview_filtered/",filename,"_overview3.csv"))
        
        
}

###################
T_Freq_all<-list()
for (i in 1:length(FilteredOverview1)){
        dat<-FilteredOverview1[[i]]
        filename<-names(FilteredOverview1)[i]
        T_Freq_all[[i]]<-dat[,c("pos","freq.Ts.ref")] 
        names(T_Freq_all)[i]<-filename
}
#assign column names for the list
for (i in 1:length(T_Freq_all)) {
        colnames(T_Freq_all[[i]])<-c("pos",paste0(names(T_Freq_all[i])))
}


# Remove the sites with >50% NA in mutation frequency in FilteredOverview files, 

TMutFreq<-T_Freq_all %>% purrr::reduce(full_join, by='pos') #8537 sites


write.csv(TMutFreq, file="Output/Q35Compare/D75002_TsMut_comparison.csv")





##### 
colMeans(TMutFreq[2:5], na.rm=T)
#       noRF     noRemap         Q30         Q35 
#0.003948235 0.003964038 0.004408516 0.003943450 

#
dat<-Overview_sum_ref[[3]]
muttypes<-dat[,c("pos","ref","Type","makesCpG","bigAAChange")]
M<-merge(TMutFreq,muttypes,by="pos")


M<-M[M$pos>=342, ]

colMeans(M[M$Type=="syn",2:5],na.rm=T) 
#       noRF     noRemap         Q30         Q35 
#0.006939931 0.007167473 0.007073138 0.006927964 

colMeans(M[M$Type=="nonsyn",2:5],na.rm=T) 
#       noRF     noRemap         Q30         Q35 
#0.002797197 0.002725361 0.003399259 0.002795260 


M2<-M[M$ref=="a"|M$ref=="t",]
colMeans(M2[M2$Type=="syn"&M2$makesCpG==0,2:5],na.rm=T) 
#      noRF    noRemap        Q30        Q35
#0.01013036 0.00996491 0.01112177 0.01006227 

colMeans(M2[M2$Type=="syn"&M2$makesCpG==1,2:5],na.rm=T) 
#      noRF   noRemap        Q30         Q35
#0.01163593 0.01211402 0.01067647 0.01163194 

colMeans(M2[M2$Type=="nonsyn"&M2$makesCpG==0,2:5],na.rm=T) 
#       noRF     noRemap         Q30         Q35 
#0.003788160 0.003759840 0.004572930 0.003786158 

colMeans(M2[M2$Type=="nonsyn"&M2$makesCpG==1,2:5],na.rm=T) 
#        noRF     noRemap         Q30         Q35 
#30.003540852 0.003411613 0.004461925 0.003539714 

#########################################################################################
##Plot across the genome


noRF_h<-hist(M$noRF,breaks=20)
q35_hist<-hist(M$Q35,breaks=20)
q30_hist<-hist(M$Q30,breaks=20)
noRemap_h<-hist(M$noRemap,breaks=20)

xtick<-seq(0,0.2, by=0.1)
pdf("Output/Q35Compare/D75002_mutfreq_histograms.pdf", width =10, height=8 )
par(mfrow=c(2,2))
plot(q30_hist$counts, xlab="Mutation frequency", ylab="counts", xaxt="n", log="y",type='h',lwd=10, lend=2, main="Q30")
mtext(xtick, side=1, line=0.8, at=c(0,10,20))
plot(q35_hist$counts, xlab="Mutation frequency", ylab="counts",xaxt="n",log="y",type='h',lwd=10, lend=2,main="Q35")
mtext(xtick, side=1, line=0.8, at=c(0,10,20))
plot(noRF_h$counts, xlab="Mutation frequency", ylab="counts",xaxt="n",log="y",type='h',lwd=10, lend=2,main="Q35_noRF")
mtext(xtick, side=1, line=0.8, at=c(0,10,20))
plot(noRemap_h$counts, xlab="Mutation frequency", ylab="counts",xaxt="n",log="y",type='h',lwd=10, lend=2,main="Q35_noRemap")
mtext(xtick, side=1, line=0.8, at=c(0,10,20))

dev.off()

####
# plot across the genome 1.
pdf(paste0("Output/Q35Compare/Tsmutfreq_across_genome.pdf"),width=18,height=10)
par(mfrow=c(2,2))

for (i in 2:5){
        plot(M[,i]~pos, data=M,t="n",log='y',yaxt='n',ylab="transition mutation frequency",
             main=paste0(colnames(M)[i]),ylim=c(0.0001,0.1),xlim=c(340,8500))
        eaxis(side = 2, at = 10^((-1):(-(4))), cex=2)
        for(i in 2:4){abline(h = 1:10 * 10^(-i), col = "gray80")}
        points(M[,i]~pos, data=M,pch=20,col="#EE667766",cex=0.5)
        
        #add rolling average
        roll100<-rollmean(M[,i], k=100, na.rm=T)
        roll100.2<-c(rep(NA, times=50),roll100,c(rep(NA, times=49)))
        lines(roll100.2~pos,data=M, col="#AA3377")
}

dev.off()



r1<-wilcox.test(TMutFreq$Q30, TMutFreq$Q35, alternative = "greater", paired =F)
#W = 31930000, p-value < 2.2e-16
r1[[3]] #8.1448e-61
      
wilcox.test(TMutFreq$noRF, TMutFreq$Q35, alternative = "greater", paired =F)                 
#W = 27648000, p-value = 0.4754
wilcox.test(TMutFreq$noRemap, TMutFreq$Q35, alternative = "less", paired =F)                 
#W = 26777000, p-value = 0.9995

dat<-Overview_sum_ref[[3]]
muttypes<-dat[,c("pos","ref","Type","Type.tv1","Type.tv2","WTAA","MUTAA","TVS1_AA","TVS2_AA","makesCpG","makesCpG.tv1","makesCpG.tv2","bigAAChange","bigAAChange.tv1","bigAAChange.tv2")]
TMutFreq2<-merge(muttypes,TMutFreq,by="pos")

wilcox.test(TMutFreq2$Q30[TMutFreq2$ref=="a"], TMutFreq2$Q35[TMutFreq2$ref=="a"], alternative = "greater", paired =F)                 
#W = 1302300, p-value < 2.2e-16
wilcox.test(TMutFreq2$Q30[TMutFreq2$ref=="t"], TMutFreq2$Q35[TMutFreq2$ref=="t"], alternative = "greater", paired =F)                 
#W = 1468900, p-value < 2.2e-16
wilcox.test(TMutFreq2$Q30[TMutFreq2$ref=="c"], TMutFreq2$Q35[TMutFreq2$ref=="c"], alternative = "greater", paired =F)                 
#W = 2936700, p-value < 2.2e-16
wilcox.test(TMutFreq2$Q30[TMutFreq2$ref=="g"], TMutFreq2$Q35[TMutFreq2$ref=="g"], alternative = "greater", paired =F)                 
#W = 2777900, p-value < 2.2e-16

###########################
Q35Ts<-read.csv("Output/MutFreqQ35/Filtered_Ts_MutFreqSummary_Q35.csv",stringsAsFactors = F) #first 50 samples
Q30Ts<-read.csv("Output/MutFreq/Filtered_Transition_MutFreqSummary.csv",stringsAsFactors = F)  #all 196 samples
Q35Ts<-Q35Ts[,-1]
Q30Ts<-Q30Ts[,-1]
Q35Ts<-Q35Ts[Q35Ts$pos>=342,]
Q30Ts<-Q30Ts[Q30Ts$pos>=342,]

Q35Ts$mean<-rowMeans(Q35Ts[2:51],na.rm=T)
Q30Ts$mean<-rowMeans(Q30Ts[2:196],na.rm=T)
wilcox.test(Q30Ts$mean, Q35Ts$mean, alternative = "greater", paired =F)
#data:  Q30Ts$mean and Q35Ts$mean
#W = 33362000, p-value < 2.2e-16

Q35Ts<-merge(muttypes,Q35Ts,by="pos")
Q30Ts<-merge(muttypes,Q30Ts,by="pos")

wilcox.test(Q30Ts$mean[Q30Ts$ref=="a"], Q35Ts$mean[Q35Ts$ref=="a"], alternative = "greater", paired =F)                 
# p-value < 2.2e-16
wilcox.test(Q30Ts$mean[Q30Ts$ref=="t"], Q35Ts$mean[Q35Ts$ref=="t"], alternative = "greater", paired =F)                 
#p-value < 2.2e-16
wilcox.test(Q30Ts$mean[Q30Ts$ref=="c"], Q35Ts$mean[Q35Ts$ref=="c"], alternative = "greater", paired =F)                 
#p-value = 1.439e-11
wilcox.test(Q30Ts$mean[Q30Ts$ref=="g"], Q35Ts$mean[Q35Ts$ref=="g"], alternative = "greater", paired =F)                 
#p-value < 2.2e-16


