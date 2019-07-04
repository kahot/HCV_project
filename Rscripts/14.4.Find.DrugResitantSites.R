library(tidyverse)
source("Rscripts/baseRscript.R")

if (FALSE){
HCVFiles_overview3<-list.files("Output1A/Overview3/",pattern="overview3.csv")
FilteredOverview<-list()
for (i in 1:length(HCVFiles_overview3)){ 
        overviews<-read.csv(paste0("Output1A/Overview3/",HCVFiles_overview3[i]),stringsAsFactors=FALSE)
        overviews<-overviews[,-1]
        FilteredOverview[[i]]<-overviews
        names(FilteredOverview)[i]<-substr(paste(HCVFiles_overview3[i]),start=1,stop=7)
        }
}
######



########################################

#specific nucleotide positions with known drug resistant mutations 
dr2<-read.csv("Data/HCV_drugresistance_NTpos.csv")
diff<-list()
diff.count<-list()
DR_mutfreq<-list()
Same_DRmutfreq<-list() 
for (i in 1: length(FilteredOverview)){
        df<-FilteredOverview[[i]]
        dname<-names(FilteredOverview)[i]
        #select specific DR sites
        DRsites<-df[df$pos %in% dr2$NTpos,]
        
        drpos<-data.frame(dr2$NTpos)
        colnames(drpos)<-'pos'
        # should be 38 sites (1 is filtered out so add that site)
        DRsites<-merge(DRsites,drpos, by='pos', all = TRUE)
        
        
        #identify the sites with mutations from the Ref gene
        DRsites$obs<-0
        for (k in 1:nrow(DRsites)){
                if (is.na(DRsites$MajNt[k])|is.na(DRsites$ref[k])) next
                if (DRsites$MajNt[k]!=DRsites$ref[k]){
                        if (dr2$Type[k]=="Ts") {mutnt<-transition(DRsites$ref[k]) 
                                if (DRsites$MajNt[k]==mutnt) DRsites$obs[k]<-1} 
                        if (dr2$Type[k]=='Tv1') {mutnt<-transv1(DRsites$ref[k])
                                if (DRsites$MajNt[k]==mutnt) DRsites$obs[k]<-1}
                        if (dr2$Type[k]=='Tv2') {mutnt<-transv2(DRsites$ref[k])
                                if (DRsites$MajNt[k]==mutnt) DRsites$obs[k]<-1} 
                        
                        if (dr2$Type[k]=='Tv') {mutnt1<-transv1(DRsites$ref[k])
                                mutnt2<-transv2(DRsites$ref[k])
                                if (DRsites$MajNt[k]==mutnt1|DRsites$MajNt[k]==mutnt2) DRsites$obs[k]<-1  }    
                        if (dr2$Type[k]=='Ts,Tv1') {mutnt1<-transition(DRsites$ref[k])
                                mutnt2<-transv1(DRsites$ref[k])    
                                if (DRsites$MajNt[k]==mutnt1|DRsites$MajNt[k]==mutnt2) DRsites$obs[k]<-1 }  
                        if (dr2$Type[k]=='Ts,Tv2') {mutnt1<-transition(DRsites$ref[k])
                                mutnt2<-transv2(DRsites$ref[k])    
                                if (DRsites$MajNt[k]==mutnt1|DRsites$MajNt[k]==mutnt2) DRsites$obs[k]<-1 }     
                        if (dr2$Type[k]=='Ts,Tv') DRsites$obs[k]<-1
                }
        }

 
        diff[[i]]<-DRsites[,c("pos","obs")]
        names(diff)[i]<-dname
        diff.count[i]<-sum(DRsites$obs==1)
        names(diff.count)[i]<-dname
        
        DRsites_same<-DRsites[which(DRsites$MajNt==DRsites$ref),]
        write.csv(DRsites,paste0("Output/DrugRes_NTsites_",dname,".csv"))
        
        
        for (k in 1:nrow(DRsites_same)){
                if (dr2$Type[k]=='Ts') DRsites_same$freq[k]<-DRsites_same$freq.Ts[k]
                if (dr2$Type[k]=='Tv1') DRsites_same$freq[k]<-DRsites_same$freq.transv1[k]
                if (dr2$Type[k]=='Tv2') DRsites_same$freq[k]<-DRsites_same$freq.transv2[k]
                if (dr2$Type[k]=='Tv') DRsites_same$freq[k]<-DRsites_same$freq.transv[k]
                if (dr2$Type[k]=='Ts,Tv1') DRsites_same$freq[k]<-(DRsites_same$freq.Ts[k]+DRsites_same$freq.transv1[k])
                if (dr2$Type[k]=='Ts,Tv2') DRsites_same$freq[k]<-(DRsites_same$freq.Ts[k]+DRsites_same$freq.transv2[k])
                if (dr2$Type[k]=='Ts,Tv') DRsites_same$freq[k]<-DRsites_same$freq.mutations[k]
        }
     
        Same_DRmutfreq[[i]]<-DRsites_same[,c("pos","freq")]
        names(Same_DRmutfreq)[i]<-dname

        
        for (k in 1:nrow(DRsites)){
                if (dr2$Type[k]=='Ts') DRsites$freq[k]<-DRsites$freq.Ts.ref[k]
                if (dr2$Type[k]=='Tv1') DRsites$freq[k]<-DRsites$freq.transv1.ref[k]
                if (dr2$Type[k]=='Tv2') DRsites$freq[k]<-DRsites$freq.transv2.ref[k]
                if (dr2$Type[k]=='Tv') DRsites$freq[k]<-DRsites$freq.transv.ref[k]
                if (dr2$Type[k]=='Ts,Tv1') DRsites$freq[k]<-(DRsites$freq.Ts.ref[k]+DRsites$freq.transv1.ref[k])
                if (dr2$Type[k]=='Ts,Tv2') DRsites$freq[k]<-(DRsites$freq.Ts.ref[k]+DRsites$freq.transv2.ref[k])
                if (dr2$Type[k]=='Ts,Tv') DRsites$freq[k]<-DRsites$freq.mutations.ref[k]
        }
        
        DR_mutfreq[[i]]<-DRsites[,c("pos","freq")]
        names(DR_mutfreq)[i]<-dname
       
}

for (i in 1:length(DR_mutfreq)) {
        colnames(DR_mutfreq[[i]])<-c("pos",paste0(names(DR_mutfreq[i])))
        colnames(diff[[i]])<-c("pos",paste0(names(diff[i])))
}

DR.mutated.counts<-do.call(rbind,diff.count)

DR_MutFreq<-DR_mutfreq%>% reduce(full_join, by='pos')
DR_diff<-diff %>% reduce(full_join, by='pos')
DR_same<-Same_DRmutfreq%>% reduce(full_join, by='pos')
write.csv(DR_MutFreq, "Output/DrugRes/DrugResistantMutationFreq_summary_filtered.csv")
write.csv(DR_same, "Output/DrugRes/DrugResistant_sameMutFreq_summary_filtered.csv")
write.csv(DR_diff, "Output/DrugRes/DrugResistant_diffMutFreq_summary_filtered.csv")


##calculate the average mutation frequecy of minor alleles at known drug resisitant sites
s<-length(FilteredOverview)
DR_MutFreq$mean<-apply(DR_MutFreq[2:(s+1)],1,mean, na.rm=T )
#plot(DR_MutFreq$mean)

#count the number of patients fixed with alternate NT

DR_diff$total<-apply(DR_diff[2:(s+1)],1,sum, na.rm=T)
DR_diff$Percent<-DR_diff$total/s*100
DR_diff$Per2<-format(round(DR_diff$Percent, 1), nsmall = 1)

#create a figure

colors=c("#44AA99","#0077BB","#CC6677" )
pdf("Output/DrugRes/Drug.resi.allele.freq_Q35.pdf", height = 8, width = 15)
ymin <- -5
plot(0, type = "n", xlim = c(1, nrow(DR_MutFreq)), ylim = c(ymin, 0), axes = FALSE, ylab = "Frequency of drug resistant allele(s)", 
     xlab = "")
axis(side = 2, at = seq(0, ymin, by=-1), labels = expression(10^0, 10^-1, 10^-2, 10^-3, 10^-4,10^-5), las = 2)

n<-seq(1, by = 2, len = (nrow(DR_MutFreq)/2))

ns3n<-nrow(dr2[dr2$Gene=="NS3",])
ns5an<-nrow(dr2[dr2$Gene=="NS5A",])

genes<-data.frame("name"= c("NS3","NS5A","NS5B"), "Begin"= c(1,ns3n+1,(ns3n+ns5an+1)),"End"= c(ns3n,(ns3n+ns5an),nrow(dr2)))
for (i in 1:3){
        nvec<-n[((genes$Begin[i]+1)/2):(genes$End[i]/2)]
        abline(v=nvec, col=paste0(colors[i],"66"),lty=1, lwd=30)
        abline(v=nvec+1, col=paste0(colors[i],"1A"),lty=1, lwd=32)
}

for (i in 1:nrow(DR_MutFreq)){
        if (i<=genes$End[1]){
                xjit <- rnorm(s, 0, .1)
                points(rep(i,s)+xjit, log10(DR_MutFreq[i,2:(s+1)]),pch='.',col=colors[1],cex = 1.5)}
        if (i<=genes$End[2]& i>=genes$Begin[2]){
                xjit <- rnorm(s, 0, .1)
                points(rep(i,s)+xjit, log10(DR_MutFreq[i,2:(s+1)]),pch='.',col=colors[2])}
        if (i<=genes$End[3]& i>=genes$Begin[3]){
                xjit <- rnorm(s, 0, .1)
                points(rep(i,s)+xjit, log10(DR_MutFreq[i,2:(s+1)]),pch='.',col=colors[3])}
}

#mtext(side = 3, at = 1:nrow(DR_MutFreq), text = paste(DR_diff$total), cex = .8)

for (i in 1:nrow(dr2)){
        if (dr2$Type[i]=='Ts') mtext(side = 3, at = i , text = paste(DR_diff$total[i]),col="#117733",cex = .8)
        if (dr2$Type[i]=='Tv1'|dr2$Type[i]=='Tv2'|dr2$Type[i]=='Tv') mtext(side = 3, at = i, text = paste(DR_diff$total[i]),col="#EE7733",cex = .8)
        if (dr2$Type[i]=='Ts,Tv2'|dr2$Type[i]=='Ts,Tv') mtext(side = 3, at = i, text = paste(DR_diff$total[i]),col="#0077BB",cex = .8)
}
                                                                          
mtext(side = 1, at = 1:nrow(DR_MutFreq), text = paste(dr2$Name), las=2,padj=0, cex = .8)

abline(v=18+.5,col='gray50', lwd=3)
abline(v=34.5,col='gray50', lwd=3 )

rect(0.5,-5,18.5,-5.19,density = NULL, angle = 45,col="white",border =colors[1])
text(9,-5.08,paste0(genes$name[1]),col="black", cex=.8)
rect(18.5,-5,34.5,-5.19,density = NULL, angle = 45,col="white",border =colors[2])
text(26.5,-5.08,paste0(genes$name[2]),col="black", cex=.8)
rect(34.5,-5,38.5,-5.19,density = NULL, angle = 45,col="white",border =colors[3])
text(36.5,-5.08,paste0(genes$name[3]),col="black", cex=.8)

dev.off()

