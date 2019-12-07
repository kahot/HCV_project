library(tidyverse)
source("Rscripts/baseRscript.R")


#specific nucleotide positions with known drug resistant mutations 
drsites<-read.csv("Data/HCV_drugresistance_Geno.csv",stringsAsFactors = F)

geno1A<-read.csv("Output_all/Filtered/Ts_summary_metadata.1A.csv",stringsAsFactors = F, row.names = 1)

#attach the original position info for all genotypes
for (i in 1:nrow(drsites)){
        drsites$pos.1A[i]<-geno1A$org.pos.1A[which(geno1A$merged.pos==drsites$merged.pos[i])]
        drsites$pos.1B[i]<-geno1A$org.pos.1B[which(geno1A$merged.pos==drsites$merged.pos[i])]
        drsites$pos.3A[i]<-geno1A$org.pos.3A[which(geno1A$merged.pos==drsites$merged.pos[i])]
}

colors=c("#44AA99","#0077BB","#CC6677" )
#########
geno<-c("1A","1B","3A")
for (j in 1:length(geno)){
        HCVFiles<-list.files(paste0("Output",geno[j],"/Overview2/"), pattern="overview2.csv")
        
        dr.sites<-list()
        diff<-list()
        diff.count<-list()
        DR_mutfreq<-list()
        
        #first deal with the DRV with one mutation 
        DR<-drsites[drsites$genotype==geno[j],]
        
        #create an id column
        for (i in 1:nrow(DR)){
                if (DR$Need_both[i]=="y") DR$ID[i]<- paste0(DR$Name[i],'.',DR$merged.pos[i])
                else DR$ID[i]<-DR$Name[i]
        }
        
        for (i in 1:length(HCVFiles)){ 
                df<-read.csv(paste0("Output",geno[j],"/Overview2/",HCVFiles[i]),stringsAsFactors=FALSE, row.names = 1)
                dname<-substr(paste(HCVFiles[i]),start=1,stop=7)
                dr<-DR
                cname<-paste0("pos.",geno[j])
                DRsites<-df[df$pos %in% dr[,cname],]

                #count the number of samples fixed with the RAVs
                dr$obs<-0
                for (k in 1:nrow(dr)){
                        pos<-dr[k,cname]
                        if (is.na(DRsites$MajNt[DRsites$pos==pos])|is.na(DRsites$ref[DRsites$pos==pos])) next
                        if (DRsites$MajNt[DRsites$pos==pos]!=DRsites$ref[DRsites$pos==pos]){
                                if (dr$Type[k]=="Ts") {mutnt<-transition(DRsites$ref[DRsites$pos==pos]) 
                                        if (DRsites$MajNt[DRsites$pos==pos]==mutnt) dr$obs[k]<-1} 
                                if (dr$Type[k]=='Tv1') {mutnt<-transv1(DRsites$ref[DRsites$pos==pos])
                                        if (DRsites$MajNt[DRsites$pos==pos]==mutnt) dr$obs[k]<-1}
                                if (dr$Type[k]=='Tv2') {mutnt<-transv2(DRsites$ref[DRsites$pos==pos])
                                        if (DRsites$MajNt[DRsites$pos==pos]==mutnt) dr$obs[k]<-1} 
                                if (dr$Type[k]=='Tv') {mutnt1<-transv1(DRsites$ref[DRsites$pos==pos])
                                        mutnt2<-transv2(DRsites$ref[DRsites$pos==pos])
                                        if (DRsites$MajNt[DRsites$pos==pos]==mutnt1|DRsites$MajNt[DRsites$pos==pos]==mutnt2) dr$obs[k]<-1  }    
                                if (dr$Type[k]=='Ts,Tv1') {mutnt1<-transition(DRsites$ref[DRsites$pos==pos])
                                        mutnt2<-transv1(DRsites$ref[DRsites$pos==pos])    
                                        if (DRsites$MajNt[DRsites$pos==pos]==mutnt1|DRsites$MajNt[DRsites$pos==pos]==mutnt2) dr$obs[k]<-1 }  
                                if (dr$Type[k]=='Ts,Tv2') {mutnt1<-transition(DRsites$ref[DRsites$pos==pos])
                                        mutnt2<-transv2(DRsites$ref[DRsites$pos==pos])    
                                        if (DRsites$MajNt[DRsites$pos==pos]==mutnt1|DRsites$MajNt[DRsites$pos==pos]==mutnt2) dr$obs[k]<-1 }     
                                if (dr$Type[k]=='Ts,Tv') dr$obs[k]<-1
                        }
                }
                
                #store the observation number in a list
                diff.count[i]<-sum(dr$obs==1)
                names(diff.count)[i]<-dname
                
                #obtain the mutation freq.
                dr$freq<-""
                for (k in 1:nrow(dr)){
                        if (dr$Type[k]=='Tv1')    dr$freq[k]<-DRsites$freq.transv1.ref[DRsites$pos==dr[k,cname]]
                        if (dr$Type[k]=='Tv2')    dr$freq[k]<-DRsites$freq.transv2.ref[DRsites$pos==dr[k,cname]]
                        if (dr$Type[k]=='Tv')     dr$freq[k]<-DRsites$freq.transv.ref [DRsites$pos==dr[k,cname]]
                        if (dr$Type[k]=='Ts')     dr$freq[k]<-DRsites$freq.Ts.ref     [DRsites$pos==dr[k,cname]]
                        if (dr$Type[k]=='Ts,Tv1') dr$freq[k]<-(DRsites$freq.Ts.ref    [DRsites$pos==dr[k,cname]]+DRsites$freq.transv1.ref[DRsites$pos==dr[k,cname]])
                        if (dr$Type[k]=='Ts,Tv2') dr$freq[k]<-(DRsites$freq.Ts.ref    [DRsites$pos==dr[k,cname]]+DRsites$freq.transv2.ref[DRsites$pos==dr[k,cname]])
                        if (dr$Type[k]=='Ts,Tv')  dr$freq[k]<-DRsites$freq.mutations.ref[DRsites$pos==dr[k,cname]]
                }
                
                #write.csv(dr, paste0("Output_all/DR/",geno[j],"/DRsites.",dname,".csv"))
                dr$freq<-as.numeric(dr$freq)
                # fixed to RAV 
                diff[[i]]<-dr[,c("ID","obs")]
                names(diff)[i]<-dname
                #mut.freq
                DR_mutfreq[[i]]<-dr[,c("ID","freq")]
                names(DR_mutfreq)[i]<-dname
        }
        
        for (i in 1:length(DR_mutfreq)) {
                colnames(DR_mutfreq[[i]])<-c("ID",paste0(names(DR_mutfreq[i])))
                colnames(diff[[i]])<-c("ID",paste0(names(diff[i])))
        }

        DR.mutated.counts<-do.call(rbind,diff.count)
        
        DR_MutFreq<-DR_mutfreq%>% purrr::reduce(full_join, by='ID')
        DR_diff<-diff %>% purrr::reduce(full_join, by='ID')
        write.csv(DR_MutFreq, paste0("Output_all/DR/RAV.MutationFreq_summary.", geno[j],".csv"))
        write.csv(DR_diff,    paste0("Output_all/DR/RAV.counts.MutFreq_summary.", geno[j],".csv"))
        
        s<-length(HCVFiles)
        #DR_MutFreq$mean<-rowMeans(DR_MutFreq[2:s+1],na.rm=T)

        ns3n <-nrow(dr[dr$Gene=="NS3",])
        ns5an<-nrow(dr[dr$Gene=="NS5A",])
        ns5bn<-nrow(dr[dr$Gene=="NS5B",])
        
        
        #count the number of patients fixed with RAV (%)
        DR_diff$total<-apply(DR_diff[2:(s+1)],1,sum, na.rm=T)
        DR_diff$Percent<-format(round(DR_diff$total/s*100, 1), nsmall=1)
        
        
        
        #create a figure
        pdf(paste0("Output_all/DR/DrugResi.allele.freq_",geno[j],".pdf"), height = 8, width = 15.5)
        ymin <- -5
        plot(0, type = "n", xlim = c(1, nrow(DR_MutFreq)), ylim = c(ymin, 0), axes = FALSE, ylab = "Frequency of RAVs", 
             xlab = "")
        axis(side = 2, at = seq(0, ymin, by=-1), labels = expression(10^0, 10^-1, 10^-2, 10^-3, 10^-4,10^-5), las = 2)
        
        n<-seq(1, by = 2, len = (nrow(DR_MutFreq)/2))

        
        genesDF<-data.frame("name"= c("NS3","NS5A","NS5B"), "Begin"= c(1,ns3n+1,(ns3n+ns5an+1)),"End"= c(ns3n,(ns3n+ns5an),nrow(DR_MutFreq)))
        for (i in 1:3){
                nvec<-seq(genesDF$Begin[i],genesDF$End[i],2)
                nvec2<-seq(genesDF$Begin[i]+1,genesDF$End[i],2)
                
                if (nvec[1]%%2==0) {v2<-nvec; v1<-nvec2}
                if (nvec[1]%%2==1) {v1<-nvec; v2<-nvec2}
                        
                abline(v=v1, col=paste0(colors[i],"66"),lty=1, lwd=16)
                abline(v=v2, col=paste0(colors[i],"1A"),lty=1, lwd=16)
        }
        
        for (i in 1:nrow(DR_MutFreq)){
                if (i<=genesDF$End[1]){
                        xjit <- rnorm(s, 0, .1)
                        points(rep(i,s)+xjit, log10(DR_MutFreq[i,2:(s+1)]),pch=16,cex=0.4,col=colors[1])}
                if (i<=genesDF$End[2]& i>=genesDF$Begin[2]){
                        xjit <- rnorm(s, 0, .1)
                        points(rep(i,s)+xjit, log10(DR_MutFreq[i,2:(s+1)]),pch=16,cex=0.4,col=colors[2])}
                if (i<=genesDF$End[3]& i>=genesDF$Begin[3]){
                        xjit <- rnorm(s, 0, .1)
                        points(rep(i,s)+xjit, log10(DR_MutFreq[i,2:(s+1)]),pch=16,cex=0.3, col=colors[3])}
        }
        
        #mtext(side = 3, at = 1:nrow(DR_MutFreq), text = paste(DR_diff$total), cex = .8)
        
        for (i in 1:nrow(dr)){
                if (dr$Type[i]=='Ts') mtext(side = 3, at = i , text = paste(DR_diff$Percent[i]),col="#117733",cex = .6,las=2)
                if (dr$Type[i]=='Tv1'|dr$Type[i]=='Tv2'|dr$Type[i]=='Tv') mtext(side = 3, at = i, text = paste(DR_diff$Percent[i]),col="#EE7733",cex = .6,las=2)
                if (dr$Type[i]=='Ts,Tv2'|dr$Type[i]=='Ts,Tv') mtext(side = 3, at = i, text = paste(DR_diff$Percent[i]),col="#0077BB",cex = .6,las=2)
        }
                                                                                  
        mtext(side = 1, at = 1:nrow(DR_MutFreq), text = paste(dr$ID), las=2,padj=0, cex = .8)
        
        abline(v=ns3n+.5,col='gray50', lwd=3)
        abline(v=ns3n+ns5an+.5,col='gray50', lwd=3 )
        
        rect(0.5,-5,ns3n+.5,-5.19,density = NULL, angle = 45,col="white",border =colors[1])
        text(ns3n/2+.5,-5.08,paste0(genesDF$name[1]),col="black", cex=.8)
        rect(ns3n+.5,-5, ns3n+ns5an+.5 ,-5.19,density = NULL, angle = 45,col="white",border =colors[2])
        text(ns3n+ns5an-ns5an/2+.5,-5.08,paste0(genesDF$name[2]),col="black", cex=.8)
        rect(ns3n+ns5an+.5,-5,nrow(dr)+.5,-5.19,density = NULL, angle = 45,col="white",border =colors[3])
        text(nrow(dr)-(nrow(dr)-ns3n-ns5an)/2+.5,-5.08,paste0(genesDF$name[3]),col="black", cex=.8)
        
        dev.off()
}
