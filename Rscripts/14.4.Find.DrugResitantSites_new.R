library(tidyverse)
source("Rscripts/baseRscript.R")
colors=c("#44AA99","#0077BB","#CC6677" )


#nucleotide positions of known drug resistant mutations 
drsites<-read.csv("Data/HCV_drugresistance_Geno.csv",stringsAsFactors = F)
#attach the original position info
geno1A<-read.csv("Output_all/Filtered/Ts_summary_metadata.1A.csv",stringsAsFactors = F, row.names = 1)
for (i in 1:nrow(drsites)) {drsites$pos.1A[i]<-geno1A$org.pos.1A[which(geno1A$merged.pos==drsites$merged.pos[i])]}
#select only 1A
DR<-drsites[drsites$genotype=="1A",]
write.csv(DR, "Output1A/DrugRes/RAV_Table.csv")



#######
#read the modified table of RAV info
DR<-read.csv("Output1A/DrugRes/RAV_Table_updated.csv",stringsAsFactors = F)
#create an id column
for (i in 1:nrow(DR)){
        if (DR$Need_both[i]=="y") { 
                if (DR$extra[i]=="n") DR$ID[i]<- paste0(DR$Name[i],'.',DR$merged.pos[i])
                if (DR$extra[i]=="y") DR$ID[i]<- paste (DR$Name[i],DR$merged.pos[i],DR$Type[i], sep = ".")}
        if (DR$Need_both[i]=="n") { 
                if (DR$extra[i]=="n") DR$ID[i]<-DR$Name[i]
                if (DR$extra[i]=="y") DR$ID[i]<- paste (DR$Name[i],DR$merged.pos[i],DR$Type[i], sep = ".")}
}

HCVFiles<-list.files("Output1A/Overview2/", pattern="overview2.csv")

DR_mutfreq<-data.frame(ID=factor(DR$ID, levels=c(DR$ID)))
Diff<-data.frame(ID=factor(DR$ID,levels=c(DR$ID)))
diff.count<-list()

for (i in 1:length(HCVFiles)){ 
        df<-read.csv(paste0("Output1A/Overview2/",HCVFiles[i]),stringsAsFactors=FALSE, row.names = 1)
        dname<-substr(paste(HCVFiles[i]),start=1,stop=7)
        dr<-DR
        cname<-"pos.1A"
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
                if (dr$Type[k]=='Ts')     dr$freq[k]<-DRsites$freq.Ts.ref     [DRsites$pos==dr[k,cname]]
        }
        
        #write.csv(dr, paste0("Output_all/DR/",geno[j],"/DRsites.",dname,".csv"))
        dr$freq<-as.numeric(dr$freq)
        
        #mut.freq
        dr1<-dr[,c("ID","freq")]
        colnames(dr1)[2]<-dname
        DR_mutfreq<-merge(DR_mutfreq,dr1, by="ID")   
        
        # no of samples fixed to RAV 
        dr2<-dr[,c("ID","obs")]
        colnames(dr2)[2]<-dname
        Diff<-merge(Diff, dr2, by="ID")
}

DR.mutated.counts<-do.call(rbind,diff.count)


DR_mutfreq<-DR_mutfreq[order(DR_mutfreq$ID),]
DR_diff<-Diff[order(Diff$ID),]

write.csv(DR_mutfreq, "Output1A/DrugRes/RAV.MutationFreq_summary.1A.csv")

s<-length(HCVFiles)
ns3n <-nrow(DR[DR$Gene=="NS3",])
ns5an<-nrow(DR[DR$Gene=="NS5A",])
ns5bn<-nrow(DR[DR$Gene=="NS5B",])

#count the number of patients fixed with RAV (%)

#count the number of non-NA per RAV
DR_diff$NonNA_count<-apply(DR_mutfreq[,2:(s+1)],1, function(x) sum(!is.na(x)))
DR_diff$total<-apply(DR_diff[2:(s+1)],1,sum, na.rm=T)
DR_diff$Percent<-format(round(DR_diff$total/DR_diff$NonNA_count*100, 1), nsmall=1)
write.csv(DR_diff, "Output1A/DrugRes/RAV.counts.MutFreq_summary.1A.csv")



#create a figure

par(mar=c(6, 4, 4, 2) + 0.1)
pdf(paste0("Output1A/DrugRes/RAVs.mf.all.pdf"), height = 8, width = 15.5)
ymin <- -5
plot(0, type = "n", xlim = c(1, nrow(DR_mutfreq)), ylim = c(ymin, 0), axes = FALSE, ylab = "Frequency of RAVs", 
     xlab = "")
axis(side = 2, at = seq(0, ymin, by=-1), labels = expression(10^0, 10^-1, 10^-2, 10^-3, 10^-4,10^-5), las = 2)

n<-seq(1, by = 2, len = (nrow(DR_mutfreq)/2))


genesDF<-data.frame("name"= c("NS3","NS5A","NS5B"), "Begin"= c(1,ns3n+1,(ns3n+ns5an+1)),"End"= c(ns3n,(ns3n+ns5an),nrow(DR_mutfreq)))
for (i in 1:3){
        nvec<-seq(genesDF$Begin[i],genesDF$End[i],2)
        nvec2<-seq(genesDF$Begin[i]+1,genesDF$End[i],2)
        
        if (nvec[1]%%2==0) {v2<-nvec; v1<-nvec2}
        if (nvec[1]%%2==1) {v1<-nvec; v2<-nvec2}
        
        abline(v=v1, col=paste0(colors[i],"66"),lty=1, lwd=16)
        abline(v=v2, col=paste0(colors[i],"1A"),lty=1, lwd=16)
}

for (i in 1:nrow(DR_mutfreq)){
        if (i<=genesDF$End[1]){
                xjit <- rnorm(s, 0, .1)
                points(rep(i,s)+xjit, log10(DR_mutfreq[i,2:(s+1)]),pch=16,cex=0.4,col=colors[1])}
        if (i<=genesDF$End[2]& i>=genesDF$Begin[2]){
                xjit <- rnorm(s, 0, .1)
                points(rep(i,s)+xjit, log10(DR_mutfreq[i,2:(s+1)]),pch=16,cex=0.4,col=colors[2])}
        if (i<=genesDF$End[3]& i>=genesDF$Begin[3]){
                xjit <- rnorm(s, 0, .1)
                points(rep(i,s)+xjit, log10(DR_mutfreq[i,2:(s+1)]),pch=16,cex=0.3, col=colors[3])}
}

#mtext(side = 3, at = 1:nrow(DR_mutfreq), text = paste(DR_diff$total), cex = .8)

for (i in 1:nrow(DR)){
        if (DR$Type[i]=='Ts') mtext(side = 3, at = i , text = paste(DR_diff$Percent[i]),col="#117733",cex = .6,las=2)
        if (DR$Type[i]=='Tv1'|DR$Type[i]=='Tv2') mtext(side = 3, at = i, text = paste(DR_diff$Percent[i]),col="#EE7733",cex = .6,las=2)
}

mtext(side = 1, at = 1:nrow(DR_mutfreq), text = paste(dr$ID), las=2,padj=0, cex = .8)

abline(v=ns3n+.5,col='gray50', lwd=3)
abline(v=ns3n+ns5an+.5,col='gray50', lwd=3 )

rect(0.5,-5,ns3n+.5,-5.19,density = NULL, angle = 45,col="white",border =colors[1])
text(ns3n/2+.5,-5.08,paste0(genesDF$name[1]),col="black", cex=.8)
rect(ns3n+.5,-5, ns3n+ns5an+.5 ,-5.19,density = NULL, angle = 45,col="white",border =colors[2])
text(ns3n+ns5an-ns5an/2+.5,-5.08,paste0(genesDF$name[2]),col="black", cex=.8)
rect(ns3n+ns5an+.5,-5,nrow(dr)+.5,-5.19,density = NULL, angle = 45,col="white",border =colors[3])
text(nrow(dr)-(nrow(dr)-ns3n-ns5an)/2+.5,-5.08,paste0(genesDF$name[3]),col="black", cex=.8)

dev.off()


#####
## Calculate selection coefficients for RAV sites
## this inlucdes transversion so can't use the SCs summary

HCVFiles3<-list.files("Output1A/Overview3/", pattern="overview3.csv")

#read the modified table of RAV info
dr1a<-read.csv("Output1A/DrugRes/RAV_Table_updated.csv",stringsAsFactors = F)

#create an id column
for (i in 1:nrow(dr1a)){
        if (dr1a$Need_both[i]=="y") { 
                if (dr1a$extra[i]=="n") dr1a$ID[i]<- paste0(dr1a$Name[i],'.',dr1a$merged.pos[i])
                if (dr1a$extra[i]=="y") dr1a$ID[i]<- paste (dr1a$Name[i],dr1a$merged.pos[i],dr1a$Type[i], sep = ".")}
        if (dr1a$Need_both[i]=="n") { 
                if (dr1a$extra[i]=="n") dr1a$ID[i]<-dr1a$Name[i]
                if (dr1a$extra[i]=="y") dr1a$ID[i]<- paste (dr1a$Name[i],dr1a$merged.pos[i],dr1a$Type[i], sep = ".")}
}


dr.mf2<-data.frame(ID=dr1a$ID)
for (i in 1:length(HCVFiles3)){ 
        df<-read.csv(paste0("Output1A/Overview3/",HCVFiles3[i]),stringsAsFactors=FALSE, row.names = 1)
        dname<-substr(paste(HCVFiles3[i]),start=1,stop=7)
        dr<-dr1a
        cname<-"pos.1A"
        DRsites<-df[df$pos %in% dr[,cname],]
        
        #obtain the mutation freq.
        dr$freq<-""
        for (k in 1:nrow(dr)){
                if (dr$Type[k]=='Tv1')    dr$freq[k]<-DRsites$freq.transv1.ref[DRsites$pos==dr[k,cname]]
                if (dr$Type[k]=='Tv2')    dr$freq[k]<-DRsites$freq.transv2.ref[DRsites$pos==dr[k,cname]]
                if (dr$Type[k]=='Ts')     dr$freq[k]<-DRsites$freq.Ts.ref     [DRsites$pos==dr[k,cname]]
        }
        
        #write.csv(dr, paste0("Output_all/DR/",geno[j],"/DRsites.",dname,".csv"))
        dr$freq<-as.numeric(dr$freq)
        dr<-dr[,c("ID","freq")]
        colnames(dr)[2]<-dname
        dr.mf2<-merge(dr.mf2,dr, by="ID")

        
}

write.csv(dr.mf2, "Output1A/DrugRes/RAV.MF_summary.1A.Filtered.csv")


dr.mf2$mean<-rowMeans(dr.mf2[,2:196],na.rm=T)
dr_sc<-dr.mf2[,c("ID","mean")]
dr_sc<-merge(dr_sc, dr1a, by="ID")

# attach the mut rates info:
mutrates<-read.csv("Data/Geller.mutation.rates.csv")

dr_sc$MR[dr_sc$ref=="a"&dr_sc$Type=="Ts"]<-mutrates$mut.rate[mutrates$mutations=="AG"]
dr_sc$MR[dr_sc$ref=="c"&dr_sc$Type=="Ts"]<-mutrates$mut.rate[mutrates$mutations=="CU"]
dr_sc$MR[dr_sc$ref=="g"&dr_sc$Type=="Ts"]<-mutrates$mut.rate[mutrates$mutations=="GA"]
dr_sc$MR[dr_sc$ref=="t"&dr_sc$Type=="Ts"]<-mutrates$mut.rate[mutrates$mutations=="UC"]

dr_sc$MR[dr_sc$ref=="a"&dr_sc$Type=="Tv1"]<-mutrates$mut.rate[mutrates$mutations=="AC"]
dr_sc$MR[dr_sc$ref=="c"&dr_sc$Type=="Tv1"]<-mutrates$mut.rate[mutrates$mutations=="CA"]
dr_sc$MR[dr_sc$ref=="g"&dr_sc$Type=="Tv1"]<-mutrates$mut.rate[mutrates$mutations=="GC"]
dr_sc$MR[dr_sc$ref=="t"&dr_sc$Type=="Tv1"]<-mutrates$mut.rate[mutrates$mutations=="UA"]

dr_sc$MR[dr_sc$ref=="a"&dr_sc$Type=="Tv2"]<-mutrates$mut.rate[mutrates$mutations=="AU"]
dr_sc$MR[dr_sc$ref=="c"&dr_sc$Type=="Tv2"]<-mutrates$mut.rate[mutrates$mutations=="CG"]
dr_sc$MR[dr_sc$ref=="g"&dr_sc$Type=="Tv2"]<-mutrates$mut.rate[mutrates$mutations=="GU"]
dr_sc$MR[dr_sc$ref=="t"&dr_sc$Type=="Tv2"]<-mutrates$mut.rate[mutrates$mutations=="UG"]


dr_sc$EstSC<-""
for (j in 1:nrow(dr_sc)){
        dr_sc$EstSC[j] <- EstimatedS(dr_sc$MR[j],dr_sc$mean[j])
}
dr_sc$EstSC<-as.numeric(dr_sc$EstSC)
dr_sc$ID<-factor(dr_sc$ID, levels = c(DR$ID))
dr_sc<-dr_sc[order(dr_sc$ID),]
write.csv(dr_sc, "Output1A/DrugRes/RAV.SC_summary.1A.csv")

dr_sc<-read.csv("Output1A/DrugRes/RAV.SC_summary.1A.csv")

#min(dr_sc$EstSC)
#plot(dr_sc$EstSC~dr_sc$ID)

pdf(paste0("Output1A/DrugRes/RAVs.Figure2.pdf"), height = 11, width = 15.5)
layout(matrix(c(1,2), byrow = TRUE), heights=c(1,3))
par(mar=c(0, 4, 4, .5))

genesDF<-data.frame("name"= c("NS3","NS5A","NS5B"), "Begin"= c(1,ns3n+1,(ns3n+ns5an+1)),"End"= c(ns3n,(ns3n+ns5an),nrow(DR_mutfreq)))

ymin1=-4
plot(0, type = "n", xlim = c(1, nrow(dr_sc)), ylim = c(ymin1, -1), axes = FALSE, ylab = "Cost", xlab = "")
axis(side = 2, at = seq(-1, ymin1, by=-1), labels = expression(10^-1, 10^-2, 10^-3, 10^-4), las = 2)
abline(v=c((0:nrow(dr_sc))+0.5), col="gray60", lty=1, lwd=.1)
#abline(h=-4, col="gray50")
segments(i,-3.9,i,log10(dr_sc$EstSC[i]), lwd=10, col=colors[1], lend=2)

for (i in 1:nrow(dr_sc)){
        if (i<=genesDF$End[1]) segments(i,ymin1,i,log10(dr_sc$EstSC[i]), lwd=10, col=colors[1],lend=2)
        if (i<=genesDF$End[2]& i>=genesDF$Begin[2]) segments(i,ymin1,i,log10(dr_sc$EstSC[i]), lwd=10, col=colors[2],lend=2)
        if (i<=genesDF$End[3]& i>=genesDF$Begin[3]) segments(i,ymin1,i,log10(dr_sc$EstSC[i]), lwd=10, col=colors[3],lend=2)
}


#for (i in 1:nrow(dr_sc)){
#        if (i<=genesDF$End[1]) points(i,log10(dr_sc$EstSC[i]),pch=15,cex=1,col=colors[1])
#        if (i<=genesDF$End[2]& i>=genesDF$Begin[2]) points(i,log10(dr_sc$EstSC[i]),pch=15,cex=1,col=colors[2])
#        if (i<=genesDF$End[3]& i>=genesDF$Begin[3]) points(i,log10(dr_sc$EstSC[i]),pch=15,cex=1,col=colors[3])
#}
par(mar=c(6, 4, 4, .5))

ymin <- -5
plot(0, type = "n", xlim = c(1, nrow(DR_mutfreq)), ylim = c(ymin, 0), axes = FALSE, ylab = "Frequency of RAVs", 
     xlab = "")
axis(side = 2, at = seq(0, ymin, by=-1), labels = expression(10^0, 10^-1, 10^-2, 10^-3, 10^-4,10^-5), las = 2)

n<-seq(1, by = 2, len = (nrow(DR_mutfreq)/2))

for (i in 1:3){
        nvec<-seq(genesDF$Begin[i],genesDF$End[i],2)
        nvec2<-seq(genesDF$Begin[i]+1,genesDF$End[i],2)
        
        if (nvec[1]%%2==0) {v2<-nvec; v1<-nvec2}
        if (nvec[1]%%2==1) {v1<-nvec; v2<-nvec2}
        
        abline(v=v1, col=paste0(colors[i],"66"),lty=1, lwd=16)
        abline(v=v2, col=paste0(colors[i],"1A"),lty=1, lwd=16)
}

for (i in 1:nrow(DR_mutfreq)){
        if (i<=genesDF$End[1]){
                xjit <- rnorm(s, 0, .1)
                points(rep(i,s)+xjit, log10(DR_mutfreq[i,2:(s+1)]),pch=16,cex=0.4,col=colors[1])}
        if (i<=genesDF$End[2]& i>=genesDF$Begin[2]){
                xjit <- rnorm(s, 0, .1)
                points(rep(i,s)+xjit, log10(DR_mutfreq[i,2:(s+1)]),pch=16,cex=0.4,col=colors[2])}
        if (i<=genesDF$End[3]& i>=genesDF$Begin[3]){
                xjit <- rnorm(s, 0, .1)
                points(rep(i,s)+xjit, log10(DR_mutfreq[i,2:(s+1)]),pch=16,cex=0.3, col=colors[3])}
}

#mtext(side = 3, at = 1:nrow(DR_mutfreq), text = paste(DR_diff$total), cex = .8)

for (i in 1:nrow(DR)){
        if (DR$Type[i]=='Ts') mtext(side = 3, at = i , text = paste(DR_diff$Percent[i]),col="#117733",cex = .6,las=2)
        if (DR$Type[i]=='Tv1'|DR$Type[i]=='Tv2') mtext(side = 3, at = i, text = paste(DR_diff$Percent[i]),col="#EE7733",cex = .6,las=2)
}

mtext(side = 1, at = 1:nrow(DR_mutfreq), text = paste(dr$ID), las=2,padj=0, cex = .8)

abline(v=ns3n+.5,col='gray50', lwd=3)
abline(v=ns3n+ns5an+.5,col='gray50', lwd=3 )

rect(0.5,-5,ns3n+.5,-5.19,density = NULL, angle = 45,col="white",border =colors[1])
text(ns3n/2+.5,-5.08,paste0(genesDF$name[1]),col="black", cex=.8)
rect(ns3n+.5,-5, ns3n+ns5an+.5 ,-5.19,density = NULL, angle = 45,col="white",border =colors[2])
text(ns3n+ns5an-ns5an/2+.5,-5.08,paste0(genesDF$name[2]),col="black", cex=.8)
rect(ns3n+ns5an+.5,-5,nrow(dr)+.5,-5.19,density = NULL, angle = 45,col="white",border =colors[3])
text(nrow(dr)-(nrow(dr)-ns3n-ns5an)/2+.5,-5.08,paste0(genesDF$name[3]),col="black", cex=.8)

dev.off()

####
#SC for RAV sites vs. none.
drsites<-dr_sc[,c("ID","pos.1A","EstSC","Type","Gene","MR")]


mean(dr_sc$EstSC[dr_sc$Type=='syn'], na.rm=T)
