library(purrr)
library(tidyverse)
library(zoo)


source("Rscripts/baseRscript.R")

filtered<-list()
fnames<-c("Ts", "Ts_NA", "Ts_zero")
for (i in 1:3){
        filename<-fnames[i]
        df<-read.csv(paste0("Output/Mut.freq.filtered/Summary_",filename,".Q35.csv"),stringsAsFactors = F)
        df<-df[,-1]
        filtered[[i]]<-df
        names(filtered)[i]<-filename
}


genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)

gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector



### Plot mutation freq. across the genome based on the mutation types 

mfs<-filtered[[3]]
maxpos<-mfs$pos[nrow(mfs)]

n<-data.frame("pos"=c(1:maxpos))
MF3<-merge(n,mfs,by="pos",all.x=T)

#pdf(paste0("Output/SummaryFig.Filtered/Ts_MutFreq_based_on_Mutation_Types.pdf"),width=15,height=7.5)
pdf(paste0("Output/SummaryFig.Filtered/Ts_NA_MutFreq_based_on_Mutation_Types.pdf"),width=15,height=7.5)
pdf(paste0("Output/SummaryFig.Filtered/Ts_zero_MutFreq_based_on_Mutation_Types.pdf"),width=15,height=7.5)

maxnuc=maxpos
par(mar = c(3,5,1,2))
#selcoeffcolumn <-SC3$mean 

plot(MF3$pos[1:maxnuc],MF3$mean[1:maxnuc],
        log="y", ylab="Mean mutation frequency",cex.lab=1.4,
        yaxt="n", xlab="",xaxt='n',
        col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-4,0.1),xlim=c(340,maxnuc))
axis(1,at=c(seq(500,8500,by=1000)), labels=c(seq(500,8500,by=1000)))
eaxis(side = 2, at = 10^((-1):(-(4))), cex=2)

for(i in 2:4){abline(h = 1:10 * 10^(-i), col = "gray80")}

for (i in 1:maxnuc){
        c=0; co = 1
        if (is.na(MF3$Type[i])==T) next
        if (MF3$Type[i]=="stop") {c=1;p=21}
        if (MF3$Type[i]=="syn") {c=cols[3];p=21}
        if (MF3$Type[i]=="nonsyn"&MF3$ref[i]%in%c("c","g")) {c=cols[4];p=21}
        if (MF3$Type[i]=="nonsyn"&MF3$ref[i]%in%c("a","t")) {c=cols[5];p=21}
        if (c!=0) points(MF3$pos[i],MF3$mean[i],pch=p,col='gray30',lwd=0.5,
                        bg=rgb(red=col2rgb(c)[1]/255,
                        green=col2rgb(c)[2]/255,
                                blue=col2rgb(c)[3]/255,
                                maxColorValue = 1,alpha=0.8),cex=.6)
                }
        #Add legend
 legpos=300; legposV=0.1
        rect(legpos, 0.42*legposV, (legpos+1000), 1.05*legposV, density = NULL, angle = 45,col=alpha("white",1))
        points((legpos+100),legposV*0.9,pch=21,bg=1,col=1,cex=1)
        text((legpos+150),legposV*0.9,"Nonsense",adj=0, cex=1)
        points((legpos+100),legposV*0.74,pch=21,bg=cols[4],col=1,cex=1)
        text((legpos+150),legposV*0.74,"Non-syn, C/G",adj=0, cex=1)
        points((legpos+100),legposV*0.6,pch=21,bg=cols[5],col=1,cex=1)
        text((legpos+150),legposV*0.6,"Non-syn, A/T",adj=0, cex=1)
        points((legpos+100),legposV*0.49,pch=21,bg=cols[3],col=1,cex=1)
        text((legpos+150),legposV*0.49,"Syn",adj=0, cex=1)
        
ylow<-0.00025
for (j in 2:(nrow(genes)-1)){
        xleft<-genes$start[j]
        xright<-genes$start[j+1]
        if (j==5){
                rect(xleft,ylow,xright,ylow*1.3,density = NULL, angle = 45,col="white",border ="#AA3377CC")
                text(xleft+80,1.15*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        }
        else if ((j==4|j==6|j==9)){
                rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="#AA3377CC")
                text(xleft+80, 1.4*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                #mtext(paste0(genes$Gene[j]),side= 1, line=-0.1, at= xleft+80, col="black", cex=0.8)
        }
        else{rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="#AA3377CC")
                text(xleft+200,1.15*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)}
}
        
        
box()
dev.off()









#2) ############
#Transversion mutations selection coefficient figure

#Format the transversion mutation types to syn, nonsyn, mixed, and nonsesne.
tv1<-Tv1.MutFreq
tv2<-Tv2.MutFreq
s<-length(Overview_fil)
tv1$mean<-apply(tv1[2:(s+1)], 1, mean, na.rm=T)
tv2$mean<-apply(tv2[2:(s+1)], 1, mean, na.rm=T)
tv1$NofNA<-apply(tv1[2:(s+1)], 1, function(x)sum(is.na(x)))
tv2$NofNA<-apply(tv2[2:(s+1)], 1, function(x)sum(is.na(x)))

tv1.1<-tv1[which(tv1$NofNA<0.5*s),]
tv2.1<-tv2[which(tv2$NofNA<0.5*s),]
tv1.2<-tv1.1[,c("pos","mean")]
tv2.2<-tv2.1[,c("pos","mean")]
colnames(tv2.2)[2]<-"mean2"

tv1.3<-merge(mutationtypes,tv1.2,by='pos')
tv2.3<-merge(mutationtypes,tv2.2,by='pos')
MFtv<-merge(tv1.3,tv2.2,by='pos')

MFtv$Type.tvs<-""
for (i in 1:nrow(MFtv)){ 
        if (is.na(MFtv$Type.tv1[i])|is.na(MFtv$Type.tv2[i])) MFtv$Type.tvs<-NA
        else 
                {if (MFtv$Type.tv1[i] == MFtv$Type.tv2[i]) MFtv$Type.tvs[i]<-MFtv$Type.tv1[i]
                 if (MFtv$Type.tv1[i] != MFtv$Type.tv2[i]) MFtv$Type.tvs[i]<-'mixed'}
}        

######
#2.1 plot the summary trnasversion mutations 

maxpos<-MFtv$pos[nrow(MFtv)]
n<-data.frame("pos"=c(1:maxpos))
MFtv.1<-merge(n,MFtv,by="pos",all.x=T)


col80<-c("#FFFFFFCC","#228833CC","#CCBB44CC","#EE6677CC","#AA3377CC")

pdf(paste0("Output/SummaryFig.Filtered/Transv_Mut_Freq_based_on_Mutation_Types.pdf"),width=15,height=7.5)
maxnuc=8600
par(mar = c(3,5,1,2))

plot(MFtv.1$pos[1:maxnuc],MFtv.1$mean[1:maxnuc],
             log="y", ylab="Mean mutation frequency",cex.lab=1.4,
             yaxt="n", xlab="",xaxt='n',
             col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-6,0.1),xlim=c(340,maxnuc))
axis(1,at=c(seq(500,8500,by=1000)), labels=c(seq(500,8500,by=1000)))
eaxis(side = 2, at = 10^((-1):(-(6))), cex=2)
        
for(i in 1:6){abline(h = 1:10 * 10^(-i), col = "gray60")}
        
for (i in 1:maxnuc){
        if (is.na(MFtv.1$Type.tvs[i])) next
        if (MFtv.1$Type.tvs[i]=="stop") {c=1;p=21}
        if (MFtv.1$Type.tvs[i]=="syn") {c=col80[3];p=21}
        if (MFtv.1$Type.tvs[i]=="nonsyn"&MFtv.1$ref[i]%in%c("c","g")) {c=col80[4];p=21}
        if (MFtv.1$Type.tvs[i]=="nonsyn"&MFtv.1$ref[i]%in%c("a","t")) {c=col80[5];p=21}
        if (MFtv.1$Type.tvs[i]=="mixed"){c=col80[2];p=21}
        points(MFtv.1$pos[i],mean(MFtv.1$mean[i],MFtv.1$mean2[i]),pch=p,
                       col='gray30',lwd = .5, bg=c,cex=.7)
       
}
        
#Add legend
legpos=7700; legposV=0.05
rect(legpos, .29*legposV, (legpos+1000), 1.7*legposV, density = NULL, angle = 45,col=alpha("white",1))
points((legpos+100),legposV/0.7,pch=21,bg=1,col=col80[1],cex=1)
text((legpos+150),legposV/0.7,"Nonsense",adj=0, cex=1)
points((legpos+100),legposV,pch=21,bg=col80[4],col=1,cex=1)
text((legpos+150),legposV,"Non-syn, C/G",adj=0, cex=1)
points((legpos+100),legposV*0.7,pch=21,bg=col80[5],col=1,cex=1)
text((legpos+150),legposV*0.7,"Non-syn, A/T",adj=0, cex=1)
points((legpos+100),legposV*0.49,pch=21,bg=col80[2],col=1,cex=1)
text((legpos+150),legposV*0.49,"Mixed",adj=0, cex=1)
points((legpos+100),legposV*0.35,pch=21,bg=cols[3],col=1,cex=1)
text((legpos+150),legposV*0.35,"Syn",adj=0, cex=1)

dev.off()

#2.2. plot trnasversion mutations separately (Type 1 and Type 2)

        
maxpos<-tv1.3$pos[nrow(tv1.3)]
n<-data.frame("pos"=c(1:maxpos))
tv1.4<-merge(n,tv1.3,by="pos",all.x=T)
tv2.4<-merge(n,tv2.3,by="pos",all.x=T)
        
col80<-c("#FFFFFFCC","#228833CC","#CCBB44CC","#EE6677CC","#AA3377CC")

pdf(paste0("Output/SummaryFig.Filtered/Transv_MutFreq_Mutation_Types_separate.pdf"),width=15,height=7.5)
maxnuc=8600
par(mar = c(3,5,1,2))
        
plot(tv1.4$pos[1:maxnuc],tv1.4$mean[1:maxnuc],
     log="y", ylab="Average mutation frequency",cex.lab=1.4,
     yaxt="n", xlab="",xaxt='n',
     col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-6,.1),xlim=c(340,maxnuc))
        axis(1,at=c(seq(500,8500,by=1000)), labels=c(seq(500,8500,by=1000)))
        eaxis(side = 2, at = 10^((-0):(-(6))), cex=2)
        
for(i in 1:6){abline(h = 1:10 * 10^(-i), col = "gray60")}
        
for (i in 1:maxnuc){
        if (is.na(tv1.4$Type.tv1[i])==T) next
        if (tv1.4$Type.tv1[i]=="stop") {c=1;p=21}
        if (tv1.4$Type.tv1[i]=="syn") {c=col80[3];p=21}
        if (tv1.4$Type.tv1[i]=="nonsyn"&tv1.4$ref[i]%in%c("c","g")) {c=col80[4];p=21}
        if (tv1.4$Type.tv1[i]=="nonsyn"&tv1.4$ref[i]%in%c("a","t")) {c=col80[5];p=21}
        if (is.na(tv2.4$Type.tv2[i])==T) next
        if (tv2.4$Type.tv2[i]=="stop") {c=1;p=21}
        if (tv2.4$Type.tv2[i]=="syn") {c=col80[3];p=21}
        if (tv2.4$Type.tv2[i]=="nonsyn"&tv2.4$ref[i]%in%c("c","g")) {c=col80[4];p=21}
        if (tv2.4$Type.tv2[i]=="nonsyn"&tv2.4$ref[i]%in%c("a","t")) {c=col80[5];p=21}

        points(tv1.4$pos[i],tv1.4$mean[i],pch=p,
              col='gray30',lwd = .5, bg=c,cex=.7)
        points(tv2.4$pos[i],tv2.4$mean[i],pch=p,
               col='gray30',lwd = .5, bg=c,cex=.7)
        }
        
        #Add legend
        legpos=7700; legposV=0.05
        rect(legpos, .29*legposV, (legpos+1000), 1.2*legposV, density = NULL, angle = 45,col=alpha("white",1))
        points((legpos+100),legposV,pch=21,bg=1,col=col80[1],cex=1)
        text((legpos+150),legposV,"Nonsense",adj=0, cex=1)
        points((legpos+100),legposV*0.7,pch=21,bg=col80[4],col=1,cex=1)
        text((legpos+150),legposV*0.7,"Non-syn, C/G",adj=0, cex=1)
        points((legpos+100),legposV*0.49,pch=21,bg=col80[5],col=1,cex=1)
        text((legpos+150),legposV*0.49,"Non-syn, A/T",adj=0, cex=1)
        points((legpos+100),legposV*0.35,pch=21,bg=cols[3],col=1,cex=1)
        text((legpos+150),legposV*0.35,"Syn",adj=0, cex=1)
dev.off()
        




##### Test if mean frequenceis differ between mutation types
#testing the mean of means
#1.Transition mutation selection coefficients 

for (i in c("t2","tv1.2","tv2.2")) {
        dat<-get(i) 
        FreqData<-merge(muttypes2,dat,by='pos')
        if (i=="t2") {type<-3;cpg<-6; filename<-"transitions"}
        if (i=="tv1.2") {type<-4;cpg<-7; filename<-"transv1"}
        if (i=="tv2.2") {type<-5;cpg<-8; filename<-"transv2"}

        types<-FreqData[type]
        typevector1<-types=="syn"
        typevector2<-types=="nonsyn"
        typevector3<-types=="stop"
        scSyn<-FreqData$mean[typevector1==T]
        scNonSyn<-FreqData$mean[typevector2==T]
        scStop<-FreqData$mean[typevector3==T]
        cpgmake<-FreqData[cpg]
        cpgvector1<-cpgmake==0
        cpgvector2<-cpgmake==1
                
        result1<-wilcox.test(scSyn, scNonSyn,alternative = "greater", paired = FALSE) #p-value < 2.2e-16
        result2<-wilcox.test(scNonSyn,scStop,alternative = "greater", paired = FALSE) # p-value = 0.952
                
        #Test if CpG making mutation frequencies are lower than non-CpGmaking 
        scSyn.CpG<-FreqData$mean[typevector1==T&cpgvector2==T]
        scSyn.nonCpG<-FreqData$mean[typevector1==T&cpgvector1==T]
        scNonSyn.CpG<-FreqData$mean[typevector2==T&cpgvector2==T]
        scNonSyn.nonCpG<-FreqData$mean[typevector2==T&cpgvector1==T]
                
        result3<-wilcox.test(scSyn.CpG, scSyn.nonCpG,alternative = "greater", paired = FALSE) #p-value = 1
        result4<-wilcox.test(scNonSyn.CpG, scNonSyn.nonCpG,alternative = "greater", paired = FALSE) #p-value = 1 
                
        #Test whether synonymous C to T  and G to A mutations are more costly than A to G 
        CostSynNoCpGA2G<-FreqData$mean[typevector1==T&cpgvector1==T&FreqData$ref=="a"]
        CostSynNoCpGC2T<-FreqData$mean[typevector1==T&cpgvector1==T&FreqData$ref=="c"]
        CostSynNoCpGG2A<-FreqData$mean[typevector1==T&cpgvector1==T&FreqData$ref=="g"]
        result5<-wilcox.test(CostSynNoCpGC2T,CostSynNoCpGA2G,alternative = "greater", paired = FALSE)  #p-value < 2.2e-16
        result6<-wilcox.test(CostSynNoCpGG2A,CostSynNoCpGA2G,alternative = "greater", paired = FALSE)   #p-value < 2.2e-16
                
        #Test whether non-synonymous C to T  muts and G to A , are more costly than A to G 
        CostnonSynNoCpGA2G<-FreqData$mean[typevector2==T&cpgvector1==T&FreqData$ref=="a"]
        CostnonSynNoCpGC2T<-FreqData$mean[typevector2==T&cpgvector1==T&FreqData$ref=="c"]
        CostnonSynNoCpGG2A<-FreqData$mean[typevector2==T&cpgvector1==T&FreqData$ref=="g"]
        result7<-wilcox.test(CostnonSynNoCpGC2T,CostnonSynNoCpGA2G,alternative = "greater", paired = FALSE) #p-value < 2.2e-16
        result8<-wilcox.test(CostnonSynNoCpGG2A,CostnonSynNoCpGA2G,alternative = "greater", paired = FALSE) #p-value < 2.2e-16
        
        WilcoxTest.results<-data.frame(matrix(ncol=2,nrow=6))
        colnames(WilcoxTest.results)<-c("test","P.value")
        for (r in 1:6){
                result<-get(paste0('result',r))
                WilcoxTest.results$test[r]<-result[[7]]
                WilcoxTest.results$P.value[r]<-result[[3]]
                
        }
       write.csv(WilcoxTest.results,paste0("Output/SummaryStats/WilcoxTestResults_MutFreq_",filename,".csv"))
        
}
