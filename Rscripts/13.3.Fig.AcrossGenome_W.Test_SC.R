library(purrr)
library(tidyverse)

#Prep data
source("Rscripts/baseRscript.R")

###### 1) ########
selcoef_Ts<-list()
selcoef_tv1<-list()
selcoef_tv2<-list()
selcoef_tvs<-list()
selcoef_all<-list()

for (i in 1:length(Overview_fil)){
        dat<-Overview_fil[[i]]
        filename<-names(Overview_fil)[i]
        
        selcoef_Ts[[i]]<-dat[,c("pos","EstSelCoeff")] #transition mutation freq. in relative to H77.
        selcoef_tv1[[i]]<-dat[,c("pos","EstSelCoeff_trans1")] #transversion mutation freq. in relative to MajNt.
        selcoef_tv2[[i]]<-dat[,c("pos","EstSelCoeff_trans2")] #trasversion mutation freq. in relative to H77.
        selcoef_tvs[[i]]<-dat[,c("pos","EstSelCoeff_transv")] #trasversion mutation freq. in relative to H77.
        #selcoef_all[[i]]<-dat[,c("pos","EstSelCoeff_all")] #trasversion mutation freq. in relative to H77.
        
        names(selcoef_Ts)[i]<-filename
        names(selcoef_tv1)[i]<-filename
        names(selcoef_tv2)[i]<-filename
        names(selcoef_tvs)[i]<-filename
        
        
}
#assign column names for the list
for (i in 1:length(selcoef_Ts)) {
        colnames(selcoef_Ts[[i]])<-c("pos",paste0(names(selcoef_Ts[i])))
        colnames(selcoef_tv1[[i]])<-c("pos",paste0(names(selcoef_tv1[i])))
        colnames(selcoef_tv2[[i]])<-c("pos",paste0(names(selcoef_tv2[i])))
        colnames(selcoef_tvs[[i]])<-c("pos",paste0(names(selcoef_tvs[i])))
        #colnames(selcoef_all[[i]])<-c("pos",paste0(names(selcoef_all[i])))
}

TsSelcoef<-selcoef_Ts%>% purrr::reduce(full_join, by='pos')
Tv1.selcoef<-selcoef_tv1 %>% purrr::reduce(full_join, by='pos')
Tv2.selcoef<-selcoef_tv2 %>% purrr::reduce(full_join, by='pos')
Tvs.selcoef<-selcoef_tvs %>%  purrr::reduce(full_join, by='pos')
#Allselcoef<-selcoef_all %>% reduce(full_join, by='pos')

write.csv(TsSelcoef, file=paste0("Output/SelCoeff/FilteredSummary_Transition_SelCoeff_",Sys.Date(),".csv"))
write.csv(Tv1.selcoef, file=paste0("Output/SelCoeff/FilteredSummary_Transv1_SelCoeff_",Sys.Date(),".csv"))
write.csv(Tv2.selcoef, file=paste0("Output/SelCoeff/FilteredSummary_Transv2_SelCoeff_",Sys.Date(),".csv"))
write.csv(Tvs.selcoef, file=paste0("Output/SelCoeff/FilteredSummary_Tranvs_SelCoeff_",Sys.Date(),".csv"))
#write.csv(Allselcoef, file=paste0("Output/SelCoeff/FilteredSummary_Mutation_frequencies_all_",Sys.Date(),".csv"))


#total number of samples 
s<-length(Overview_fil)

#1)
SC<-TsSelcoef
#average 
SC$mean<-apply(SC[2:(s+1)], 1, mean, na.rm=T)
range(SC$mean, na.rm=T) 

nrow(SC)
n<-data.frame("pos"=c(1:nrow(SC)))
SC2<-merge(n,SC,by="pos",all.x=T)

#Trans.mut.freq[,2:ncol(Trans.mut.freq)][Trans.mut.freq[,2:ncol(Trans.mut.freq)]>=0.2] <-NA 

cols <- c("#66CCEE","#228833","#CCBB44","#EE6677","#AA3377","#4477AA","#BBBBBB")

#plot the average SNP frequency across the genome (based on H77)
pdf(paste0("Output/SummaryFig.Filtered/SelCoeff_Transitions.pdf"),width=15,height=7.5)
plot(mean~pos, data=SC2,t="n",log='y',yaxt='n',xlab='Genome position (H77)',ylab="Average selection coefficient",
     main="Selection coefficient",ylim=c(0.001,1),xlim=c(340,8500))
eaxis(side = 2, at = 10^((-0):(-(3))), cex=2)
for(i in 1:3){abline(h = 1:10 * 10^(-i), col = "gray80")}
points(mean~pos, data=SC2,pch=20,col="#EE667766",cex=0.5)

#add rolling average
roll100<-rollmean(SC2$mean, k=100)
SC2$roll100<-c(rep(NA, times=50),roll100,c(rep(NA, times=49)))
lines(roll100~pos,data=SC2, col="#AA3377")

dev.off()

########################################
### Plot SelectionCoefficients across the genome based on different mutation types 
dat<-Overview_fil[[3]]
mutationtypes<-dat[,c("pos","ref","Type","Type.tv1","Type.tv2")]
muttypes2<-dat[,c("pos","ref","Type","Type.tv1","Type.tv2","makesCpG","makesCpG.tv1","makesCpG.tv2")]

SCs<-merge(mutationtypes,SC,by='pos')

maxpos<-SC$pos[nrow(SCs)]
n<-data.frame("pos"=c(1:maxpos))
SC3<-merge(n,SCs,by="pos",all.x=T)

pdf(paste0("Output/SummaryFig.Filtered/SelCoef_based_on_Mutation_Types_1Transitions.pdf"),width=15,height=7.5)
maxnuc=8600
par(mar = c(3,5,1,2))
#selcoeffcolumn <-SC3$mean 

plot(SC3$pos[1:maxnuc],SC3$mean[1:maxnuc],
        log="y", ylab="Mean of estimated selection coefficient",cex.lab=1.4,
        yaxt="n", xlab="",xaxt='n',
        col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-5,1),xlim=c(340,maxnuc))
axis(1,at=c(seq(500,8500,by=1000)), labels=c(seq(500,8500,by=1000)))
eaxis(side = 2, at = 10^((-0):(-(5))), cex=2)

for(i in 1:5){abline(h = 1:10 * 10^(-i), col = "gray41")}

for (i in 1:maxnuc){
        c=0; co = 1
        if (is.na(SC3$Type[i])==T) next
        if (SC3$Type[i]=="stop") {c=1;p=21}
        if (SC3$Type[i]=="syn") {c=cols[3];p=21}
        if (SC3$Type[i]=="nonsyn"&SC3$ref[i]%in%c("c","g")) {c=cols[4];p=21}
        if (SC3$Type[i]=="nonsyn"&SC3$ref[i]%in%c("a","t")) {c=cols[5];p=21}
        if (c!=0) points(SC3$pos[i],SC3$mean[i],pch=p,col='gray30',lwd=0.5,
                        bg=rgb(red=col2rgb(c)[1]/255,
                        green=col2rgb(c)[2]/255,
                                blue=col2rgb(c)[3]/255,
                                maxColorValue = 1,alpha=0.8),cex=.6)
                }
        #Add legend
 legpos=300; legposV=0.53
        rect(legpos, 0.42*legposV, (legpos+1200), 1.6*legposV, density = NULL, angle = 45,col=alpha("white",1))
        points((legpos+100),legposV/0.7,pch=21,bg=1,col=1,cex=1)
        text((legpos+150),legposV/0.7,"Nonsense",adj=0, cex=1)
        points((legpos+100),legposV,pch=21,bg=cols[4],col=1,cex=1)
        text((legpos+150),legposV,"Non-syn, C/G",adj=0, cex=1)
        points((legpos+100),legposV*0.7,pch=21,bg=cols[5],col=1,cex=1)
        text((legpos+150),legposV*0.7,"Non-syn, A/T",adj=0, cex=1)
        points((legpos+100),legposV*0.49,pch=21,bg=cols[3],col=1,cex=1)
        text((legpos+150),legposV*0.49,"Syn",adj=0, cex=1)

dev.off()



#2) ############
#Transversion mutations selection coefficient figure

#Format the transversion mutation types to syn, nonsyn, mixed, and nonsesne.
SCtv1<-Tv1.selcoef
SCtv2<-Tv2.selcoef
s<-length(Overview_fil)
SCtv1$mean<-apply(SCtv1[2:(s+1)], 1, mean, na.rm=T)
SCtv2$mean<-apply(SCtv2[2:(s+1)], 1, mean, na.rm=T)
SCtv1<-SCtv1[,c("pos","mean")]
SCtv2<-SCtv2[,c("pos","mean")]
colnames(SCtv2)[2]<-"mean2"

SCtv1.1<-merge(mutationtypes,SCtv1,by='pos')
SCtv2.1<-merge(mutationtypes,SCtv2,by='pos')
SCtv<-merge(SCtv1.1,SCtv2,by='pos')

SCtv$Type.tvs<-""
for (i in 1:nrow(SCtv)){ 
        if (is.na(SCtv$Type.tv1[i])|is.na(SCtv$Type.tv2[i])) SCtv$Type.tvs<-NA
        else 
                {if (SCtv$Type.tv1[i] == SCtv$Type.tv2[i])  SCtv$Type.tvs[i]<-SCtv$Type.tv1[i]
                 if (SCtv$Type.tv1[i] != SCtv$Type.tv2[i]) SCtv$Type.tvs[i]<-'mixed'}
}        

######
#2.1 plot the summary trnasversion mutations 

maxpos<-SCtv$pos[nrow(SCtv)]
n<-data.frame("pos"=c(1:maxpos))
SCtv.1<-merge(n,SCtv,by="pos",all.x=T)


col80<-c("#FFFFFFCC","#228833CC","#CCBB44CC","#EE6677CC","#AA3377CC")

pdf(paste0("Output/SummaryFig.Filtered/SelCoef_2Transversions_mixed.pdf"),width=15,height=7.5)
maxnuc=8600
par(mar = c(3,5,1,2))

plot(SCtv.1$pos[1:maxnuc],SCtv.1$mean[1:maxnuc],
             log="y", ylab="Mean of estimated selection coefficient",cex.lab=1.4,
             yaxt="n", xlab="",xaxt='n',
             col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-3,1),xlim=c(340,maxnuc))
axis(1,at=c(seq(500,8500,by=1000)), labels=c(seq(500,8500,by=1000)))
eaxis(side = 2, at = 10^((-0):(-(3))), cex=2)
        
for(i in 1:3){abline(h = 1:10 * 10^(-i), col = "gray60")}
        
for (i in 1:maxnuc){
        if (is.na(SCtv.1$Type.tvs[i])) next
        if (SCtv.1$Type.tvs[i]=="stop") {c=1;p=21}
        if (SCtv.1$Type.tvs[i]=="syn") {c=col80[3];p=21}
        if (SCtv.1$Type.tvs[i]=="nonsyn"&SCtv.1$ref[i]%in%c("c","g")) {c=col80[4];p=21}
        if (SCtv.1$Type.tvs[i]=="nonsyn"&SCtv.1$ref[i]%in%c("a","t")) {c=col80[5];p=21}
        if (SCtv.1$Type.tvs[i]=="mixed"){c=col80[2];p=21}
        points(SCtv.1$pos[i],mean(SCtv.1$mean[i],SCtv.1$mean2[i]),pch=p,
                       col='gray30',lwd = .5, bg=c,cex=.7)
       
}
        
#Add legend
legpos=7500; legposV=0.01
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
        
maxpos<-SCtv1$pos[nrow(SCtv1.1)]
n<-data.frame("pos"=c(1:maxpos))
SCtv1.2<-merge(n,SCtv1.1,by="pos",all.x=T)
SCtv2.2<-merge(n,SCtv2.1,by="pos",all.x=T)
        
col80<-c("#FFFFFFCC","#228833CC","#CCBB44CC","#EE6677CC","#AA3377CC")

pdf(paste0("Output/SummaryFig.Filtered/SelCoef_based_on_Mutation_Types_2Transversions.pdf"),width=15,height=7.5)
maxnuc=8600
par(mar = c(3,5,1,2))
        
plot(SCtv1.2$pos[1:maxnuc],SCtv1.2$mean[1:maxnuc],
     log="y", ylab="Mean of estimated selection coefficient",cex.lab=1.4,
     yaxt="n", xlab="",xaxt='n',
     col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-3,1),xlim=c(340,maxnuc))
        axis(1,at=c(seq(500,8500,by=1000)), labels=c(seq(500,8500,by=1000)))
        eaxis(side = 2, at = 10^((-0):(-(3))), cex=2)
        
for(i in 1:3){abline(h = 1:10 * 10^(-i), col = "gray60")}
        
for (i in 1:maxnuc){
        if (is.na(SCtv1.2$Type.tv1[i])==T) next
        if (SCtv1.2$Type.tv1[i]=="stop") {c=1;p=21}
        if (SCtv1.2$Type.tv1[i]=="syn") {c=col80[3];p=21}
        if (SCtv1.2$Type.tv1[i]=="nonsyn"&SCtv1.2$ref[i]%in%c("c","g")) {c=col80[4];p=21}
        if (SCtv1.2$Type.tv1[i]=="nonsyn"&SCtv1.2$ref[i]%in%c("a","t")) {c=col80[5];p=21}
        if (is.na(SCtv2.2$Type.tv2[i])==T) next
        if (SCtv2.2$Type.tv2[i]=="stop") {c=1;p=21}
        if (SCtv2.2$Type.tv2[i]=="syn") {c=col80[3];p=21}
        if (SCtv2.2$Type.tv2[i]=="nonsyn"&SCtv2.2$ref[i]%in%c("c","g")) {c=col80[4];p=21}
        if (SCtv2.2$Type.tv2[i]=="nonsyn"&SCtv2.2$ref[i]%in%c("a","t")) {c=col80[5];p=21}

        points(SCtv1.2$pos[i],SCtv1.2$mean[i],pch=p,
              col='gray30',lwd = .5, bg=c,cex=.7)
        points(SCtv2.2$pos[i],SCtv2.2$mean[i],pch=p,
               col='gray30',lwd = .5, bg=c,cex=.7)
        }
        
        #Add legend
legpos=7500; legposV=0.01
        rect(legpos, .29*legposV, (legpos+1000), 1.1*legposV, density = NULL, angle = 45,col=alpha("white",1))
        points((legpos+100),legposV,pch=21,bg=1,col=col80[1],cex=1)
        text((legpos+150),legposV,"Nonsense",adj=0, cex=1)
        points((legpos+100),legposV*0.7,pch=21,bg=col80[4],col=1,cex=1)
        text((legpos+150),legposV*0.7,"Non-syn, C/G",adj=0, cex=1)
        points((legpos+100),legposV*0.49,pch=21,bg=col80[5],col=1,cex=1)
        text((legpos+150),legposV*0.49,"Non-syn, A/T",adj=0, cex=1)
        points((legpos+100),legposV*0.35,pch=21,bg=cols[3],col=1,cex=1)
        text((legpos+150),legposV*0.35,"Syn",adj=0, cex=1)
        
dev.off()
        




##### Test if selection coefficients /mean frequenceis differ between mutation types
#testing the mean of means

for (i in c("SC","SCtv1","SCtv2")) {
        dat<-get(i) 
        SC4<-merge(muttypes2,dat,by='pos')
        if (i=="SC") {type<-3;cpg<-6; filename<-"transitions"}
        if (i=="SCtv1") {type<-4;cpg<-7; filename<-"transv1"}
        if (i=="SCtv2") {type<-5;cpg<-8; filename<-"transv2"}

        types<-SC4[type]
        typevector1<-types=="syn"
        typevector2<-types=="nonsyn"
        typevector3<-types=="stop"
        scSyn<-SC4$mean[typevector1==T]
        scNonSyn<-SC4$mean[typevector2==T]
        scStop<-SC4$mean[typevector3==T]
        cpgmake<-SC4[cpg]
        cpgvector1<-cpgmake==0
        cpgvector2<-cpgmake==1
                
        result1<-wilcox.test(scSyn, scNonSyn,alternative = "greater", paired = FALSE) #p-value < 2.2e-16
        result2<-wilcox.test(scNonSyn,scStop,alternative = "greater", paired = FALSE) # p-value = 0.952
                
        #Test if CpG making mutation frequencies are lower than non-CpGmaking 
        scSyn.CpG<-SC4$mean[typevector1==T&cpgvector2==T]
        scSyn.nonCpG<-SC4$mean[typevector1==T&cpgvector1==T]
        scNonSyn.CpG<-SC4$mean[typevector2==T&cpgvector2==T]
        scNonSyn.nonCpG<-SC4$mean[typevector2==T&cpgvector1==T]
                
        result3<-wilcox.test(scSyn.CpG, scSyn.nonCpG,alternative = "greater", paired = FALSE) #p-value = 1
        result4<-wilcox.test(scNonSyn.CpG, scNonSyn.nonCpG,alternative = "greater", paired = FALSE) #p-value = 1 
                
        #Test whether synonymous C to T  and G to A mutations are more costly than A to G 
        CostSynNoCpGA2G<-SC4$mean[typevector1==T&cpgvector1==T&SC4$ref=="a"]
        CostSynNoCpGC2T<-SC4$mean[typevector1==T&cpgvector1==T&SC4$ref=="c"]
        CostSynNoCpGG2A<-SC4$mean[typevector1==T&cpgvector1==T&SC4$ref=="g"]
        result5<-wilcox.test(CostSynNoCpGC2T,CostSynNoCpGA2G,alternative = "greater", paired = FALSE)  #p-value < 2.2e-16
        result6<-wilcox.test(CostSynNoCpGG2A,CostSynNoCpGA2G,alternative = "greater", paired = FALSE)   #p-value < 2.2e-16
                
        #Test whether non-synonymous C to T  muts and G to A , are more costly than A to G 
        CostnonSynNoCpGA2G<-SC4$mean[typevector2==T&cpgvector1==T&SC4$ref=="a"]
        CostnonSynNoCpGC2T<-SC4$mean[typevector2==T&cpgvector1==T&SC4$ref=="c"]
        CostnonSynNoCpGG2A<-SC4$mean[typevector2==T&cpgvector1==T&SC4$ref=="g"]
        result7<-wilcox.test(CostnonSynNoCpGC2T,CostnonSynNoCpGA2G,alternative = "greater", paired = FALSE) #p-value < 2.2e-16
        result8<-wilcox.test(CostnonSynNoCpGG2A,CostnonSynNoCpGA2G,alternative = "greater", paired = FALSE) #p-value < 2.2e-16
        
        WilcoxTest.results<-data.frame(matrix(ncol=2,nrow=6))
        colnames(WilcoxTest.results)<-c("test","P.value")
        for (r in 1:6){
                result<-get(paste0('result',r))
                WilcoxTest.results$test[r]<-result[[7]]
                WilcoxTest.results$P.value[r]<-result[[3]]
                
        }
       write.csv(WilcoxTest.results,paste0("Output/SummaryStats/WilcoxTestResults_",filename,".csv"))
        
}
