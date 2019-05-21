library(purrr)
library(tidyverse)
library(zoo)


#Script to analyse the frequency data and associate with features. 

source("Rscripts/baseRscript.R")

HCVFiles_overview2<-list.files("Output/Overview_output_Ref/",pattern="overview2.csv")

#for Overview with Reference Sequence:
Overview_sum_ref<-list()
for (i in 1:length(HCVFiles_overview2)){ 
        overviews<-read.csv(paste0("Output/Overview_output_Ref/",HCVFiles_overview2[i]),stringsAsFactors=FALSE)
        overviews<-overviews[,-1]
        Overview_sum_ref[[i]]<-overviews
        names(Overview_sum_ref)[i]<-substr(paste(HCVFiles_overview2[i]),start=1,stop=7)
}



#############
## Create figures (for individual samples)
# 1) Selection Coefficients across the genome
j=1

for (j in 1:length(HCVFiles)){  
        id<-substr(paste(HCVFiles_overview[j]),start=1,stop=7)
        print(id)
        OverviewDF<-Overview_sum_ref[[j]]

        
        #Make a figure with the selection coefficients across genome
                pdf(paste0("Output/Ref/SelCoef_",id,".pdf"),width=15,height=7.5)
                par(mfrow=c(1,1))
                maxnuc=8500
                par(mar = c(3,5,1,2))
                selcoeffcolumn = which(names(OverviewDF)=="EstSelCoeff")

                plot(OverviewDF$pos[1:maxnuc],OverviewDF[1:maxnuc,selcoeffcolumn],
                        log="y", ylab="Estimated Selection Coefficient (cost)",cex.lab=1.4,
                        yaxt="n", xlab="",xaxt='n',
                        col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-5,1),xlim=c(340,maxnuc))
                axis(1,at=c(seq(500,8500,by=1000)), labels=c(seq(500,8500,by=1000)))
                eaxis(side = 2, at = 10^((-0):(-(5))), cex=2)

                for(i in 1:5){abline(h = 1:10 * 10^(-i), col = "gray41")}

        for (i in 1:maxnuc){
                c=0; co = 1
                if (is.na(OverviewDF$Type[i])==T) next
                if (OverviewDF$Type[i]=="stop"&OverviewDF$ref[i]%in%c("g","c")) {c=1;p=21}
                if (OverviewDF$Type[i]=="syn"&OverviewDF$ref[i]%in%c("g","c")) {c=cols[3];p=21}
                if (OverviewDF$Type[i]=="syn"&OverviewDF$ref[i]%in%c("a","t")) {c=cols[3];p=21}
                if (OverviewDF$Type[i]=="nonsyn"&OverviewDF$ref[i]%in%c("c","g")) {c=cols[4];p=21}
                if (OverviewDF$Type[i]=="nonsyn"&OverviewDF$ref[i]%in%c("a","t")) {c=cols[5];p=21}
                if (c!=0) points(OverviewDF$pos[i],OverviewDF[i,selcoeffcolumn],pch=p,col='gray30',lwd=0.5,
                                bg=rgb(red=col2rgb(c)[1]/255,
                                green=col2rgb(c)[2]/255,
                                blue=col2rgb(c)[3]/255,
                                maxColorValue = 1,alpha=0.8),cex=.6)
                }

        #Add legend
        legpos=7400; legposV=0.53
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
        }


########  2) ############
# 2) use Selection Coefficient estimated using HIV mutation rates

for (j in 1:length(HCVFiles_overview2)){   
                id<-substr(paste(HCVFiles_overview2[j]),start=1,stop=7)
                print(id)
                OverviewDF<-Overview_sum_ref[[j]]
                pdf(paste0("Output/EstSelCoeff_HCV_HIV1_",id,".pdf"),width=15,height=7.5)
                par(mfrow=c(1,1))
                maxnuc=8500
                par(mar = c(3,5,1,2))
                selcoeffcolumn = which(names(OverviewDF)=="EstSelCoeff_hiv")
                
                
                plot(OverviewDF$pos[1:maxnuc],OverviewDF[1:maxnuc,selcoeffcolumn],
                     log="y", ylab="Estimated Selection Coefficient (cost)",cex.lab=1.4,
                     yaxt="n", xlab="",xaxt='n',
                     col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-4,1),xlim=c(340,maxnuc))
                axis(1,at=c(seq(500,8500,by=1000)), labels=c(seq(500,8500,by=1000)))
                eaxis(side = 2, at = 10^((-0):(-(5))), cex=2)
                
                for(i in 1:5){abline(h = 1:10 * 10^(-i), col = "gray60")}
                
                for (i in 1:maxnuc){
                        c=0; co = 1
                        if (is.na(OverviewDF$Type[i])==T) next
                        if (OverviewDF$Type[i]=="stop"&OverviewDF$ref[i]%in%c("g","c")) {c=1;p=21}
                        if (OverviewDF$Type[i]=="syn"&OverviewDF$ref[i]%in%c("g","c")) {c=cols[3];p=21}
                        if (OverviewDF$Type[i]=="syn"&OverviewDF$ref[i]%in%c("a","t")) {c=cols[3];p=21}
                        if (OverviewDF$Type[i]=="nonsyn"&OverviewDF$ref[i]%in%c("c","g")) {c=cols[4];p=21}
                        if (OverviewDF$Type[i]=="nonsyn"&OverviewDF$ref[i]%in%c("a","t")) {c=cols[5];p=21}
                        if (c!=0) points(OverviewDF$pos[i],OverviewDF[i,selcoeffcolumn],pch=p,
                                         col='gray40',lwd=0.5,
                                         bg=rgb(red=col2rgb(c)[1]/255,
                                                green=col2rgb(c)[2]/255,
                                                blue=col2rgb(c)[3]/255,
                                                maxColorValue = 1,alpha=0.8),cex=.7)
                }
 
                #Add legend
                legpos=7400; legposV=0.53
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
        }


##############################################################
##### Test if mean frequenceis differ between mutation types

for (j in 1:length(HCVFiles_overview2)){  
                id<-substr(paste(HCVFiles_overview2[j]),start=1,stop=7)
                print(id)
                OverviewDF<-Overview_sum_ref[[j]]
                
                #Test whether non syn muts, syn muts and nonsense muts are different in freq
                FreqsSyn<-OverviewDF$freq.tranv[OverviewDF$Type.=="syn"]
                FreqsNonSyn<-OverviewDF$freq.maj[OverviewDF$Type=="nonsyn"]
                FreqsStop<-OverviewDF$freq.maj[OverviewDF$Type=="stop"]
                
                print(wilcox.test(FreqsSyn, FreqsNonSyn,alternative = "greater", paired = FALSE))
                print(wilcox.test(FreqsNonSyn,FreqsStop,alternative = "greater", paired = FALSE))
                
                #Test if CpG making mutation frequencies are lower than non-CpGmaking 
                FreqsSyn.CpG<-OverviewDF$freq.maj[OverviewDF$Type=="syn"&OverviewDF$makesCpG==1]
                FreqsSyn.nonCpG<-OverviewDF$freq.maj[OverviewDF$Type=="syn"&OverviewDF$makesCpG==0]
                FreqsNonSyn.CpG<-OverviewDF$freq.maj[OverviewDF$Type=="nonsyn"&OverviewDF$makesCpG==1]
                FreqsNonSyn.nonCpG<-OverviewDF$freq.maj[OverviewDF$Type=="nonsyn"&OverviewDF$makesCpG==0]
                
                print(wilcox.test(FreqsSyn.CpG, FreqsSyn.nonCpG,alternative = "greater", paired = FALSE))
                print(wilcox.test(FreqsNonSyn.CpG, FreqsNonSyn.nonCpG,alternative = "greater", paired = FALSE))                
                

                #Test if CpG making mutation frequencies are lower than non-CpGmaking 
                FreqsSyn.CpG2<-OverviewDF$freq.maj[OverviewDF$Type=="syn"&OverviewDF$makesCpG_all==1]
                FreqsSyn.nonCpG2<-OverviewDF$freq.maj[OverviewDF$Type=="syn"&OverviewDF$makesCpG_all==0]
                FreqsNonSyn.CpG2<-OverviewDF$freq.maj[OverviewDF$Type=="nonsyn"&OverviewDF$makesCpG_all==1]
                FreqsNonSyn.nonCpG2<-OverviewDF$freq.maj[OverviewDF$Type=="nonsyn"&OverviewDF$makesCpG_all==0]
                
                print(wilcox.test(FreqsSyn.CpG2, FreqsSyn.nonCpG2,alternative = "greater", paired = FALSE))
                print(wilcox.test(FreqsNonSyn.CpG2, FreqsNonSyn.nonCpG2, alternative = "greater", paired = FALSE))                
                
                
                
                #Test whether synonymous C to T  muts and G to A , are more costly than A to G 
                CostSynNoCpGA2G<-OverviewDF$EstSelCoeff[OverviewDF$Type=="syn"&OverviewDF$ref=="a"&OverviewDF$makesCpG==0]
                CostSynNoCpGC2T<-OverviewDF$EstSelCoeff[OverviewDF$Type=="syn"&OverviewDF$ref=="c"&OverviewDF$makesCpG==0]
                CostSynNoCpGG2A<-OverviewDF$EstSelCoeff[OverviewDF$Type=="syn"&OverviewDF$ref=="g"&OverviewDF$makesCpG==0]
                print(wilcox.test(CostSynNoCpGC2T,CostSynNoCpGA2G,alternative = "greater", paired = FALSE))
                print(wilcox.test(CostSynNoCpGG2A,CostSynNoCpGA2G,alternative = "greater", paired = FALSE))
                
                #Test whether non-synonymous C to T  muts and G to A , are more costly than A to G 
                CostnonSynNoCpGA2G<-OverviewDF$EstSelCoeff[OverviewDF$Type=="nonsyn"&OverviewDF$ref=="a"&OverviewDF$makesCpG==0]
                CostnonSynNoCpGC2T<-OverviewDF$EstSelCoeff[OverviewDF$Type=="nonsyn"&OverviewDF$ref=="c"&OverviewDF$makesCpG==0]
                CostnonSynNoCpGG2A<-OverviewDF$EstSelCoeff[OverviewDF$Type=="nonsyn"&OverviewDF$ref=="g"&OverviewDF$makesCpG==0]
                print(wilcox.test(CostnonSynNoCpGC2T,CostnonSynNoCpGA2G,alternative = "greater", paired = FALSE))
                print(wilcox.test(CostnonSynNoCpGG2A,CostnonSynNoCpGA2G,alternative = "greater", paired = FALSE))
                
}    
# output saved in a pages file.      



####### MUTATION FREQUENCY  ####################################
######## plot mutation rates in a the way as selection coefficient.

for (j in 1:length(HCVFiles_overview2)){ 
        
        id<-substr(paste(HCVFiles_overview2[j]),start=1,stop=7)
        print(id)
        OverviewDF<-Overview_sum_ref[[j]]
        
        #Make a figure with the selection coefficients across genome
        pdf(paste0("Output/Mut.Rates_HCV_",id,".pdf"),width=15,height=7.5)
        par(mfrow=c(1,1))
        maxnuc=8500
        par(mar = c(3,5,1,2))
        mutcolumn = which(names(OverviewDF)=="freq.maj")
        
        
        plot(OverviewDF$pos[1:maxnuc],OverviewDF[1:maxnuc,mutcolumn],
             log="y", ylab="Observed mutation frequency",cex.lab=1.4,
             yaxt="n", xlab="",xaxt='n',
             col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-5,1),xlim=c(340,maxnuc))
        axis(1,at=c(seq(500,8500,by=1000)), labels=c(seq(500,8500,by=1000)))
        eaxis(side = 2, at = 10^((-0):(-(5))), cex=2)
        
        for(i in 1:5){abline(h = 1:10 * 10^(-i), col = "gray60")}
        
        for (i in 1:maxnuc){
                c=0; co = 1
                if (is.na(OverviewDF$Type[i])==T) next
                if (OverviewDF$Type[i]=="stop"&OverviewDF$ref[i]%in%c("g","c")) {c=1;p=21}
                if (OverviewDF$Type[i]=="syn"&OverviewDF$ref[i]%in%c("g","c")) {c=cols[3];p=21}
                if (OverviewDF$Type[i]=="syn"&OverviewDF$ref[i]%in%c("a","t")) {c=cols[3];p=21}
                if (OverviewDF$Type[i]=="nonsyn"&OverviewDF$ref[i]%in%c("c","g")) {c=cols[4];p=21}
                if (OverviewDF$Type[i]=="nonsyn"&OverviewDF$ref[i]%in%c("a","t")) {c=cols[5];p=21}
                if (c!=0) points(OverviewDF$pos[i],OverviewDF[i,mutcolumn],pch=p,
                                 col='gray30',lwd = .5,
                                 bg=rgb(red=col2rgb(c)[1]/255,
                                        green=col2rgb(c)[2]/255,
                                        blue=col2rgb(c)[3]/255,
                                        maxColorValue = 1,alpha=0.8),cex=.7)
        }
        
        ##Add "Protease" and "RT" words
        #rect(0, 0.000001, 1200, 3.5*10^-4, density = NULL, angle = 45,col=1,border = NA)
        #text(55*3,2.9*10^-4,"PROTEASE",col="white", cex=1.2)
        #rect(297.5, 0.000001, 1200, 3.5*10^-4, density = NULL, angle = 45,col="grey40",border = NA)
        #text(220*3,2.9*10^-4,"REVERSE TRANSCRIPTASE",col="white", cex=1.2)
        
        #Add legend
        legpos=7400; legposV=0.53
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
}


col2rgb(c("#FF0000", "blue"))
           
####### MUTATION FREQUENCY  ####################################
## plot transverion mutation frequencies
#Figure out the colors
rgb(1,1,1,alpha=0.9,max=1)
#"#FF0000CC"

col80<-c("#FFFFFFCC","#228833CC","#CCBB44CC","#EE6677CC","#AA3377CC")
for (j in 1:length(HCVFiles_overview2)){ 
         
        id<-substr(paste(HCVFiles_overview2[j]),start=1,stop=7)
        print(id)
        OverviewDF<-Overview_sum_ref[[j]]
        
        for (i in 1:nrow(OverviewDF)){ 
                if (is.na(OverviewDF$Type.tv1[i])|is.na(OverviewDF$Type.tv2[i])) next
                else {if (OverviewDF$Type.tv1[i] == OverviewDF$Type.tv2[i])  OverviewDF$Type.tvs[i]<-OverviewDF$Type.tv1[i]
                if (OverviewDF$Type.tv1[i] != OverviewDF$Type.tv2[i]) OverviewDF$Type.tvs[i]<-'mixed'}
        }        
        
        
        #Make a figure with the selection coefficients across genome
        pdf(paste0("Output/Transv.Mut.Rates_HCV_",id,".pdf"),width=15,height=7.5)
        par(mfrow=c(1,1))
        par(mar = c(3,5,1,2))
        
        maxnuc=8500
        mutcolumn = which(names(OverviewDF)=="freq.transv")
        
        
        plot(OverviewDF$pos[1:maxnuc],OverviewDF[1:maxnuc,mutcolumn],
             log="y", ylab="Observed mutation frequency",cex.lab=1.4,
             yaxt="n", xlab="",xaxt='n',
             col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-5,1),xlim=c(340,maxnuc))
        axis(1,at=c(seq(500,8500,by=1000)), labels=c(seq(500,8500,by=1000)))
        eaxis(side = 2, at = 10^((-0):(-(5))), cex=2)
        
        for(i in 1:5){abline(h = 1:10 * 10^(-i), col = "gray60")}
        
        for (i in 1:maxnuc){
                if (is.na(OverviewDF$Type.tvs[i])==T) next
                if (OverviewDF$Type.tvs[i]=="stop") {c=col80[1];p=21}
                if (OverviewDF$Type.tvs[i]=="syn") {c=col80[3];p=21}
                if (OverviewDF$Type.tvs[i]=="nonsyn"&OverviewDF$ref[i]%in%c("c","g")) {c=col80[4];p=21}
                if (OverviewDF$Type.tvs[i]=="nonsyn"&OverviewDF$ref[i]%in%c("a","t")) {c=col80[5];p=21}
                if (OverviewDF$Type.tvs[i]=="mixed"){c=col80[2];p=21}
                points(OverviewDF$pos[i],OverviewDF[i,mutcolumn],pch=p,
                                 col='gray30',lwd = .5,
                                 bg=c,cex=.7)
        }
        
        ##Add "Protease" and "RT" words
        #rect(0, 0.000001, 1200, 3.5*10^-4, density = NULL, angle = 45,col=1,border = NA)
        #text(55*3,2.9*10^-4,"PROTEASE",col="white", cex=1.2)
        #rect(297.5, 0.000001, 1200, 3.5*10^-4, density = NULL, angle = 45,col="grey40",border = NA)
        #text(220*3,2.9*10^-4,"REVERSE TRANSCRIPTASE",col="white", cex=1.2)
        
        #Add legend
        legpos=180; legposV=0.53
        rect(legpos, .29*legposV, (legpos+1000), 1.9*legposV, density = NULL, angle = 45,col=alpha("white",1))
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
}











##############################################################

# Not Finished Yet (8/17/18)

#Make a figure with single site frequency spectra for HCV genome
#Currently figure 1 in the revision P Genetics Sept 2017
if (TRUE){
        pdf("Output/SingleSiteFrequencySpectra_HCV_8.2018.pdf",width=7,height=7)
        
        zerobar=50; h2=22; x1=0.25
        #cols <- c(0,brewer.pal(6, "Set2")[c(2, 1)])
        layout(matrix(c(1,2,3, 4, 5, 6, 7,7,7), 3, 3, byrow = TRUE), 
               widths=c(1,1,1), heights=c(1,1,2))
        par(mar = c(1,3,4,2))
        for (i in 172:174){
                #first create empty plot with title
                if (i == 172){
                        t=paste("nonsense mutation",sep="")
                        hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),yaxt="n",
                             col=1,border=0,
                             #    main= bquote(paste(.(t),(C %->% T ))), cex=1.3,
                             main="",cex=1.2,
                             xlab="Frequency", ylab="Count",cex.lab=1.4)
                        title(t,cex=1.2,line=0)
                        text(x1,30,"observed data",cex=1.3)
                        text(x1,h2,"(C172T)",cex=1.3)
                        mtext(text = "A", side = 3, at=-.0, cex=1.7,col=1)
                }
                if (i == 173){
                        t=paste("         non-synonymous mutation",sep="")
                        #    t=paste("Protease: site ", i,"\n non-synonymous mutation",sep="")
                        hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),yaxt="n",
                             col=cols[5],border=0,
                             #    main = bquote(paste(.(t),(A %->% G ))), cex=1.3,
                             main= "", cex=1.2,
                             xlab="Frequency", ylab="Count",cex.lab=1.4)
                        title(t,cex=1.2,line=0)
                        text(x1,30,"observed data",cex=1.3)
                        text(x1,h2,"(A173G)",cex=1.3)
                        
                }
                if (i == 174){
                        t=paste("   synonymous mutation",sep="")
                        #    t=paste("Protease: site ", i,"\n synonymous mutation",sep="")
                        hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),yaxt="n",
                             col=cols[3],border=0,
                             #    main = bquote(paste(.(t),(G %->% A ))), cex=1.3,
                             main= "", cex=1.2,
                             xlab="Frequency", ylab="Count",cex.lab=1.4)
                        title(t,cex=1.2,line=0)
                        text(x1,30,"observed data",cex=1.3)
                        text(x1,h2,"(G174A)",cex=1.3)
                        
                }
                #Next, show  0 bar
                if (i == 172){
                        hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),
                             yaxt="n",col=1,add=T)}
                if (i == 173){
                        hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),
                             yaxt="n",col=cols[5],add=T)}
                if (i == 174){
                        hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),
                             yaxt="n",col=cols[3],add=T)}
                
                #next show all data (unfiltered), but only until 50 for 0 cat
                if (i == 172){
                        hist(c(rep(0,min(zerobar-10,length(which(freqPatTs[1,i]<0.02)))),freqPatTs[1,i][which(freqPatTs[1,i]>0)]),
                             breaks=seq(0,1,by=0.02),add=T,
                             col=1)}
                if (i == 173){
                        hist(c(rep(0,min(zerobar-10,length(which(freqPatTs0[,i]<0.02)))),freqPatTs0[,i][which(freqPatTs0[,i]>0)]),
                             breaks=seq(0,1,by=0.02),add=T,
                             col=cols[5])}
                if (i == 174){
                        hist(c(rep(0,min(zerobar-10,length(which(freqPatTs[,i]<0.02)))),freqPatTs0[,i][which(freqPatTs0[,i]>0)]),
                             breaks=seq(0,1,by=0.02),add=T,
                             col=cols[3])}
                
                axis(2,labels = c(10,20,30,max(zerobar,length(which(freqPatTs0[,i]<0.02)))), 
                     at = c(10,20,30,zerobar), las=1)
                if (length(which(freqPatTs0[,i]<0.02))>=zerobar){
                        axis.break(axis=2,breakpos=zerobar-10,bgcol="white",breakcol="black",style="slash",brw=0.02)
                        points(c(0.01,0.02),c(zerobar-10,zerobar-10),pch=15,cex=2.5,col="white")
                }else{axis(2,labels = zerobar-10,at=zerobar-10,las=1)}
                
        }
        #next, show simulated frequencies
        #    par(fig=c(2/3,1,0,1), new = TRUE)
        par(mar = c(4.5,3,1,2))
        for (i in 172:174){
                if (i ==172)Freqs<-read.csv("Output/SimFreqs172.csv",row.names=1)[1][,1]
                if (i ==173)Freqs<-read.csv("Output/SimFreqs173.csv",row.names=1)[1][,1]
                if (i ==174)Freqs<-read.csv("Output/SimFreqs174.csv",row.names=1)[1][,1]
                t=paste("simulated data",sep="")
                hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),yaxt="n",
                     col=cols[3],border=0,
                     main="",cex=1.2,
                     xlab="Frequency", ylab="Count",cex.lab=1.4)
                #title(t,cex=1.2,line=0)
                text(x1,30,"simulated data",cex=1.3)
                if (i ==172)text(x1,h2,"(s=1)",cex=1.3)
                if (i ==173)text(x1,h2,paste("(s=",round(OverviewDF$EstSelCoeff[173],3),")",sep=""),cex=1.3)
                if (i ==174)text(x1,h2,paste("(s=",round(OverviewDF$EstSelCoeff[174],3),")",sep=""),cex=1.3)
                
                hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),
                     yaxt="n",col=cols[5],add=T)
                hist(c(rep(0,
                           min(zerobar,length(which(Freqs<0.02)))
                ),
                Freqs[which(Freqs>=0.02)]),
                breaks=seq(0,1,by=0.02),add=T,
                col=cols[1])
                axis(2,labels = c(10,20,30,max(zerobar,length(which(Freqs<0.02)))), 
                     at = c(10,20,30,zerobar), las=1)
                if (length(which(Freqs<0.02))>=zerobar){
                        axis.break(axis=2,breakpos=zerobar-10,bgcol="white",breakcol="black",style="slash",brw=0.02)
                        points(c(0.01,0.02),c(zerobar-10,zerobar-10),pch=15,cex=2.5,col="white")
                }else{axis(2,labels = zerobar-10,at=zerobar-10,las=1)}
        }
        
        par(mar = c(1,3,4,4))
        plotter("Bacheler") #this plotter function comes from "ranking_Ordered_Figure1.R"
        mtext(text = "B", side = 3, at=-50, cex=1.7,col=1)
        
        dev.off()
}

