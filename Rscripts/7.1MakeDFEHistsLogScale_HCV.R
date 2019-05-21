
source("Rscripts/baseRscript.R")

# read the files saved in Overview_output:

HCVFiles_overview<-list.files("Output/Overview_output/",pattern="overview.csv")

Overview_summary<-list()
for (i in 1:length(HCVFiles_overview)){ 
        overviews<-read.csv(paste0("Output/Overview_output/",HCVFiles_overview[i]),stringsAsFactors=FALSE)
        overviews<-overviews[,-1]
        Overview_summary[[i]]<-overviews
        names(Overview_summary)[i]<-substr(paste(HCVFiles_overview[i]),start=1,stop=7)
}

##############

breaks=seq(-6,0,by=.2)

for (i in 1:length(Overview_summary)){
        OverviewData<-Overview_summary[[i]]

        filename<-paste0("Output/SelCoeff/DFE-",names(Overview_summary)[i],".pdf")
        pdf(filename, width = 12, height = 7)
        par(mfcol=c(2,4))
        par(mar = c(4,4.3,2.5,1))
        
        for (typeofsite in c("syn", "nonsyn")){
                dat<-OverviewData
                
        k=1
        mycol <- rgb(190, 190, 190, max = 255, alpha = 100, names = "grey50")
        for (wtnt in c("a", "t", "c", "g")){
                datavector<- dat$EstSelCoeff[dat$Type==typeofsite & dat$MajNt==wtnt]
                datavector<-na.omit(datavector)
                
                hist(log10(datavector), breaks = breaks, xlim=c(-6,-0),col=cols[k], 
                main = paste0(wtnt, " : ", typeofsite), ylab="Count", xlab="",
                xaxt="n", cex.main=1.7, las=1, cex.axis =1.2,cex.lab=1.6)
                        
                #datavectorNonCpG<- dat$EstSelCoeff[dat$Type==typeofsite & dat$MajNt==wtnt & dat$makesCpG==0]
                #datavectorNonCpG<-as.numeric(datavector)
                #datavectorNonCpG<-na.omit(datavector)
                #if (k<0) 
                
                #hist(log10(datavectorNonCpG), breaks = breaks, xlim=c(-5,-0),col=mycol, 
                 #       main = paste0(wtnt, " : ", typeofsite), ylab="Count", xlab="",xaxt="n",add=TRUE)
                ##This last line highlights the effect of CpG creating mutations, but I didn't use it. 
                  abline(v=median(log10(datavector), na.rm = TRUE),col="white",lwd=3,lty=1)
                       abline(v=median(log10(datavector), na.rm = TRUE),lwd=3,lty=2)
            
            x=seq(-6,0,by=1)
            axis(1, at=x,labels=paste("10^",x,sep=""), col.axis="black", las=1, line=-.5, cex.axis=1.2)
            mtext(text = "Selection Coefficient", side = 1, line = 1.8, cex=1.1)
            
            k=k+1
        }
    }
        dev.off()

}

################################################

# calculate the average selection coefficient for each mutation type (transition)
#
s<-data.frame()
s_nonCpG<-data.frame()
sel_CpG<-data.frame()
s_nonCpG2<-data.frame()
sel_CpG2<-data.frame()

sel.coef<-list()
sel.coef.CpG<-list()
sel.coef.nonCpG<-list()
sel.coef.CpG2<-list()
sel.coef.nonCpG2<-list()

for (i in 1:length(Overview_summary)){
        dat<-Overview_summary[[i]]
        filename<-names(Overview_summary)[i]
        print(filename)
        for (typeofsite in c("syn", "nonsyn")){
                for (wtnt in c("a", "t", "c", "g")){
                        selcoeff<- dat$EstSelCoeff[dat$Type==typeofsite & dat$MajNt==wtnt]
                        selcoeff<-selcoeff[!is.na(selcoeff)]
                        ave.sel.coef<-mean(selcoeff)
                        s[typeofsite,wtnt]<-ave.sel.coef
                        
                        s_NonCpG<- dat$EstSelCoeff[dat$Type==typeofsite & dat$MajNt==wtnt & dat$makesCpG==0]
                        s_NonCpG<-as.numeric(s_NonCpG)
                        s_NonCpG<-s_NonCpG[!is.na(s_NonCpG)]
                        s_nonCpG[typeofsite,wtnt]<-mean(s_NonCpG)
                        
                        s_CpG<-dat$EstSelCoeff[dat$Type==typeofsite & dat$MajNt==wtnt & dat$makesCpG==1]
                        s_CpG<-as.numeric(s_CpG)
                        s_CpG<-s_CpG[!is.na(s_CpG)]
                        sel_CpG[typeofsite,wtnt]<-mean(s_CpG)
                        
                        #s_NonCpG2<- dat$EstSelCoeff[dat$Type==typeofsite & dat$MajNt==wtnt & dat$makesCpG_all==0]
                        #s_NonCpG2<-as.numeric(s_NonCpG2)
                        #s_NonCpG2<-s_NonCpG2[!is.na(s_NonCpG2)]
                        #s_nonCpG2[typeofsite,wtnt]<-mean(s_NonCpG2)
                        
                        #s_CpG2<-dat$EstSelCoeff[dat$Type==typeofsite & dat$MajNt==wtnt & dat$makesCpG_all==1]
                        #s_CpG2<-as.numeric(s_CpG2)
                        #s_CpG2<-s_CpG2[!is.na(s_CpG2)]
                        #sel_CpG2[typeofsite,wtnt]<-mean(s_CpG2)
                }
        }
        sel.coef[[i]]<-s
        sel.coef.CpG[[i]]<-sel_CpG
        sel.coef.nonCpG[[i]]<-s_nonCpG
        sel.coef.CpG2[[i]]<-sel_CpG2
        sel.coef.nonCpG2[[i]]<-s_nonCpG2
        names(sel.coef)[i]<-filename
        names(sel.coef.CpG)[i]<-filename
        names(sel.coef.nonCpG)[i]<-filename
        names(sel.coef.CpG2)[i]<-filename
        names(sel.coef.nonCpG2)[i]<-filename
        
        #write.csv(s,paste0("Output/SelCoeff/SelCoeffSummary-",names(Overview_summary[i]),".csv"))
}


## Compare Selection Coeffficient of syn vs nonsyn transition mutations
syn<-do.call(rbind, lapply(sel.coef, function(x) x[1,]))
nonsyn<-do.call(rbind, lapply(sel.coef, function(x) x[2,]))
write.csv(syn,"SummaryStats/Sel.coef.syn.summary_8.24.18.csv")
#D75005, D75011 have high sel.coefficient values (3, 6 rows)

colMeans(syn)->syn.ave
syn.ave
#         a          t          c          g 
#0.02841626 0.02406681 0.07002692 0.08010234 

write.csv(nonsyn,"SummaryStats/Sel.coef.nonsyn.summary_8.24.18.csv")
colMeans(nonsyn)->nonsyn.ave
nonsyn.ave
#a          t          c          g 
#0.04220739 0.04588764 0.11487751 0.11108594


CpG<-do.call(rbind, lapply(sel.coef.CpG, function(x) x[1,1:2]))
write.csv(CpG,"SummaryStats/Sel.coef.CpGmaking.summary_8.24.18.csv")
colMeans(CpG)->CpG.ave
CpG.ave
#          a           t 
#0.02878017 0.02593796


nonCpG<-do.call(rbind, lapply(sel.coef.nonCpG, function(x) x[1,]))
write.csv(nonCpG,"SummaryStats/Sel.coef.nonCpGmaking.summary_8.24.18.csv")
colMeans(nonCpG)->nonCpG.ave
nonCpG.ave
#          a           t           c           g 
#0.02827314 0.02296779 0.07002692 0.08010234

### not finished ####
#all CpG makign mutations (not just transition) --create mutation frequency for all mutations
CpG2<-do.call(rbind, lapply(sel.coef.CpG2, function(x) x[1,1:2]))
write.csv(CpG2,"SummaryStats/SummaryStats/Sel.coef.CpGmaking_all.summary_8.24.18.csv")
colMeans(CpG2)
#          a           t 
#0.03133402 0.02923428

nonCpG2<-do.call(rbind, lapply(sel.coef.nonCpG2, function(x) x[1,]))
write.csv(nonCpG,"SummaryStats/SummaryStats/Sel.coef.nonCpGmaking_all.summary_8.24.18.csv")
colMeans(nonCpG2)
#a           t           c           g 
#0.03491145 0.02716698 0.07950244 0.08745244 
#####


####################################################
# calculate the average selection coefficient for each mutation type (transversion)
# Transversion 1 (to A or C)
s<-data.frame()
s_nonCpG<-data.frame()
sel_CpG<-data.frame()
sel.coef<-list()
sel.coef.CpG<-list()
sel.coef.nonCpG<-list()

for (i in 1:length(Overview_summary)){
        dat<-Overview_summary[[i]]
        filename<-names(Overview_summary)[i]
        for (typeofsite in c("syn", "nonsyn")){
                for (wtnt in c("a", "t", "c", "g")){
                        selcoeff<- dat$EstSelCoeff_T1[dat$Type.tv1==typeofsite & dat$MajNt==wtnt]
                        selcoeff<-selcoeff[!is.na(selcoeff)]
                        ave.sel.coef<-mean(selcoeff)
                        s[typeofsite,wtnt]<-ave.sel.coef
                        
                        s_NonCpG<- dat$EstSelCoeff_T1[dat$Type.tv1==typeofsite & dat$MajNt==wtnt & dat$makesCpG.tv1==0]
                        s_NonCpG<-as.numeric(s_NonCpG)
                        s_NonCpG<-s_NonCpG[!is.na(s_NonCpG)]
                        s_nonCpG[typeofsite,wtnt]<-mean(s_NonCpG)
                        
                        s_CpG<-dat$EstSelCoeff_T1[dat$Type.tv1==typeofsite & dat$MajNt==wtnt & dat$makesCpG.tv1==1]
                        s_CpG<-as.numeric(s_CpG)
                        s_CpG<-s_CpG[!is.na(s_CpG)]
                        sel_CpG[typeofsite,wtnt]<-mean(s_CpG)
                }
        }
        sel.coef[[i]]<-s
        sel.coef.CpG[[i]]<-sel_CpG
        sel.coef.nonCpG[[i]]<-s_nonCpG
        names(sel.coef)[i]<-filename
        names(sel.coef.CpG)[i]<-filename
        names(sel.coef.nonCpG)[i]<-filename
        #write.csv(s,paste0("Output/SelCoeffSummary_T1-",names(Overview_summary[i]),".csv"))
}

## Compare Selection Coeffficient of syn vs nonsyn transversion mutations
syn1<-do.call(rbind, lapply(sel.coef, function(x) x[1,]))
nonsyn1<-do.call(rbind, lapply(sel.coef, function(x) x[2,]))
write.csv(syn1,"SummaryStats/Sel.coef.syn.tv1.summary_8.24.18.csv")
#D75005, D75011 have high sel.coefficient values (3, 6 rows)
colMeans(syn1)
colMeans(syn1)->syn1.ave
#         a          t          c          g 
#0.4609333 0.2300671 0.2070777 0.6536702 

write.csv(nonsyn1,"SummaryStats/Sel.coef.nonsyn.tv1.summary_8.24.18.csv")
colMeans(nonsyn1)->nonsyn1.ave
nonsyn1.ave
#a          t          c          g 
#0.5000507 0.3201027 0.3553874 0.7359491 


CpG1<-do.call(rbind, lapply(sel.coef.CpG, function(x) x[1,]))
write.csv(CpG1,"SummaryStats/Sel.coef.CpGmaking.tv1.summary_8.24.18.csv")
colMeans(CpG1)->CpG1.ave
#        a         t         c         g 
#      NaN 0.2021844 0.1826322       NaN 

nonCpG1<-do.call(rbind, lapply(sel.coef.nonCpG, function(x) x[1,]))
write.csv(nonCpG1,"SummaryStats/Sel.coef.nonCpGmaking.tv1.summary_8.24.18.csv")
colMeans(nonCpG1)->nonCpG1.ave
#        a         t         c         g 
#0.4609333 0.2604837 0.2290468 0.6536702


####################################################
# calculate the average selection coefficient for each mutation type (transversion)
# Transversion 2 (to T or G)
s<-data.frame()
s_nonCpG<-data.frame()
sel_CpG<-data.frame()
sel.coef<-list()
sel.coef.CpG<-list()
sel.coef.nonCpG<-list()

for (i in 1:length(Overview_summary)){
        dat<-Overview_summary[[i]]
        filename<-names(Overview_summary)[i]
        for (typeofsite in c("syn", "nonsyn")){
                for (wtnt in c("a", "t", "c", "g")){
                        selcoeff<- dat$EstSelCoeff_T2[dat$Type.tv2==typeofsite & dat$MajNt==wtnt]
                        selcoeff<-selcoeff[!is.na(selcoeff)]
                        ave.sel.coef<-mean(selcoeff)
                        s[typeofsite,wtnt]<-ave.sel.coef
                        
                        
                        s_NonCpG<- dat$EstSelCoeff_T2[dat$Type.tv2==typeofsite & dat$MajNt==wtnt & dat$makesCpG.tv2==0]
                        s_NonCpG<-as.numeric(s_NonCpG)
                        s_NonCpG<-s_NonCpG[!is.na(s_NonCpG)]
                        s_nonCpG[typeofsite,wtnt]<-mean(s_NonCpG)
                        
                        s_CpG<-dat$EstSelCoeff_T2[dat$Type.tv2==typeofsite & dat$MajNt==wtnt & dat$makesCpG.tv2==1]
                        s_CpG<-as.numeric(s_CpG)
                        s_CpG<-s_CpG[!is.na(s_CpG)]
                        sel_CpG[typeofsite,wtnt]<-mean(s_CpG)
                }
        }
        sel.coef[[i]]<-s
        sel.coef.CpG[[i]]<-sel_CpG
        sel.coef.nonCpG[[i]]<-s_nonCpG
        names(sel.coef)[i]<-filename
        names(sel.coef.CpG)[i]<-filename
        names(sel.coef.nonCpG)[i]<-filename
        #write.csv(s,paste0("Output/SelCoeffSummary_T1-",names(Overview_summary[i]),".csv"))
}


## Compare Selection Coeffficient of syn vs nonsyn transversion mutations
syn2<-do.call(rbind, lapply(sel.coef, function(x) x[1,]))
nonsyn2<-do.call(rbind, lapply(sel.coef, function(x) x[2,]))
write.csv(syn2,"SummaryStats/Sel.coef.syn.tv2.summary_8.24.18.csv")
#D75005, D75011 have high sel.coefficient values (3, 6 rows)

colMeans(syn2)->syn2.ave
syn2.ave
#         a          t          c          g 
#0.2777536 0.4280874 0.6924311 0.2799253

write.csv(nonsyn2,"SummaryStats/Sel.coef.nonsyn.tv2.summary_8.24.18.csv")
colMeans(nonsyn2)->nonsyn2.ave
nonsyn2.ave
#a          t          c          g 
#0.3396992 0.4805756 0.7382025 0.3308271



CpG2<-do.call(rbind, lapply(sel.coef.CpG, function(x) x[1,]))
write.csv(CpG2,"SummaryStats/Sel.coef.CpGmaking.tv2.summary_8.24.18.csv")
colMeans(CpG2)->CpG2.ave
CpG2.ave
#        a         t         c         g 
#0.2072462       NaN       NaN 0.1827888 


nonCpG2<-do.call(rbind, lapply(sel.coef.nonCpG, function(x) x[1,]))
write.csv(nonCpG2,"SummaryStats/Sel.coef.nonCpGmaking.tv2.summar_8.24.18y.csv")
colMeans(nonCpG2)->nonCpG2.ave
nonCpG2.ave
#         a        t          c         g 
#0.3214179 0.4280874 0.6924311 0.3469681



##################
#calculate means for tranversion mutation frequencies
mean(as.numeric(nonsyn1.ave,nonsyn2.ave)) #0.4778725
mean(as.numeric(syn1.ave,syn2.ave)) # 0.3879371

#Average CpGmaking mutation frequency
as.numeric(CpG1.ave[2:3])
mean(c(as.numeric(CpG1.ave[2:3]),as.numeric(CpG2.ave[c(1,4)]))) #0.1937129
mean(as.numeric(nonCpG1.ave,nonCpG2.ave)) #[1] 0.4010335


types<-paste0(rep(c("syn", "nonsyn","CpG","nonCpG"), times = 3), rep(c('', 1, 2), each = 4),rep('.ave', each=12))

SelCoeff_average<- data.frame(sapply(types,get))
colnames(SelCoeff_average)<-c("Transition_syn","Transition_nonsyn","CpG-making","nonCpG-making","Transvs1_syn","Transvs1_nonsyn","Tvs1_CpG-making","Tvs1_nonCpG-making","Transvs2_syn","Transvs2_nonsyn","Tvs2_CpG-making","Tvs2_nonCpG-making")
write.csv(SelCoeff_average,"Output/SelCoeff_averages_by_mutationTypes.csv")
SelCoeff_average<-t(data.frame(syn.ave))
SelCoeff_average<-rbind(SelCoeff_average,nonsyn.ave)
SelCoeff_average<-combine(SelCoeff_average, CpG.ave,nonCpG.ave)


##############
## Overlay the non-CpG making mutation frequency distribution -> base a and t shows the difference: shaded area = CpG making frequencies

for (i in 1:length(Overview_summary)){
        OverviewData<-Overview_summary[[i]]
        filename<-paste0("Output/SelCoeff/DFE_-NonCpGmaking",names(Overview_summary)[i],".pdf", sep = "")
        pdf(filename, width = 12, height = 7)
        par(mfcol=c(2,4))
        par(mar = c(4,4.3,2.5,1))
        
        for (typeofsite in c("syn", "nonsyn")){
                dat<-OverviewData
                k=1
                mycol <- rgb(190, 190, 190, max = 255, alpha = 100, names = "grey50")
                for (wtnt in c("a", "t", "c", "g")){
                        datavector<- dat$EstSelCoeff[dat$Type==typeofsite & dat$MajNt==wtnt]
                        datavector<-na.omit(datavector)
                        
                        hist(log10(datavector), breaks = breaks, xlim=c(-6,-0),col=cols[k], 
                             main = paste0(wtnt, " : ", typeofsite), ylab="Count", xlab="",
                             xaxt="n", cex.main=1.7, las=1, cex.axis =1.2,cex.lab=1.6,boarder='gray')
                        
                        datavectorNonCpG<- dat$EstSelCoeff[dat$Type==typeofsite & dat$MajNt==wtnt & dat$makesCpG==0]
                        datavectorNonCpG<-na.omit(datavectorNonCpG)
        
                        
                        hist(log10(datavectorNonCpG), breaks = breaks, xlim=c(-6,-0),col='black', angle=45,density=20,
                               main = paste0(wtnt, " : ", typeofsite), ylab="Count", xlab="",xaxt="n",add=T)
                        ##This last line highlights the effect of CpG creating mutations, but I didn't use it. 
                        abline(v=median(log10(datavector), na.rm = TRUE),col="white",lwd=2,lty=1)
                        abline(v=median(log10(datavector), na.rm = TRUE),lwd=2,lty=2)
                        
                        x=seq(-6,0,by=1)
                        axis(1, at=x,labels=paste("10^",x,sep=""), col.axis="black", las=1, line=-.5, cex.axis=1.2)
                        mtext(text = "Selection Coefficient", side = 1, line = 1.8, cex=1.1)
                        
                        k=k+1
                }
        }
        dev.off()
        
}


## Plot CpGmaking and non-CpGmaking together
for (i in 1:length(Overview_summary)){
        OverviewData<-Overview_summary[[i]]
        filename<-paste0("Output/SelCoeff/DFE_CpG.vs.nonCpG",names(Overview_summary)[i],".pdf", sep = "")
        pdf(filename, width = 12, height = 7)
        par(mfcol=c(2,4))
        par(mar = c(4,4.3,2.5,1))
        
        for (typeofsite in c("syn", "nonsyn")){
                dat<-OverviewData
                #colnames(dat)[3]<-"TypeOfSite"
                
                
                mycol <- rgb(190, 190, 190, max = 255, alpha = 100, names = "grey50")
                for (wtnt in c("a", "t")){
                        datavectorNonCpG<- dat$EstSelCoeff[dat$Type==typeofsite & dat$MajNt==wtnt & dat$makesCpG==0]
                        datavectorNonCpG<-na.omit(datavectorNonCpG)
                        
                        hist(log10(datavectorNonCpG), breaks = breaks, xlim=c(-6,-0),col=c("pink"),
                             main = paste0(wtnt, " : ", typeofsite,"-nonCpG"), ylab="Count", xlab="",xaxt="n")
                        ##This last line highlights the effect of CpG creating mutations, but I didn't use it. 
                        abline(v=median(log10(datavectorNonCpG), na.rm = TRUE),col="white",lwd=2,lty=1)
                        abline(v=median(log10(datavectorNonCpG), na.rm = TRUE),lwd=2,lty=2)
                        
                        x=seq(-6,0,by=1)
                        axis(1, at=x,labels=paste("10^",x,sep=""), col.axis="black", las=1, line=-.5, cex.axis=1.2)
                        mtext(text = "Selection Coefficient", side = 1, line = 1.8, cex=1.1)
                        
                        
                        datavectorCpG<- dat$EstSelCoeff[dat$Type==typeofsite & dat$MajNt==wtnt & dat$makesCpG==1]
                        datavectorCpG<-na.omit(datavectorCpG)
                        
                        
                        hist(log10(datavectorCpG), breaks = breaks, xlim=c(-6,-0),col='gray',
                             main = paste0(wtnt, " : ", typeofsite,'-CpG'), ylab="Count", xlab="",xaxt="n")
                        ##This last line highlights the effect of CpG creating mutations, but I didn't use it. 
                        abline(v=median(log10(datavectorCpG), na.rm = TRUE),col="white",lwd=2,lty=1)
                        abline(v=median(log10(datavectorCpG), na.rm = TRUE),lwd=2,lty=2)
                        
                        x=seq(-6,0,by=1)
                        axis(1, at=x,labels=paste("10^",x,sep=""), col.axis="black", las=1, line=-.5, cex.axis=1.2)
                        mtext(text = "Selection Coefficient", side = 1, line = 1.8, cex=1.1)
                        
                        
                }
        }
        dev.off()
        
}
