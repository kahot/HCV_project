library(ggplot2)
library(reshape)
library(ggpubr)
library(ggthemes)
library(plotrix)
library(grid)

source("Rscripts/baseRscript.R")

# read the files saved in Overview_output:
HCVFiles_overview<-list.files("Output1A/Overview1/",pattern="overview.csv")

Overview_summary<-list()
for (i in 1:length(HCVFiles_overview)){ 
        overviews<-read.csv(paste0("Output1A/Overview1/",HCVFiles_overview[i]),stringsAsFactors=FALSE)
        overviews<-overviews[,-1]
        Overview_summary[[i]]<-overviews
        names(Overview_summary)[i]<-substr(paste(HCVFiles_overview[i]),start=1,stop=7)
}


#Data Prep: Calculate the average mutation frequency for each mutation type
#Run the following scripts to get the summary frequency data:
source("Rscripts/MutationFreqSum.R")

#the results are saved in Output1A/MutFreq/ directory
################################################

# Data Prep: Calculate the average mutation frequency for each mutation type

#Plot summary of 1) Transistion mutations
TransFiles<-list.files("Output1A/MutFreq/Maj/",pattern="Transition.csv")
TransFiles<-TransFiles[-6]

NofSamples<-length(Overview_summary)

for (i in 1:length(TransFiles)){
        mdata<-read.csv(paste0("Output1A/MutFreq/Maj/",TransFiles[i]),stringsAsFactors=FALSE,row.names=1)
        mdata<-mdata[1:NofSamples,]
        colnames(mdata)<-c("A","T","C","G")
        Dat<-melt(mdata)
        Dat<-Dat[!(is.na(Dat$value)),]
        
        filename<-sub(".csv$","",paste(TransFiles[i]))

        MFplot<-ggplot(Dat,aes(x=variable,y=value,fill=variable))+geom_boxplot(alpha=0.5)+labs(x="Nucleotide",y="Mutation frequency")+
                scale_fill_manual(values=c("#66CCEECC","#228833CC","#CCBB44CC","#EE6677CC")) + theme_classic()+
                guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
                theme(axis.text.y = element_text(size =10))+
                ggtitle(paste("Ave. freq of ",filename," mutation")) + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
                theme(plot.title = element_text(hjust = 0.5))
        if (i==1|i==2){
                ggsave(filename=paste0("Output1A/MutFreq/Maj/",filename,".pdf"),width=5, height=7, units='in',device='pdf', plot=MFplot)}
        else{
                ggsave(filename=paste0("Output1A/MutFreq/Maj/",filename,".pdf"), width=9, height=7, units='in',device='pdf',plot=MFplot)
        }
        
}

########### Same Transition Mutation Figures but two plots side by side #########
#plots Syn vs NonSyn each other:
for (i in 5:6){ 
        mdata<-read.csv(paste0("Output1A/MutFreq/Maj/",TransFiles[i]),stringsAsFactors=FALSE,row.names=1)
        mdata<-mdata[1:NofSamples,]
        colnames(mdata)<-c("A","T","C","G")
        Dat<-melt(mdata)
        Dat<-Dat[!(is.na(Dat$value)),]
 
        filename<-sub(".csv$","",paste(TransFiles[i]))
        
        MFplot<-ggplot(Dat,aes(x=variable,y=value,fill=variable))+geom_boxplot(alpha=0.5)+labs(x="Nucleotide",y="Mutation frequency")+
                scale_fill_manual(values=c("#66CCEECC","#228833CC","#CCBB44CC","#EE6677CC")) + theme_classic()+
                guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
                theme(axis.text.y = element_text(size =10))+
                ggtitle(paste(filename)) + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
                theme(plot.title = element_text(hjust = 0.5))
        plotname<-paste0("plot_",i)

        assign(plotname,MFplot)
        }
require(gridExtra)
pdf("Output1A/MutFreq/Maj/SynvsNonsynTs.pdf",width=10, height=7)
grid.arrange(plot_5, plot_6, ncol=2)
dev.off()

#### #compare CpG creating vs. non-CpGcreating separately for syn & nonsyn

for (i in 1:2){
        if (i==1) type<-"nonsynonymous"
        if (i==2) type<-"synonymous"
        #read CpG-syn data
        mdata1<-read.csv(paste0("Output1A/MutFreq/Maj/",TransFiles[i]),stringsAsFactors=FALSE,row.names=1)
        mdata1<-mdata1[1:NofSamples,1:2] #remove SE values
        colnames(mdata1)<-c("A(CpG)","T(CpG)")
        #read nonCpg-syn data
        mdata2<-read.csv(paste0("Output1A/MutFreq/Maj/",TransFiles[i+2]),stringsAsFactors=FALSE,row.names=1)
        mdata2<-mdata2[1:NofSamples,] #remove SE values
        colnames(mdata2)<-c("A","T","C","G")

        trans<-cbind(mdata1,mdata2)
        transm<-melt(trans)
        transm<-transm[!(is.na(transm$value)),]

        MFplot1<-ggplot(transm,aes(x=variable,y=value,fill=variable))+geom_boxplot(alpha=0.5)+labs(x="Nucleotide",y="Mutation frequency")+
                scale_fill_manual(values=c("#66CCEE","#228833","#66CCEE99","#22883399","#CCBB44CC","#EE6677CC")) + theme_classic()+
                guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
                theme(axis.text.y = element_text(size =10))+
                ggtitle(paste0("CpG vs. nonCpG ", type, " mutations")) + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
                theme(plot.title = element_text(hjust = 0.5))+
                geom_vline(xintercept = 2.5)

        ggsave(filename=paste0("Output1A/MutFreq/Maj/CpGvsNonCpG.",type,".Ts.pdf"),width=6.3, height=5, units='in',device='pdf', plot=MFplot1)
}


#################
#Plot summary of 2) Transversion mutations         
Tv1Files<-list.files("Output1A/MutFreq/Maj/",pattern="Transversion_1")
Tv2Files<-list.files("Output1A/MutFreq/Maj/",pattern="Transversion_2")

#CpGmaking Transversion plots
for (i in 1:2){
        mdata1<-read.csv(paste0("Output1A/MutFreq/Maj/",Tv1Files[i]),stringsAsFactors=FALSE,row.names=1)
        mdata2<-read.csv(paste0("Output1A/MutFreq/Maj/",Tv2Files[i]),stringsAsFactors=FALSE,row.names=1)
        mdata1<-mdata1[1:NofSamples,]
        mdata2<-mdata2[1:NofSamples,]
        mdata<-cbind(mdata1[,c(1,4)],mdata2[,c(2,3)])
        mdata<-mdata[,c(1,3,4,2)]
        colnames(mdata)<-c("A","T","C","G")
        dat<-melt(mdata)
        
        filename<-sub("_1.csv$","",paste(Tv1Files[i]))

        MFplot<-ggplot(dat,aes(x=variable,y=value,fill=variable))+geom_boxplot(alpha=0.5)+labs(x="Nucleotide",y="Mutation frequency")+
                scale_fill_manual(values=c("#66CCEECC","#228833CC","#CCBB44CC","#EE6677CC")) + theme_classic()+
                guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
                theme(axis.text.y = element_text(size =10))+
                ggtitle(paste0("Ave. freq of ",filename," mutation")) + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
                theme(plot.title = element_text(hjust = 0.5))
        ggsave(filename=paste0("Output1A/MutFreq/Maj/",filename,".pdf"),width=9, height=7, units='in',device='pdf',plot=MFplot)
        
        #dev.off()
}

#non-CpGmaking Transversion plots
for (i in c(3,4,5,7)){
        mdata1<-read.csv(paste0("Output1A/MutFreq/Maj/",Tv1Files[i]),stringsAsFactors=FALSE,row.names=1)
        mdata2<-read.csv(paste0("Output1A/MutFreq/Maj/",Tv2Files[i]),stringsAsFactors=FALSE,row.names=1)
        mdata1<-mdata1[1:NofSamples,]
        mdata2<-mdata2[1:NofSamples,]
        mdata<-mdata1+mdata2
        colnames(mdata)<-c("A","T","C","G")
        dat<-melt(mdata)
        
        filename<-sub("_1.csv$","",paste(Tv1Files[i]))

        MFplot<-ggplot(dat,aes(x=variable,y=value,fill=variable))+geom_boxplot(alpha=0.5)+labs(x="Nucleotide",y="Mutation frequency")+
                scale_fill_manual(values=c("#66CCEECC","#228833CC","#CCBB44CC","#EE6677CC")) + theme_classic()+
                guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
                theme(axis.text.y = element_text(size =10))+
                ggtitle(paste0("Ave. freq of ",filename," mutation")) + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
                theme(plot.title = element_text(hjust = 0.5))
        ggsave(filename=paste0("Output1A/MutFreq/Maj/",filename,".pdf"),width=9, height=7, units='in',plot=MFplot)
        
        
}        

####################
### #compare CpG creating vs. non-CpGcreating Plot A vs A, G vs G etc.


for (i in 1:2){
        if (i==1) type<-"nonsynonymous"
        if (i==2) type<-"synonymous"

        #CpG-syn vs. nonCpG-syn transversion data
        mdata1<-read.csv(paste0("Output1A/MutFreq/Maj/",Tv1Files[i]),stringsAsFactors=FALSE,row.names=1)
        mdata2<-read.csv(paste0("Output1A/MutFreq/Maj/",Tv2Files[i]),stringsAsFactors=FALSE,row.names=1)
        mdata1<-mdata1[1:NofSamples,]
        mdata2<-mdata2[1:NofSamples,]
        cpg1<-cbind(mdata1[,c(1,4)],mdata2[,c(2,3)])
        cpg1<-cpg1[,c(1,3,4,2)]
        colnames(cpg1)<-c("A(CpG)","T(CpG)","C(CpG)","G(CpG)")

        mdata2.2<-read.csv(paste0("Output1A/MutFreq/Maj/",Tv2Files[i+2]),stringsAsFactors=FALSE,row.names=1)
        mdata1.2<-read.csv(paste0("Output1A/MutFreq/Maj/",Tv1Files[i+2]),stringsAsFactors=FALSE,row.names=1)
        mdata1.2<-mdata1.2[1:NofSamples,]
        mdata2.2<-mdata2.2[1:NofSamples,]
        noncpg1<-cbind(mdata1.2[,c(1,4)],mdata2.2[,c(2,3)])
        noncpg1<-noncpg1[,c(1,3,4,2)]
        colnames(noncpg1)<-c("A","T","C","G")

        cpgTv<-cbind(cpg1, noncpg1)
        dat1<-melt(cpgTv)

        MFplot3<-ggplot(dat1,aes(x=variable,y=value,fill=variable))+geom_boxplot(alpha=0.5)+labs(x="Nucleotide",y="Mutation frequency")+
                scale_fill_manual(values=c("#66CCEE","#228833","#CCBB44CC","#EE6677CC","#66CCEE99","#22883399","#CCBB44CC","#EE6677CC"))+
                theme_classic()+
                guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
                theme(axis.text.y = element_text(size =10))+
                ggtitle(paste0("CpG vs. nonCpG ",type, " transversion")) + 
                theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
                theme(plot.title = element_text(hjust = 0.5))+
                geom_vline(xintercept = 4.5)

        ggsave(filename=paste0("Output1A/MutFreq/Maj/CpGvsNonCpG.",type, ".Transv.pdf"),width=6.3, height=5, units='in',device='pdf', plot=MFplot3)
}


#####
#Average all bases and plot together
#TransFiles<-list.files("Output1A/MutFreq/Maj/",pattern="Transition.csv")
#TransFiles<-TransFiles[-6]

datalist<-list()
for (i in 1:length(TransFiles)){
        if (i==1|i==2|i==3|i==4){
                mdata<-read.csv(paste0("Output1A/MutFreq/Maj/",TransFiles[i]),stringsAsFactors=FALSE,row.names=1)
                mdata<-mdata[1:NofSamples,1:2]  #use A and T only for nonCpGmaking mutations
                dat<-melt(mdata)
                dat<-dat[!(is.na(dat$value)),]
                
                if (i==1) dat$variable<-c(rep("NonSyn_CpG", n=length(dat$variable)))
                if (i==2) dat$variable<-c(rep("Syn_CpG", n=length(dat$variable)))
                if (i==3) dat$variable<-c(rep("NonSyn_nonCpG", n=length(dat$variable)))
                if (i==4) dat$variable<-c(rep("Syn_nonCpG", n=length(dat$variable)))
                }
        if (i==5|i==6){
                mdata<-read.csv(paste0("Output1A/MutFreq/Maj/",TransFiles[i]),stringsAsFactors=FALSE,row.names=1)
                mdata<-mdata[1:NofSamples,]
                dat<-melt(mdata)
                dat<-dat[!(is.na(dat$value)),]
                if (i==5) dat$variable<-c(rep("NonSyn", n=length(dat$variable)))
                if (i==6) dat$variable<-c(rep("Syn", n=length(dat$variable)))
                }
        datalist[[i]]<-dat
        }

Transitions<-do.call(rbind,datalist)        
title<-"Transition Mutations"
MFplot<-ggplot(Transitions,aes(x=variable,y=value,fill=variable))+geom_boxplot(alpha=0.5)+labs(x="Mutation Type",y="Mutation frequency")+
        scale_fill_manual(values=c("#66CCEECC","#228833CC","#CCBB44CC","#EE6677CC","#AA3377CC","#4477AACC")) + theme_classic()+
        guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10))+
        ggtitle(paste0(title)) + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
        theme(plot.title = element_text(hjust = 0.5))
ggsave(filename=paste0("Output1A/MutFreq/Maj/",title,"_Summary(A,T only).pdf"),width=9, height=7, units='in',plot=MFplot)




####Transversion

#Tv1Files<-list.files("Output/MutFreq/maj/",pattern="Transversion_1.csv")
#Tv2Files<-list.files("Output/MutFreq/maj/",pattern="Transversion_2.csv")
Tv1Files<-Tv1Files[-6]
Tv2Files<-Tv2Files[-6]

datalist.tv<-list()
for (i in 1:length(Tv1Files)){
        if (i==1|i==2){
                mdata1<-read.csv(paste0("Output1A/MutFreq/Maj/",Tv1Files[i]),stringsAsFactors=FALSE,row.names=1)
                mdata2<-read.csv(paste0("Output1A/MutFreq/Maj/",Tv2Files[i]),stringsAsFactors=FALSE,row.names=1)
                mdata1<-mdata1[1:NofSamples,]
                mdata2<-mdata2[1:NofSamples,]
                mdata<-cbind(mdata1[,c(1,4)],mdata2[,c(2,3)])
                dat<-melt(mdata)
                if (i==1) dat$variable<-c(rep("NonSyn_CpG", n=length(dat$variable)))
                if (i==2) dat$variable<-c(rep("Syn_CpG", n=length(dat$variable)))
        }
        
        else {
                mdata1<-read.csv(paste0("Output1A/MutFreq/Maj/",Tv1Files[i]),stringsAsFactors=FALSE,row.names=1)
                mdata2<-read.csv(paste0("Output1A/MutFreq/Maj/",Tv2Files[i]),stringsAsFactors=FALSE,row.names=1)
                mdata1<-mdata1[1:NofSamples,]
                mdata2<-mdata2[1:NofSamples,]
                mdata<-mdata1+mdata2
                dat<-melt(mdata)
                if (i==3) dat$variable<-c(rep("NonSyn_nonCpG", n=length(dat$variable)))
                if (i==4) dat$variable<-c(rep("Syn_nonCpG", n=length(dat$variable)))
                if (i==5) dat$variable<-c(rep("NonSyn", n=length(dat$variable)))
                if (i==6) dat$variable<-c(rep("Syn", n=length(dat$variable)))
        }
        datalist.tv[[i]]<-dat
}

Transvs<-do.call(rbind,datalist.tv)        
title<-"Transversion Mutations"
colors<-rep(c("#66CCEECC","#228833CC","#CCBB44CC","#EE6677CC"), times=2)
MFplot2<-ggplot(Transvs,aes(x=variable,y=value,fill=variable))+geom_boxplot(alpha=0.5)+labs(x="Mutation Type",y="Mutation frequency")+
        scale_fill_manual(values=c("#66CCEECC","#228833CC","#CCBB44CC","#EE6677CC","#AA3377CC","#4477AACC")) + theme_classic()+
        guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10))+
        ggtitle(paste0(title)) + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
        theme(plot.title = element_text(hjust = 0.5))
ggsave(filename=paste0("Output1A/MutFreq/Maj/",title,"_Summary.pdf"),width=9, height=7, units='in',plot=MFplot2)


##################################
# 1) Histogram of all transition mutation frequencies

#Plot with Median
datanames<-c("A_syn","T_syn","C_syn","G_syn","A_nonsyn","T_nonsyn","C_nonsyn","G_nonsyn")
colors<-c("#66CCEECC", "#228833CC" ,"#CCBB44CC", "#EE6677CC","#66CCEECC", "#228833CC", "#CCBB44CC", "#EE6677CC")
filename<-paste0("Output1A/MutFreq/Maj/Transition_Mut_Freq_Histograms_median.pdf")
pdf(filename, width = 12, height = 7)
par(mfrow=c(2,4))
par(mar = c(4,4.3,2.5,1))
breaks=seq(-6,0,by=.2)

for (i in 1:8){
        datavector<-get(datanames[i])
        nt<-substr(paste(datanames[i]),start=1,stop=1)
        typ<-substr(paste(datanames[i]),start=3, stop=8)
        hist(log10(datavector), breaks = breaks, xlim=c(-5,-0),col=colors[i], 
             main = paste0(nt, " : ", typ), ylab="Count", xlab="",
             xaxt="n", cex.main=1.7, las=1, cex.axis =1,cex.lab=1.4)
        abline(v=median(log10(datavector), na.rm = TRUE),col="white",lwd=2,lty=1)
        abline(v=median(log10(datavector), na.rm = TRUE),col="red",lwd=2,lty=2)
        x=seq(-5,0,by=1)
        labels <- sapply(x,function(y)as.expression(bquote(10^ .(y))))
        axis(1, at=x,labels=labels, col.axis="black", line=-.7, cex.axis=1)
        mtext(text = "Transition Mutation Frequency", side = 1, line = 1.5, cex=.9)
}

dev.off()

#Same plot with Mean
datanames<-c("A_syn","T_syn","C_syn","G_syn","A_nonsyn","T_nonsyn","C_nonsyn","G_nonsyn")
colors<-c("#66CCEECC", "#228833CC" ,"#CCBB44CC", "#EE6677CC","#66CCEECC", "#228833CC", "#CCBB44CC", "#EE6677CC")
filename<-paste0("Output1A/MutFreq/Maj/Transition_Mut_Freq_Histograms_mean.pdf")
pdf(filename, width = 12, height = 7)
par(mfrow=c(2,4))
par(mar = c(4,4.3,2.5,1))
breaks=seq(-6,0,by=.2)

for (i in 1:8){
        datavector<-get(datanames[i])
        nt<-substr(paste(datanames[i]),start=1,stop=1)
        typ<-substr(paste(datanames[i]),start=3, stop=8)
        hist(log10(datavector), breaks = breaks, xlim=c(-5,-0),col=colors[i], 
             main = paste0(nt, " : ", typ), ylab="Count", xlab="",
             xaxt="n", cex.main=1.7, las=1, cex.axis =1,cex.lab=1.4)
        abline(v=log10(mean(datavector, na.rm = TRUE)),col="white",lwd=2,lty=1)
        abline(v=log10(mean(datavector, na.rm = TRUE)),lwd=2,lty=2)
        x=seq(-5,0,by=1)
        labels <- sapply(x,function(y)as.expression(bquote(10^ .(y))))
        axis(1, at=x,labels=labels, col.axis="black", line=-.7, cex.axis=1)
        mtext(text = "Transition Mutation Frequency", side = 1, line = 1.5, cex=.9)
}

dev.off()

##### 2) Histogram of tranversion mutation frequencies
datanames1<-c("A_tv1_syn","T_tv1_syn","C_tv1_syn","G_tv1_syn","A_tv1_nonsyn","T_tv1_nonsyn","C_tv1_nonsyn","G_tv1_nonsyn")
datanames2<-c("A_tv2_syn","T_tv2_syn","C_tv2_syn","G_tv2_syn","A_tv2_nonsyn","T_tv2_nonsyn","C_tv2_nonsyn","G_tv2_nonsyn")



filename<-paste0("Output1A/MutFreq/Maj/Tranversion_Mut_Freq_Histograms_mean.pdf")
pdf(filename, width = 12, height = 7)
par(mfrow=c(2,4))
par(mar = c(4,4.3,2.5,1))
breaks=seq(-6,0,by=.2)

for (i in 1:8){
        datavector1<-get(datanames1[i])
        datavector2<-get(datanames2[i])
        datavector1<-c(datavector1,datavector2)
        nt<-substr(paste(datanames1[i]),start=1,stop=1)
        typ<-substr(paste(datanames1[i]),start=7, stop=12)
        hist(log10(datavector1), breaks = breaks, xlim=c(-5,-0),col=colors[i], 
             main = paste0(nt, " : ", typ), ylab="Count", xlab="",
             xaxt="n", cex.main=1.7, las=1, cex.axis =1,cex.lab=1.5)
        abline(v=log10(mean(datavector1, na.rm = TRUE)),col="white",lwd=2,lty=1)
        abline(v=log10(mean(datavector1, na.rm = TRUE)),lwd=2,lty=2)
        x=seq(-5,0,by=1)
        labels <- sapply(x,function(y)as.expression(bquote(10^ .(y))))
        axis(1, at=x,labels=labels, col.axis="black", line=-.7, cex.axis=1)
        mtext(text = "Transversion Mutation Frequency", side = 1, line = 1.5, cex=0.9)
}

dev.off()

filename<-paste0("Output1A/MutFreq/Maj/Tranversion_Mut_Freq_Histograms_median.pdf")
pdf(filename, width = 12, height = 7)
par(mfrow=c(2,4))
par(mar = c(4,4.3,2.5,1))
breaks=seq(-6,0,by=.2)
for (i in 1:8){
        datavector1<-get(datanames1[i])
        datavector2<-get(datanames2[i])
        datavector1<-c(datavector1,datavector2)
        nt<-substr(paste(datanames1[i]),start=1,stop=1)
        typ<-substr(paste(datanames1[i]),start=7, stop=12)
        hist(log10(datavector1), breaks = breaks, xlim=c(-5,-0),col=colors[i], 
             main = paste0(nt, " : ", typ), ylab="Count", xlab="",
             xaxt="n", cex.main=1.7, las=1, cex.axis =1,cex.lab=1.5)
        abline(v=median(log10(datavector1), na.rm = TRUE),col="white",lwd=2,lty=1)
        abline(v=median(log10(datavector1), na.rm = TRUE),lwd=2,lty=2,col='red')
        x=seq(-5,0,by=1)
        labels <- sapply(x,function(y)as.expression(bquote(10^ .(y))))
        axis(1, at=x,labels=labels, col.axis="black", line=-.7, cex.axis=1)
        mtext(text = "Transversion Mutation Frequency", side = 1, line = 1.5, cex=0.9)
}
dev.off()



