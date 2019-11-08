#Read in the data file and convert the first col to rownames
library(dplyr)
source("Rscripts/baseRscript.R")

########################################
plotter <- function(d){
        par(mar = c(4.5, 4.5, 2, 2))
        main.dat <- d
        wheresthebreak <- 5
        remap = 1       
        dataset <- main.dat[complete.cases(main.dat), ]
        toPlot <- dataset[,"Freq"]
        toPlot <- toPlot[!is.na(toPlot)]
        colVect <- rep(0, nrow(dataset))
        colVect[dataset$TypeOfSite == "nonsyn"] <- cols[5]
        colVect[dataset$TypeOfSite == "syn"] <- cols[3]
        colVect[dataset$TypeOfSite == "stop"] <- "black"
        plot(5, type = "n", log = "y", axes = FALSE, xlim = c(0, length(toPlot)), 
             ylim = c(10^-(wheresthebreak),0.1),  
             ylab = "Average mutation frequency", xlab = "Mutations ordered by mean mutation frequency",
             cex.lab = 1.3)
        for(i in 1:wheresthebreak){
                abline(h = 1:10 * 10^(-i), col = "gray70")
        }
        abline(h = 10^-(wheresthebreak), col = "gray70")
        if(remap == 1){
                eaxis(side = 2, at = 10^((-1):(-(wheresthebreak-1))))
                axis(side = 2, at = c(10^-(wheresthebreak)), label = c(0), las = 2)
                box()
                axis.break(2,2*10^-(wheresthebreak),style="slash")
        }else{
                eaxis(side = 2, at = 10^((0):(-(wheresthebreak))))
                axis(side = 2, at = 1, label = 1, las =2)
                box()
        }
        cexval <- 1.5
        toPlot[toPlot == 0] <- 10^-(wheresthebreak)
        points(1:length(toPlot), sort(toPlot), col = colVect[order(toPlot)], pch = "|", cex = cexval)
        axis(1)
        legend("bottomright", c("Synonymous", "Non-synonymous", "Nonsense"), col = c(cols[3], cols[5], "black"), pch = "|", bg = "white", pt.cex = cexval)
}

####

filedate<-"2018-10-27"
TsMutFreq<-read.csv(paste0("Output/MutFreq/FilteredSummary_Transition_mutation_frequencies_",filedate,".csv"))
TsMutFreq<-TsMutFreq[,-1]
t1<-TsMutFreq
t1$mean<-apply(t1[2:(s+1)], 1, mean, na.rm=T)
range(t1$mean, na.rm=T) #[1] 0.0000125248 0.0504393357

#Count the number of NAs per row
t1$NofNA<-apply(t1[2:(s+1)], 1, function(x)sum(is.na(x)))
#Remove the rows that has more than 50% NAs
t1<-t1[which(t1$NofNA<0.5*s),]
t2<-t1[,c("pos","mean")]


dat<-read.csv("Output/Overview_filtered/D75002-_overview3.csv",stringsAsFactors=FALSE)
types<-dat[,c("pos","ref","Type","WTAA","MUTAA","bigAAChange")]

MF<-merge(types,t2,by='pos')



# create figures for transition mutations
        OverviewData<-MF
        colnames(OverviewData)[colnames(OverviewData)=="Type"]<-"TypeOfSite"
        colnames(OverviewData)[colnames(OverviewData)=="mean"]<-"Freq"
        pdf(paste0("Output/SummaryFig.Filtered/Fig.1_MutFreq-ordered.pdf"), height = 6, width = 9)
        plotter(OverviewData)
        dev.off() 
      


######### for combined transversion plots, skip to the below section ####
# create figures for transversion mutations #1
for (i in 1:length(Overview_summary)){
        OverviewData<-Overview_summary[[i]]
        colnames(OverviewData)[colnames(OverviewData)=="Type.tv1"]<-"TypeOfSite"
        colnames(OverviewData)[colnames(OverviewData)=="freq.transv.1"]<-"Freq"
        pdf(paste0("Output/Fig.1/MutFreq-ordered-Transv1_", names(Overview_summary)[i], ".pdf"), height = 6, width = 9)
        plotter(OverviewData)
        dev.off() 
}  

# create figures for transversion mutations #2
for (i in 1:length(Overview_summary)){
        OverviewData<-Overview_summary[[i]]
        colnames(OverviewData)[colnames(OverviewData)=="Type.tv2"]<-"TypeOfSite"
        colnames(OverviewData)[colnames(OverviewData)=="freq.transv.2"]<-"Freq"
        pdf(paste0("Output/Fig.1/MutFreq-ordered-Transv2_", names(Overview_summary)[i], ".pdf"), height = 6, width = 9)
        plotter(OverviewData)
        dev.off() 
}

######################################################

# Create figures for all transversion mutations
plotter.Tvs <- function(d){
        par(mar = c(4.5, 4.5, 2, 2))
        main.dat <- d
        wheresthebreak <- 6
        remap = 1       
        dataset <- main.dat[complete.cases(main.dat), ]
        toPlot <- dataset[,"Freq"]
        toPlot <- toPlot[!is.na(toPlot)]
        colVect <- rep(0, nrow(dataset))
        colVect[dataset$TypeOfSite == "nonsyn"] <- cols[5]
        colVect[dataset$TypeOfSite == "syn"] <- cols[3]
        colVect[dataset$TypeOfSite == "stop"] <- "black"
        colVect[dataset$TypeOfSite == "mixed"] <- "red"
        plot(5, type = "n", log = "y", axes = FALSE, xlim = c(0, length(toPlot[!is.na(toPlot)])), 
             ylim = c(10^-(wheresthebreak), max(toPlot, na.rm = TRUE)),  
             ylab = "Mean mutation frequency", xlab = "Mutations ordered by mean mutation frequency",
             cex.lab = 1.3)
        for(i in 1:wheresthebreak){
                abline(h = 1:10 * 10^(-i), col = "gray70")
        }
        abline(h = 10^-(wheresthebreak), col = "gray70")
        if(remap == 1){
                eaxis(side = 2, at = 10^((-1):(-(wheresthebreak-1))))
                axis(side = 2, at = c(1, 10^-(wheresthebreak)), label = c(1, 0), las = 2)
                box()
                axis.break(2,2*10^-(wheresthebreak),style="slash")
        }else{
                eaxis(side = 2, at = 10^((0):(-(wheresthebreak))))
                axis(side = 2, at = 1, label = 1, las =2)
                box()
        }
        cexval <- 1.2
        toPlot[toPlot == 0] <- 10^-(wheresthebreak)
        points(1:length(toPlot), sort(toPlot), col = colVect[order(toPlot)], pch = "|", cex = cexval)
        axis(1)
        legend("bottomright", c("Synonymous", "Non-synonymous","Mixed", "Nonsense"), col = c(cols[3], cols[5], "red", "black"), pch = "|", bg = "white", pt.cex = cexval)
}



for (i in 1:length(Overview_summary)){
        OverviewData<-Overview_summary[[i]]
        for (j in 1:nrow(OverviewData)){ 
                if (is.na(OverviewData$Type.tv1[j])|is.na(OverviewData$Type.tv2[j])) next
                else {if (OverviewData$Type.tv1[j] == OverviewData$Type.tv2[j])  OverviewData$TypeOfSite[j]<-OverviewData$Type.tv1[j]
                if (OverviewData$Type.tv1[j] != OverviewData$Type.tv2[j]) OverviewData$TypeOfSite[j]<-'mixed'}
        }        
        colnames(OverviewData)[colnames(OverviewData)=="freq.transv"]<-"Freq"
        
        pdf(paste0("Output/Fig.1/MutFreq-ordered-Transv_", names(Overview_summary)[i], ".pdf"), height = 6, width = 9)
        plotter.Tvs(OverviewData)
        dev.off() 
}




