#make AA transition figure
#Just NS

source("Rscripts/baseRscript.R")


computeColInds <- function(wtnt){
        colsToRet <- wtnt
        colsToRet[wtnt == "a"] <- 1
        colsToRet[wtnt == "t"] <- 2
        colsToRet[wtnt == "c"] <- 3
        colsToRet[wtnt == "g"] <- 4
        return(as.numeric(colsToRet))
}

#####
HCVFiles_overview3<-list.files("Output/Overview_filtered/",pattern="overview3.csv")
s<-length(HCVFiles_overview3)


filtered<-list()
fnames<-c("Ts", "Ts_NA", "Ts_zero")
for (i in 1:3){
        filename<-fnames[i]
        df<-read.csv(paste0("Output/Mut.freq.filtered/Summary_",filename,".Q35.csv"),stringsAsFactors = F)
        df<-df[,-1]
        filtered[[i]]<-df
        names(filtered)[i]<-filename
}




#create a figure  ylim =-10^-5

        hcvdat<-filtered[[1]]

        NS <- hcvdat[hcvdat$Type == "nonsyn",]
        stotab <- table(NS$WTAA, NS$MUTAA)
        allComps <- matrix(data = NA, nrow = sum(stotab > 0), ncol = 2) #Get a list of all AA substitutions in the dataset
        
        for(i in 1:nrow(stotab)){
                for(j in 1:ncol(stotab)){
                        if(stotab[i,j] != 0){ allComps[min(which(is.na(allComps[,1]))), ] <- c(rownames(stotab)[i], colnames(stotab)[j] ) }
                }
        }

        opac = 50
        pchs <- c(16, 17)
        ylimval <- -5
       
        pdf(paste0("Output/SummaryFig.Filtered/AAchanges_Ts.pdf"), height = 5, width = 9)
        par(mar = c(5, 4, 1, 5))
        
        plot(0, type = "n", xlim = c(1, nrow(allComps)), ylim = c(ylimval, -1), axes = FALSE, ylab = "Average mutation frequency", xlab = "Amino acid change")
        axis(side = 2, at = seq(-1,ylimval, by=-1), labels = expression(10^-1, 10^-2, 10^-3, 10^-4,10^-5), las = 2)
        labsize <- .9
        linepos=5.4
        yesnoswitch=1
        
        for(i in unique(allComps[,1])){
                rangevals <- range(which(allComps[,1] == i))
                segments(rangevals[1] - .25, -linepos, rangevals[2] + .25, -linepos, col="black", xpd=NA)
                if (yesnoswitch==1) {polygon(c(rangevals[1] - .25,rangevals[2] + .25, rangevals[2] + .25,rangevals[1] - .25),
                                               c(-linepos-.01,-linepos-.01,-linepos-.2,-linepos-.2),col="#FFE7DE",border=0,xpd=NA)
                                for (j in rangevals[1]:rangevals[2]){
                                        abline(v=j,col = "#FFE7DE", lty=1, lwd=14)}}
                                        
                                if (yesnoswitch==1)  mtext(paste(i), 1, line = 1, at = mean(rangevals), cex = labsize,col=1)
                if (yesnoswitch==-1) mtext(paste(i), 1, line = 1, at = mean(rangevals), cex = labsize,col=1)
                yesnoswitch=yesnoswitch*-1
        }
        abline(v = 1:nrow(allComps), col = "grey90", lty=c(1,2,3))
        abline(v = 1:nrow(allComps), col = "grey90", lty=1, lwd=3)
        abline(h = 0, col = "grey90", lwd=2)
        box()
        
        for(i in 1:nrow(allComps)){
                relInds <- intersect(which(NS$WTAA == allComps[i, 1]), which(NS$MUTAA == allComps[i, 2]))
                relDat <- NS[relInds,]
                xjit <- rnorm(length(relInds), 0, .1)
                points(rep(i, length(relInds)) + xjit, log10(relDat$mean), 
                       pch = pchs[relDat$bigAAChange[1]+1], col = cols[computeColInds(relDat$ref)], cex = .75)
        }
        
        mtext(side = 1, at = 1:length(allComps[,2]), text = paste(allComps[,2]), cex = labsize*.8)
        mtext("Mutant", 1, line = 0, at = -1, cex = labsize, adj = 1)
        
        segments(-6 , -linepos, -0.5, -linepos, col="black", xpd=NA)
        mtext("Ancestral", 1, line = 1, at = -1, cex = labsize, adj = 1)
        
        leg.location<-nrow(allComps)+3
        legend(leg.location, -1, c("A", "T", "C", "G"), col = c(cols, "black", "black"), pch = c(15, 15, 15, 15), xpd = NA, box.lty = "blank", title = "Ancestral\nnucleotide")
        legend(leg.location, -2.5, c("No", "Yes"), pch = pchs, xpd = NA, border = "white", box.lwd = 0, title = "Drastic AA\nchange", box.lty = "blank")
        
dev.off()





### plot 10^-1 to 10^-4
#create a figure

hcvdat<-filtered[[3]]


opac = 50
pchs <- c(16, 17)
ylimval <- -4
pdf(paste0("Output/SummaryFig.Filtered/AAchanges_Ts.pdf"), height = 5, width = 9)
pdf(paste0("Output/SummaryFig.Filtered/AAchanges_Ts_NA.pdf"), height = 5, width = 9)
pdf(paste0("Output/SummaryFig.Filtered/AAchanges_Ts_zero.pdf"), height = 5, width = 9)


par(mar = c(5, 4, 1, 5))


plot(0, type = "n", xlim = c(1, nrow(allComps)), ylim = c(ylimval, -1), axes = FALSE, ylab = "Average mutation frequency", xlab = "Mutation")
axis(side = 2, at = seq(-1,ylimval, by=-1), labels = expression(10^-1, 10^-2, 10^-3, 10^-4), las = 2)
labsize <- .9
linepos=4.3
yesnoswitch=1
for(i in unique(allComps[,1])){
        rangevals <- range(which(allComps[,1] == i))
        segments(rangevals[1] - .25, -linepos, rangevals[2] + .25, -linepos, col="black", xpd=NA)
        if (yesnoswitch==1) {polygon(c(rangevals[1] - .25,rangevals[2] + .25, rangevals[2] + .25,rangevals[1] - .25),
                                     c(-linepos-.01,-linepos-.01,-linepos-.2,-linepos-.2),col="#FFE7DE",border=0,xpd=NA)
                for (j in rangevals[1]:rangevals[2]){
                        abline(v=j,col = "#FFE7DE", lty=1, lwd=14)}}
        
        if (yesnoswitch==1)  mtext(paste(i), 1, line = 1, at = mean(rangevals), cex = labsize,col=1)
        if (yesnoswitch==-1) mtext(paste(i), 1, line = 1, at = mean(rangevals), cex = labsize,col=1)
        yesnoswitch=yesnoswitch*-1
}
abline(v = 1:nrow(allComps), col = "grey90", lty=c(1,2,3))
abline(v = 1:nrow(allComps), col = "grey90", lty=1, lwd=3)
abline(h = 0, col = "grey90", lwd=2)
box()

for(i in 1:nrow(allComps)){
        relInds <- intersect(which(NS$WTAA == allComps[i, 1]), which(NS$MUTAA == allComps[i, 2]))
        relDat <- NS[relInds,]
        xjit <- rnorm(length(relInds), 0, .1)
        points(rep(i, length(relInds)) + xjit, log10(relDat$mean), 
               pch = pchs[relDat$bigAAChange[1]+1], col = cols[computeColInds(relDat$ref)], cex = .75)
}

mtext(side = 1, at = 1:length(allComps[,2]), text = paste(allComps[,2]), cex = labsize*.8)
mtext("Mutant", 1, line = 0, at = -1, cex = labsize, adj = 1)

segments(-5 , -linepos, -0.5, -linepos, col="black", xpd=NA)
mtext("Ancestral", 1, line = 1, at = -1, cex = labsize, adj = 1)

leg.location<-nrow(allComps)+3
legend(leg.location, -1, c("A", "T", "C", "G"), col = c(cols, "black", "black"), pch = c(15, 15, 15, 15), xpd = NA, box.lty = "blank", title = "Ancestral\nnucleotide")
legend(leg.location, -2.5, c("No", "Yes"), pch = pchs, xpd = NA, border = "white", box.lwd = 0, title = "Drastic AA\nchange", box.lty = "blank")

dev.off()




        
###################
#plot with SelCoeff as Y axis
d<-FilteredOverview[[1]]
d1<-d[d$pos>=342,c("pos","TSmutrate")]
SC<-merge(TsMutFreq,d1, by="pos")

for (i in 1:nrow(SC)){
        SC$EstSC[i]<-EstimatedS(SC$TSmutrate[i],SC$mean[i])
}

sc_g<-SC[SC$Type=="nonsyn"&SC$ref=="g",]

######################
#create a figure


SC1 <- SC[SC$Type == "nonsyn",]
stotab <- table(SC1$WTAA, SC1$MUTAA)
allComps <- matrix(data = NA, nrow = sum(stotab > 0), ncol = 2) #Get a list of all AA substitutions in the dataset

for(i in 1:nrow(stotab)){
        for(j in 1:ncol(stotab)){
                if(stotab[i,j] != 0){ allComps[min(which(is.na(allComps[,1]))), ] <- c(rownames(stotab)[i], colnames(stotab)[j] ) }
        }
}

opac = 50
pchs <- c(16, 17)
ylimval <- -4

pdf(paste0("Output/SummaryFig.Filtered/AAchanges_transition_selCoeff.pdf"), height = 5, width = 9)
par(mar = c(5, 4, 1, 5))

plot(0, type = "n", xlim = c(1, nrow(allComps)), ylim = c(ylimval, 0), axes = FALSE, ylab = "Estimated selection coefficient", xlab = "Amino acid mutation")
axis(side = 2, at = seq(0,ylimval, by=-1), labels = expression(10^-1, 10^-2, 10^-3, 10^-4,10^-5), las = 2)
labsize <- .9
linepos=4.4
yesnoswitch=1
for(i in unique(allComps[,1])){
        rangevals <- range(which(allComps[,1] == i))
        segments(rangevals[1] - .25, -linepos, rangevals[2] + .25, -linepos, col="black", xpd=NA)
        if (yesnoswitch==1) { polygon(c(rangevals[1] - .25,rangevals[2] + .25, rangevals[2] + .25,rangevals[1] - .25),
                                       c(-linepos-.01,-linepos-.01,-linepos-.2,-linepos-.2),col="#FFE7DE",border=0,xpd=NA)
                for (j in rangevals[1]:rangevals[2]){
                        abline(v=j,col = "#FFE7DE", lty=1, lwd=14)}}
        if (yesnoswitch==1)  mtext(paste(i), 1, line = 1, at = mean(rangevals), cex = labsize,col=1)
        if (yesnoswitch==-1) mtext(paste(i), 1, line = 1, at = mean(rangevals), cex = labsize,col=1)
        yesnoswitch=yesnoswitch*-1
}
abline(v = 1:nrow(allComps), col = "grey90", lty=c(1,2,3))
abline(v = 1:nrow(allComps), col = "grey90", lty=1, lwd=3)
abline(h = 0, col = "grey90", lwd=2)
box()

for(i in 1:nrow(allComps)){
        relInds <- intersect(which(SC$WTAA == allComps[i, 1]), which(SC$MUTAA == allComps[i, 2]))
        relDat <- SC[relInds,]
        xjit <- rnorm(length(relInds), 0, .1)
        points(rep(i, length(relInds)) + xjit, log10(relDat$EstSC), 
               pch = pchs[relDat$bigAAChange[1]+1], col = cols[computeColInds(relDat$ref)], cex = .75)
}

mtext(side = 1, at = 1:length(allComps[,2]), text = paste(allComps[,2]), cex = labsize*.8)
mtext("Mutant", 1, line = 0, at = -1, cex = labsize, adj = 1)

segments(-6 , -linepos, -0.5, -linepos, col="black", xpd=NA)
mtext("Ancestral", 1, line = 1, at = -1, cex = labsize, adj = 1)

leg.location<-nrow(allComps)+3
legend(leg.location, -1, c("A", "T", "C", "G"), col = c(cols, "black", "black"), pch = c(15, 15, 15, 15), xpd = NA, box.lty = "blank", title = "Ancestral\nnucleotide")
legend(leg.location, -2.5, c("No", "Yes"), pch = pchs, xpd = NA, border = "white", box.lwd = 0, title = "Drastic AA\nchange", box.lty = "blank")

dev.off()



###################
#Eliminate HVR1 from the dataset:
filedate<-"2018-10-27"
TsMutFreq<-read.csv(paste0("Output/MutFreq/FilteredSummary_Transition_mutation_frequencies_",filedate,".csv"))
TsMutFreq<-TsMutFreq[,-1]

genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
TsMutFreq2<-TsMutFreq[-c(genes$start[4]:genes$end[4]),]

t1<-TsMutFreq2
t1$mean<-apply(t1[2:(s+1)], 1, mean, na.rm=T)
range(t1$mean, na.rm=T) 

#Count the number of NAs per row
t1$NofNA<-apply(t1[2:(s+1)], 1, function(x)sum(is.na(x)))
#Remove the rows that has more than 50% NAs
t1<-t1[which(t1$NofNA<0.5*s),]
t2<-t1[,c("pos","mean")]


dat<-Overview_fil[[1]]
types<-dat[,c("pos","ref","Type","WTAA","MUTAA","bigAAChange")]

MF<-merge(types,t2,by='pos')



opac = 50
pchs <- c(16, 17)
ylimval <- -4
pdf(paste0("Output/SummaryFig.Filtered/AAchanges_transition_noHVR1.pdf"), height = 5, width = 9)
par(mar = c(5, 4, 1, 5))
plot(0, type = "n", xlim = c(1, nrow(allComps)), ylim = c(ylimval, 0), axes = FALSE, ylab = "Average mutation frequency", xlab = "Mutation")
axis(side = 2, at = seq(0,ylimval, by=-1), labels = expression(10^-0,10^-1, 10^-2, 10^-3, 10^-4), las = 2)
labsize <- .9
linepos=4.4
yesnoswitch=1
for(i in unique(allComps[,1])){
        rangevals <- range(which(allComps[,1] == i))
        segments(rangevals[1] - .25, -linepos, rangevals[2] + .25, -linepos, col="black", xpd=NA)
        if (yesnoswitch==1) {polygon(c(rangevals[1] - .25,rangevals[2] + .25, rangevals[2] + .25,rangevals[1] - .25),
                                     c(-linepos-.01,-linepos-.01,-linepos-.2,-linepos-.2),col="#FFE7DE",border=0,xpd=NA)
                for (j in rangevals[1]:rangevals[2]){
                        abline(v=j,col = "#FFE7DE", lty=1, lwd=14)}}
        
        if (yesnoswitch==1)  mtext(paste(i), 1, line = 1, at = mean(rangevals), cex = labsize,col=1)
        if (yesnoswitch==-1) mtext(paste(i), 1, line = 1, at = mean(rangevals), cex = labsize,col=1)
        yesnoswitch=yesnoswitch*-1
}
abline(v = 1:nrow(allComps), col = "grey90", lty=c(1,2,3))
abline(v = 1:nrow(allComps), col = "grey90", lty=1, lwd=3)
abline(h = 0, col = "grey90", lwd=2)
box()

for(i in 1:nrow(allComps)){
        relInds <- intersect(which(NS$WTAA == allComps[i, 1]), which(NS$MUTAA == allComps[i, 2]))
        relDat <- NS[relInds,]
        xjit <- rnorm(length(relInds), 0, .1)
        points(rep(i, length(relInds)) + xjit, log10(relDat$mean), 
               pch = pchs[relDat$bigAAChange[1]+1], col = cols[computeColInds(relDat$ref)], cex = .75)
}

mtext(side = 1, at = 1:length(allComps[,2]), text = paste(allComps[,2]), cex = labsize*.8)
mtext("Mutant", 1, line = 0, at = -1, cex = labsize, adj = 1)

segments(-5 , -linepos, -0.5, -linepos, col="black", xpd=NA)
mtext("Ancestral", 1, line = 1, at = -1, cex = labsize, adj = 1)

leg.location<-nrow(allComps)+3
legend(leg.location, -1, c("A", "T", "C", "G"), col = c(cols, "black", "black"), pch = c(15, 15, 15, 15), xpd = NA, box.lty = "blank", title = "Ancestral\nnucleotide")
legend(leg.location, -2.5, c("No", "Yes"), pch = pchs, xpd = NA, border = "white", box.lwd = 0, title = "Drastic AA\nchange", box.lty = "blank")

dev.off()

   
###################
#Eliminate HVR1 from the dataset: Plot with Selection Coefficients:
filedate<-"2018-10-27"
Tsselcoef<-read.csv(paste0("Output/SelCoeff/FilteredSummary_Transition_SelCoeff_",filedate,".csv"))
Tsselcoef<-Tsselcoef[,-1]

genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
Tsselcoef2<-Tsselcoef[-c(genes$start[4]:genes$end[4]),]



t1<-Tsselcoef2
t1$mean<-apply(t1[2:(s+1)], 1, mean, na.rm=T)
range(t1$mean, na.rm=T) 


#Count the number of NAs per row
t1$NofNA<-apply(t1[2:(s+1)], 1, function(x)sum(is.na(x)))
#Remove the rows that has more than 50% NAs
t1<-t1[which(t1$NofNA<0.5*s),]
t2<-t1[,c("pos","mean")]


dat<-Overview_fil[[1]]
types<-dat[,c("pos","ref","Type","WTAA","MUTAA","bigAAChange")]

SC<-merge(types,t2,by='pos')



opac = 50
pchs <- c(16, 17)
ylimval <- -4
pdf(paste0("Output/SummaryFig.Filtered/AAchanges_transition_noHVR2.pdf"), height = 5, width = 9)
par(mar = c(5, 4, 1, 5))
plot(0, type = "n", xlim = c(1, nrow(allComps)), ylim = c(ylimval, 0), axes = FALSE, ylab = "Average selection coefficient", xlab = "Mutation")
axis(side = 2, at = seq(0,ylimval, by=-1), labels = expression(10^-0,10^-1, 10^-2, 10^-3, 10^-4), las = 2)
labsize <- .9
linepos=4.4
yesnoswitch=1
for(i in unique(allComps[,1])){
        rangevals <- range(which(allComps[,1] == i))
        segments(rangevals[1] - .25, -linepos, rangevals[2] + .25, -linepos, col="black", xpd=NA)
        if (yesnoswitch==1) {polygon(c(rangevals[1] - .25,rangevals[2] + .25, rangevals[2] + .25,rangevals[1] - .25),
                                     c(-linepos-.01,-linepos-.01,-linepos-.2,-linepos-.2),col="#FFE7DE",border=0,xpd=NA)
                for (j in rangevals[1]:rangevals[2]){
                        abline(v=j,col = "#FFE7DE", lty=1, lwd=14)}}
        
        if (yesnoswitch==1)  mtext(paste(i), 1, line = 1, at = mean(rangevals), cex = labsize,col=1)
        if (yesnoswitch==-1) mtext(paste(i), 1, line = 1, at = mean(rangevals), cex = labsize,col=1)
        yesnoswitch=yesnoswitch*-1
}
abline(v = 1:nrow(allComps), col = "grey90", lty=c(1,2,3))
abline(v = 1:nrow(allComps), col = "grey90", lty=1, lwd=3)
abline(h = 0, col = "grey90", lwd=2)
box()

for(i in 1:nrow(allComps)){
        relInds <- intersect(which(NS$WTAA == allComps[i, 1]), which(NS$MUTAA == allComps[i, 2]))
        relDat <- NS[relInds,]
        xjit <- rnorm(length(relInds), 0, .1)
        points(rep(i, length(relInds)) + xjit, log10(relDat$mean), 
               pch = pchs[relDat$bigAAChange[1]+1], col = cols[computeColInds(relDat$ref)], cex = .75)
}

mtext(side = 1, at = 1:length(allComps[,2]), text = paste(allComps[,2]), cex = labsize*.8)
mtext("Mutant", 1, line = 0, at = -1, cex = labsize, adj = 1)

segments(-5 , -linepos, -0.5, -linepos, col="black", xpd=NA)
mtext("Ancestral", 1, line = 1, at = -1, cex = labsize, adj = 1)

leg.location<-nrow(allComps)+3
legend(leg.location, -1, c("A", "T", "C", "G"), col = c(cols, "black", "black"), pch = c(15, 15, 15, 15), xpd = NA, box.lty = "blank", title = "Ancestral\nnucleotide")
legend(leg.location, -2.5, c("No", "Yes"), pch = pchs, xpd = NA, border = "white", box.lwd = 0, title = "Drastic AA\nchange", box.lty = "blank")

dev.off()

        
        
        
        