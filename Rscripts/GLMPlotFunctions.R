nucord <- c("a", "t", "c", "g")

t.col <- function(color, percent = 50, name = NULL) {
        rgb.val <- col2rgb(color)
        t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, alpha = (100-percent)*255/100)
        return(t.col )
}

makePlot.axisbreak <- function(main){
        
        plot(0, type = "n", xlim = c(.5, 4.5), ylim = c(0.00004, .1), axes = FALSE, ylab = "Mutation frequency", xlab = "Mutation type", main = main, log = "y")
        
        col.par <- "gray95"
        #    abline(h = seq(0, .15, by = .025), col = col.par)
        abline(h = (1:9)*10^(-1), col = col.par)
        abline(h = (1:9)*10^(-2), col = col.par)
        abline(h = (1:9)*10^(-3), col = col.par)
        abline(h = (1:9)*10^(-4), col = col.par)
        abline(h = (1:2)*10^(-5), col = col.par)
        abline(h = (1)*10^(-4.4), col = col.par)
        box()
        
        axislabs <- c(expression("A" %->% "G"), expression("T" %->% "C"), expression("C" %->% "T"), expression("G" %->% "A"))
        axis(1, at = 1:4, axislabs)
        eaxis(2, at = 10^c(-1*(0:4)))
        axis(2, at = 10^(-4.4), label = c("0"), las = 2)
        axis.break(2,2*10^-(4.5),style="slash")
        
        
}


makeDataFrameToModify <- function(NsOrNo = 0, CpGorNo = 0, bigAAChangeOrNo = 0){
        
        atcg.mat <- as.data.frame(matrix(data = 0, ncol = ncol(glmData)-1, nrow = 4))
        #ATCG elements
        diag(atcg.mat[1:4, 1:4]) <- 1
        #Reserve the first column for intercept
        atcg.mat[,1] <- 1
        #synonymous or nonsynonymous mutation?
        atcg.mat[,6] <- NsOrNo
        #bigAA change?
        atcg.mat[,7] <- bigAAChangeOrNo * atcg.mat[,6]
        #CpG mutation or not
        atcg.mat[1:2,5] <- CpGorNo
        atcg.mat[3:4,5] <- 0
        
        #nonysynonymous interactions with a t c g
        atcg.mat[8:10] <- atcg.mat[,2:4] * atcg.mat[,6]
        
        #names
        names(atcg.mat) <- c(rownames(modcoef) )
        
        return(as.matrix(atcg.mat))
} 


## Function to plot the observed mut freq data 
plotDat <- function(NsOrNo, CpGorNo, bigAAChangeOrNo, colVal, offsetval){
        
        pchpar <- 16
        cexpar <- .6
        
        nsinds <- which(d$Type == c("syn", "nonsyn")[NsOrNo + 1])
        cpginds <- which(d$makesCpG == CpGorNo)
        aachangeinds <- which(d$bigAAChange == bigAAChangeOrNo)
        
        datInds <- intersect(intersect(nsinds, cpginds), aachangeinds)
        
        xinds <- as.character(d$ref[datInds])
        for(i in 1:4){
                xinds[xinds == nucord[i]] <- i
        }
        xinds <- as.numeric(xinds)
        yinds <- d$mean[datInds]
        xjit <- rnorm(length(xinds), offsetval, .03) #jitter
        
        yinds[yinds == 0] <- 10^(-4.4)
        
        points(xinds + xjit, yinds, col = t.col(colVal),  pch = pchpar, cex = cexpar)
}

#Function to plot the estimated mut freq from GLM/beta regreassion
plotVals <- function(NsOrNo, CpGorNo, bigAAChangeOrNo, colVal, offset){
        
        pchpar <- 21
        cexpar <- 1
        arrowlwd <- 1.5
        arrowlen <-.08
        
        setUpDat <- makeDataFrameToModify(NsOrNo, CpGorNo, bigAAChangeOrNo)[,rownames(modcoef)]
        
        #compute
        #%*% = matrix manipulation
        pointToPlot <- exp(setUpDat %*% coef.vals)
        arrowToPlotLow <- exp(setUpDat %*% coef.pSE.vals)
        arrowToPlotHigh <- exp(setUpDat %*% coef.mSE.vals)
        
        #don't plot estimates for CpG sites at Cs or Gs
        restrict <- 1:4
        if(CpGorNo == 1){
                restrict <- 1:2
        }
        
        #plot
        arrows((1:4)[restrict] + offset, arrowToPlotLow[restrict], (1:4)[restrict] + offset, arrowToPlotHigh[restrict], length = arrowlen, code = 3, angle = 90, lwd = arrowlwd, col = "black")
        points((1:4)[restrict] + offset, pointToPlot[restrict], col ="black", bg = substr(colVal,start=1, stop=7),  pch = pchpar, cex = .8)
        #arrows((1:4)[restrict] + offset, arrowToPlotLow[restrict], (1:4)[restrict] + offset, arrowToPlotHigh[restrict], length = arrowlen, code = 3, angle = 90, lwd = arrowlwd, col = "black")
        
}
