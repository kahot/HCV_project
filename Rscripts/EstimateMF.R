# Function to estimate the mutation frequency at a site
EstimateMF<- function(NsOrNo = 0, CpGorNo = 0, bigAAChangeOrNo = 0,CoreorNo=0,E1orNo=0, HVRorNo=0,E2orNo=0,NS1orNo=0,NS2orNo=0,NS4AorNo=0,NS5BorNo=0){
        atcg.mat <- as.data.frame(matrix(data = 0, ncol = 18, nrow = 4))
        #ATCG elements
        diag(atcg.mat[1:4, 1:4]) <- 1
        #Reserve the first column for intercept
        atcg.mat[,1] <- 1
        #CpG mutation or not
        atcg.mat[1:2,5] <- CpGorNo
        atcg.mat[3:4,5] <- 0
        #synonymous or nonsynonymous mutation?
        atcg.mat[,6] <- NsOrNo
        #bigAA change?
        atcg.mat[,7] <- bigAAChangeOrNo * atcg.mat[,6]
        
        #nonysynonymous interactions with a t c g
        atcg.mat[8:10] <- atcg.mat[,2:4] * atcg.mat[,6]
        #genes
        atcg.mat[,11] <- CoreorNo
        atcg.mat[,12] <- E1orNo
        atcg.mat[,13] <- HVRorNo
        atcg.mat[,14] <- E2orNo
        atcg.mat[,15] <- NS1orNo
        atcg.mat[,16] <- NS2orNo
        atcg.mat[,17] <- NS4AorNo
        atcg.mat[,18] <- NS5BorNo
        
        #names
        names(atcg.mat) <- c(rownames(modcoef) )
        
        return(as.matrix(atcg.mat))
} 
