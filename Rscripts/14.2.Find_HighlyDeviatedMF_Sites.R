####
HCVFiles_overview3<-list.files("Output/Overview_filtered/",pattern="overview3.csv")
FilteredOverview<-list()
for (i in 1:length(HCVFiles_overview3)){ 
        overviews<-read.csv(paste0("Output/Overview_filtered/",HCVFiles_overview3[i]),stringsAsFactors=FALSE)
        overviews<-overviews[,-1]
        FilteredOverview[[i]]<-overviews
        names(FilteredOverview)[i]<-substr(paste(HCVFiles_overview3[i]),start=1,stop=7)
}

############

data<-read.csv("Output/GLM/GlmdataFull.Ts.Q35.csv")
data<-data[,-1]
s<-length(FilteredOverview)


# from 10.3_Filtered OveriewData_GLM.R
#Read data file
GlmData <-read.csv("Output/GLM/GlmdataFull.Ts.Q35.csv")
GlmData1<-GlmData[,-1]

modcoef<-read.csv("Output/GLM/BetaReg_mod.g2_Ts.Q35.csv",stringsAsFactors = F)
rownames(modcoef)<-modcoef$X
modcoef<-modcoef[,-1]
coef.vals <- modcoef[,1]


#from GLMPlotFunctions.R
makeDataFrameToModify1 <- function(NsOrNo = 0, CpGorNo = 0, bigAAChangeOrNo = 0,CoreorNo=0,E1orNo=0, HVRorNo=0,E2orNo=0,NS1orNo=0,NS2orNo=0,NS4AorNo=0,NS5BorNo=0){
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


data$EstiamtedMF<-0
for (i in 1:nrow(data)){
        NsOrNo=data$Nonsyn[i]
        CpGorNo = data$CpG[i]
        bigAAChangeOrNo =data$bigAAChange[i]
        CoreorNo<-data$Core[i]
        E1orNo<-data$E1[i]
        HVRorNo = data$HVR1[i]
        E2orNo<-data$E2[i]
        NS1orNo<-data$NS1[i]
        NS2orNo<-data$NS2[i]
        NS4AorNo<-data$NS4A[i]
        NS5BorNo<-data$NS5B[i]
        
        setUpDat <- makeDataFrameToModify1(NsOrNo, CpGorNo, bigAAChangeOrNo, CoreorNo,E1orNo, HVRorNo,E2orNo,NS1orNo,NS2orNo,NS4AorNo,NS5BorNo)[,rownames(modcoef)]
        
        if (data$a[i]==1) rown =1
        if (data$t[i]==1) rown =2
        if (data$c[i]==1) rown =3
        if (data$g[i]==1) rown =4
        
        data$EstiamtedMF[i] <- exp(setUpDat[rown,] %*% coef.vals)
}


#visualize the results
plot(data$mean[500:2000], pch=21, col="#66CCEE", cex=.5, ylab='mutation frequency')
points(data$EstiamtedMF[500:2000],pch=21, col="blue",cex=0.5)
legend('topright', legend=c('Ovserved', 'Estimated'), col=c("#66CCEE",'blue'), pch= 21,cex=.7)

plot(data$mean, pch=21, col="#66CCEE", cex=.5, ylab='mutation frequency')
points(data$EstiamtedMF,pch=21, col="blue",cex=0.5)
legend('topright', legend=c('Ovserved', 'Estimated'), col=c("#66CCEE",'blue'), pch= 21,cex=.7)


# find the differences betweeen [observed - estimated]
data$diff<- data$mean-data$EstiamtedMF
plot(data$diff,pch=21, col="#66CCEE", cex=.5, ylab="Difference between Obs & Est MF",xlab="genome position")

data$over0.008<-sapply(data$diff,function(x){ if (abs(x)>0.008) x=x
                            else x<-NA})
plot(data$over0.008,pch=21, col="#66CCEE", cex=.5)

data_all<-data
data_s<-data[data$Stop==0,]
mean(data_all$diff) #Most mutation frequency is lower than estiamted frequency. 


# the differece is small, so focus on the lowest 5%

#####  Look at the results by gene
###

ByGenes<-list()
for (i in 2:12){
        df<-data[data$gene==i,]
        n=i-1
        ByGenes[[(i-1)]]<-df
        names(ByGenes)[(i-1)]<-paste0(genes$Gene[i])
}

pdf(file="Output/GLM/Compare_ObservedvsEstimatedMF.pdf", width=8,height=11.5)
par(mfrow = c(6,2),mar=c(2,2,1,1))
for (i in 1:length(ByGenes)){
        df<-ByGenes[[i]]
        plot(df$over0.008, pch=16, col="#EE6677", cex=.8, main = names(ByGenes)[i], xlab='Genome position', ylab= '[Observed - Estimated] mutation Frequency')  
        abline(h=0, col="gray")
}
dev.off()



pdf(file="Output/GLM/Compare_ObservedvsEstimatedMF_bottom10.pdf", width=8,height=11.5)
par(mfrow = c(6,2),mar=c(2,2,1,1))
for (i in 1:length(ByGenes)){
        df<-ByGenes[[i]]
        df<-df[df$Stop==0,]
        df_bottom10<-head(df[order(df$diff, decreasing= F),], n=10)
        plot(df_bottom10$pos,df_bottom10$diff, pch=16, col="#EE6677", cex=.8, main = names(ByGenes)[i], xlab='Genome position', ylab= '[Observed - Estimated] mutation Frequency')  
}
dev.off()


#################
#Attach the reference sequence info
Transitions<-read.csv("Output/Mut.freq.filtered/Summary_Ts.Q35.csv")
Transitions<-Transitions[,-1]
data<-merge(data, Transitions[,1:2], by ="pos")


#Select only the highly deviated sites
data2<-data[!is.na(data$over0.008),]#151 sites
table(data2$ref) 
# a  c  g  t 
#28 64 24 35  


#Sites with higher than expected mutation freq
df3<-data2[(!is.na(data2$over0.008)& data2$over0.008>0),] #137
#Sites with lower than expected mutation freq
df4<-data2[(!is.na(data2$over0.008)& data2$over0.008<0),] #14

table(df3$ref)
s#a  c  g  t 
#25 64 24 24 

table(df4$ref)
#a   c   g   t 
# 3  0  0 11 

## how many of the lower than expected sites are CpG creating sits?
table(df4$CpG)
#0   1 
#11  3

table(df3$CpG)
# 0  1 
#124  13

#CpG sites have lower mean frequency? -> No
mean(df4$mean[df4$CpG==0])
#[1] 0.004098213
mean(df4$mean[df4$CpG==1])
#[1] 0.004157806
mean(df3$mean[df3$CpG==0])
#[1] 0.01906609
mean(df3$mean[df3$CpG==1])
#[1] 0.01938992

mean(data2$mean[data2$CpG==0]) #[1] 0.01784648
mean(data2$mean[data2$CpG==1]) #[1] 0.0165339




###################
genes<-read.csv("Data/HCV_annotations2.csv", stringsAsFactors = F)
genes[6,1]<-"NS1(P7)"
# Add the codon position info

#proportion of deviated sites:
PercentDeviates<-data.frame("Gene"= paste0(genes$Gene[2:12]), "Proportion"=rep('',times=11), stringsAsFactors = FALSE)
ByGene2<-list()
for (i in 2:12){
        vec<-data.frame('pos'=c(genes$start[i]:((genes$end[i]-genes$start[i])+genes$start[i])))
        df<-data[data$gene==i,]
        dfm<-merge(df, vec, by="pos", all.y=TRUE)
        dfm$codon<-rep(1:3, times=nrow(dfm)/3) 
        
        dfm<-dfm[!is.na(dfm$mean),]
        #filter to the most deviated sites (>0.008)
        df2<-dfm[!is.na(dfm$over0.008),]
        
        PercentDeviates$Proportion[(i-1)] <- nrow(df2)/(nrow(dfm))*100
        PercentDeviates$DevCounts[(i-1)]<-(nrow(df2))
        PercentDeviates$Length[i-1]<-nrow(dfm)
        PercentDeviates$higher[(i-1)]<-nrow(df2[df2$over0.008>0,])
        PercentDeviates$lower[(i-1)]<-nrow(df2[df2$over0.008<0,])
        
        ByGene2[[i]]<-dfm
        names(ByGene2)[[i]]<-paste0(genes$Gene[i])
}

write.csv(PercentDeviates, "Output/GLM/Observed.vs.Estimated.MuFreq.byGenes_noStop.csv")

### Lookat the bottom 10 sites (lower than expected mutation frequency = conserved sites) -including the stop sites
Summary<-data.frame("Gene"=genes$Gene[2:12], 'a'= rep(0, times=11), 't'= rep(0, times=11),'c'= rep(0, times=11),'g'= rep(0, times=11))
for (i in 2:12){
        dat<-ByGene2[[i]]
        dat$ref<-as.character(dat$ref)
        dat10<-head(dat[order(dat$diff, decreasing= F),], n=10)
        k=i-1
        Summary$max_diff[k]<-max(dat10$diff)
        Summary$Nonsyn[k]<-nrow(dat10[dat10$Nonsyn==1,])
        Summary$Stop[k]<-nrow(dat10[dat10$Stop==1,])
        Summary$CpG[k]<-nrow(dat10[dat10$CpG==1,])
        Summary$bigAAChange[k]<-nrow(dat10[dat10$bigAAChange==1,])
        nt<-data.frame(table(dat10$ref),stringsAsFactors = F)
        for (j in c('a','t','c','g')){
                if (j %in% levels(nt$Var1)) Summary[k,j]<-nt$Freq[nt$Var1==j] 
                else Summary[k,j]<-0 
        }
                        
        #the max diff site
        Summary$lowest_pos[k]<- dat10$pos[which.min(dat10$diff)]
        Summary$lowest_base[k]<-dat10$ref[which.min(dat10$diff)]
        Summary$lowest_codon[k]<-dat10$codon[which.min(dat10$diff)]

        #codon positions
        Summary$No_codon1[k]<-nrow(dat10[dat10$codon==1,])
        Summary$No_codon2[k]<-nrow(dat10[dat10$codon==2,])
        Summary$No_codon3[k]<-nrow(dat10[dat10$codon==3,])
                                   
}

write.csv(Summary,"Output/GLM/Summary_lowest10_byGene_all.csv")


#Filter out the stop sites
Summary2<-data.frame("Gene"=genes$Gene[2:12], 'a'= rep(0, times=11), 't'= rep(0, times=11),'c'= rep(0, times=11),'g'= rep(0, times=11))
for (i in 2:12){
        dat<-ByGene2[[i]]
        dat$ref<-as.character(dat$ref)
        dat<-dat[dat$Stop==0,]
        dat10<-head(dat[order(dat$diff, decreasing= F),], n=10)
        k=i-1
        Summary2$max_diff[k]<-max(dat10$diff)
        Summary2$Nonsyn[k]<-nrow(dat10[dat10$Nonsyn==1,])
        Summary2$Stop[k]<-nrow(dat10[dat10$Stop==1,])
        Summary2$CpG[k]<-nrow(dat10[dat10$CpG==1,])
        Summary2$bigAAChange[k]<-nrow(dat10[dat10$bigAAChange==1,])
        nt<-data.frame(table(dat10$ref),stringsAsFactors = F)
        for (j in c('a','t','c','g')){
                if (j %in% levels(nt$Var1)) Summary2[k,j]<-nt$Freq[nt$Var1==j] 
                else Summary2[k,j]<-0 
        }
        
        #the max diff site
        Summary2$lowest_pos[k]<- dat10$pos[which.min(dat10$diff)]
        Summary2$lowest_base[k]<-dat10$ref[which.min(dat10$diff)]
        Summary2$lowest_codon[k]<-dat10$codon[which.min(dat10$diff)]
        
        #codon positions
        Summary2$No_codon1[k]<-nrow(dat10[dat10$codon==1,])
        Summary2$No_codon2[k]<-nrow(dat10[dat10$codon==2,])
        Summary2$No_codon3[k]<-nrow(dat10[dat10$codon==3,])
        
}

write.csv(Summary2,"Output/GLM/Summary_lowest10_byGene_all_noStop.csv")



#Look at the drug resistant sites
dr2<-read.csv("Data/HCV_drugresistance_NTpos.csv")
drpos<-data$pos %in% dr2$NTpos 
DRsites_all<-data[drpos,]
#Ave mut freq of known drug resistant sites
mean(DRsites_all$mean)# 0.004984544

mean(DRsites_all$mean[DRsites_all$Nonsyn==1])  #[1] 0.004442822
mean(DRsites_all$mean[DRsites_all$Nonsyn==0])  #0.0066699

table(DRsites_all$Nonsyn) # More Nonsyn drug resistant sites than syn DR sites
#0  1 
#9 28

######
# Look at only the highly deviated drug resistant sites: 
DRsites_h<-data2$pos%in% dr2$NTpos
# None of the sites matched to the drug resistant sites


####
## Look at the 3 genes with known drug resistance
names(ByGene2)

#add the codon positions to the dataframe
#genes$end[8]-genes$start[8]
names(ByGene2)
ByGene3<-list()
ByGene3[[1]]<-ByGene2[[8]]
ByGene3[[2]]<-ByGene2[[11]]
ByGene3[[3]]<-ByGene2[[12]]

NS3<-ByGene3[[1]]

Summary1<-data.frame("Gene"=c('NS3','NS5A','NS5B'))
for (i in 1:3){
        dat<-ByGene3[[i]]
        dat$ref<-as.character(dat$ref)
        #summary of highly deviated sites form expected mut freq
        dat<-dat[!is.na(dat$over0.008),]
        dat1<-dat[dat$diff<0,]
        dat2<-dat[dat$diff>0,]
        
        Summary1$Total_sites[i]<-nrow(dat)
        Summary1$Total_lower[i]<-nrow(dat1)
        Summary1$Total_higher[i]<-nrow(dat2)
        
        Summary1$mean_diff[i]<-mean(dat$diff)
        Summary1$meanMF_syn[i]<-mean(dat$mean[dat$Nonsyn==0])
        Summary1$menaMF_ns[i]<-mean(dat$mean[dat$Nonsyn==1])
        
        nt<-data.frame(table(dat$ref))
        Summary1$a[i]<-nt$Freq[1]
        Summary1$t[i]<-nt$Freq[4]
        Summary1$c[i]<-nt$Freq[2]
        Summary1$g[i]<-nt$Freq[3]
        
        Summary1$CpG[i]<-nrow(dat[dat$CpG==1,])
        Summary1$noCpG[i]<-nrow(dat[dat$CpG==0,])
        
        #the max diff site
        if (i==2) {
                Summary1$lowest_pos[i]<- NA
                Summary1$lowest_base[i]<-NA
                Summary1$lowest_codon[i]<-NA}
                
        else  {      
        Summary1$lowest_pos[i]<- dat1$pos[which.min(dat1$diff)]
        Summary1$lowest_base[i]<-dat1$ref[which.min(dat1$diff)]
        Summary1$lowest_codon[i]<-dat1$codon[which.min(dat1$diff)]}
        
        Summary1$highest_pos[i]<- dat2$pos[which.max(dat2$diff)]
        Summary1$highest_base[i]<-dat2$ref[which.max(dat2$diff)]
        Summary1$highest_codon[i]<-dat2$codon[which.max(dat2$diff)]
             
        #codon positions
        if (i==2) {
                Summary1$Low_codon1[i]<-NA
                Summary1$Low_codon2[i]<-NA
                Summary1$Low_codon3[i]<-NA}
        else {
                Summary1$Low_codon1[i]<-nrow(dat1[dat1$codon==1,])
                Summary1$Low_codon2[i]<-nrow(dat1[dat1$codon==2,])
                Summary1$Low_codon3[i]<-nrow(dat1[dat1$codon==3,])}
        
        Summary1$high_codon1[i]<-nrow(dat2[dat2$codon==1,])
        Summary1$high_codon2[i]<-nrow(dat2[dat2$codon==2,])
        Summary1$high_codon3[i]<-nrow(dat2[dat2$codon==3,])
        
        
        
}

write.csv(Summary1,"Output/GLM/Summary_DeviatedMF_sites_3genes.csv")



#NS3
mean(NS3$mean[NS3$Nonsyn==1]) #[1] 0.003708104
mean(NS3$mean[NS3$Nonsyn==0]) #[1] 0.0094703

#plot(NS3$diff[NS3$Nonsyn==1])

mean(NS3_com$diff,na.rm=T) #-2.507588e-05 


# The most higher than expected site:
NS3_com[which.max(NS3_com$diff),]
#     pos a t c g Syn Nonsyn Stop CpG bigAAChange       mean gene hvr EstiamtedMF       diff   over0.01 ref codon
#546 3965 0 0 0 1   1      0    0   0           0 0.03046536    8   0 0.006177991 0.02428737 0.02428737   g     3

# the most lower than expected site:
NS3_com[which.min(NS3_com$diff),]
#      pos a t c g Syn Nonsyn Stop CpG bigAAChange        mean gene hvr EstiamtedMF       diff   over0.01 ref codon
#1591 5010 0 1 0 0   0      1    0   0           0 0.001806651    8   0  0.01412755 -0.0123209 -0.0123209   t     1

## Which are highly deviated sites within NS3
NS3s<-NS3_com[!is.na(NS3_com$over0.008),] #76 sites
table(NS3s$codon)
#1  2  3 
#25 40 11

NS3s_low<-NS3s[NS3s$diff<0,] #65 sites
t1<-data.frame(table(NS3s_low$codon))   # Most are 1st and 2nd positions in the codon
#1  2  3 
#24 40  1 
table(NS3s_low$Nonsyn)  # all are nonsyn mutations
#1 
#65

table(NS3s_low$CpG)  # most are non-CpG sites
# 0  1 
#57  8 
table(NS3s_low$bigAAChange)  # most are Big AA cahnge sites
# 0  1 
#12 53 



###
#2. NS5A

vec2<-data.frame('pos'=c(genes$start[11]:((genes$end[11]-genes$start[11])+genes$start[11])))

NS5A_com<-merge(NS5A, vec2, by="pos", all.y=TRUE)
NS5A_com$codon <- rep(1:3, times=nrow(NS5A_com)/3)

mean(NS5A_com$diff,na.rm=T)  #-0.003766768
NS5A_com[which.max(NS5A_com$diff),]
#     pos a t c g Syn Nonsyn Stop CpG bigAAChange       mean gene hvr EstiamtedMF       diff   over0.01 ref codon
#927 7184 0 0 1 0   1      0    0   0           0 0.02834871   11   0 0.007553138 0.02079558 0.02079558   c     3
NS5A_com[which.min(NS5A_com$diff),]
#    pos a t c g Syn Nonsyn Stop CpG bigAAChange        mean gene hvr EstiamtedMF        diff    over0.01 ref codon
#35 6292 0 1 0 0   0      1    0   0           1 0.001878365   11   0  0.01412755 -0.01224919 -0.01224919   t     2

####################

data

