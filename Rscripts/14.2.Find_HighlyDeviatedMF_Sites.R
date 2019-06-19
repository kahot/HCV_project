library(ggplot2)
library(miscTools)
library(tidyverse)
library(janitor)
library(reshape2)


cols2<-c("#66CCEE","#EE6677" ,"#228833")
####
HCVFiles_overview3<-list.files("Output1A/Overview3/",pattern="overview3.csv")

data<-read.csv("Output1A/GLM/GlmdataFull.Ts.Q35.csv")
data<-data[,-1]
s<-length(HCVFiles_overview3)


# from 10.3_Filtered OveriewData_GLM.R
#Read data file
GlmData <-read.csv("Output1A/GLM/GlmdataFull.Ts.Q35.csv")
GlmData1<-GlmData[,-1]

modcoef<-read.csv("Output1A/GLM/BetaReg_mod.g2_Ts.Q35.csv",stringsAsFactors = F)
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
plot(data$mean[500:2000], pch=16, col="#66CCEE", cex=.5, xaxt='n', ylab='Mutation frequency',xlab='Genome position')
points(data$EstiamtedMF[500:2000],pch=16, col="blue",cex=0.5)
xticks<-seq(0,1500, by=500)
axis(side=1, at=xticks,labels=seq(500,2000, by = 500))
legend('topright', legend=c('Ovserved', 'Estimated'), col=c("#66CCEE",'blue'), pch= 16,cex=.7)

pdf("Output1A/GLM/Expected.vs.observed.mf.pdf", width=10, height = 6)
plot(data$mean, pch=21, col="#66CCEE", cex=.3, ylab='Mutation frequency',xlab='Genome position')
points(data$EstiamtedMF,pch=21, col="blue",cex=0.3)
legend('topright', legend=c('Ovserved', 'Estimated'), col=c("#66CCEE",'blue'), pch= 21,cex=.7)
dev.off()



# find the differences betweeen [observed - estimated]
data$diff<- data$mean-data$EstiamtedMF

pdf("Output1A/GLM/Expected.vs.observed.mf-Difference.pdf", width=10, height = 6)
plot(data$diff,pch=21, col="#66CCEE", cex=.3, ylab="Difference between Obs & Est MF",xlab="Genome position")
dev.off()

#Find the sites that expected and ovserved MF differed > 0.008
data$over0.008<-sapply(data$diff,function(x){ if (abs(x)>0.008) x=x
                            else x<-NA})
plot(data$over0.008,pch=16, col="#66CCEE", cex=.5)



#Average difference in Exp-Obs Mut Freq.

MFdiff<-data.frame("Mutation Type"= c("Overall","Without Stop", "Stop", "Syn","Nonsyn"))
MFdiff$Mean.diff<-c(mean(data$diff),mean(data$diff[data$Stop==0]),mean(data$diff[data$Stop==1]),mean(data$diff[data$Syn==1]),mean(data$diff[data$Nonsyn==1]))
MFdiff$Mutation.Type<-factor(MFdiff$Mutation.Type, c("Overall","Without Stop", "Stop", "Syn","Nonsyn"))
write.csv(MFdiff, "Output1A/GLM/Expected.vs.observed.mf-Difference.csv")

ggplot(MFdiff,aes(x=Mutation.Type,y=Mean.diff))+
        geom_bar(stat="identity", color="#EE6677",fill="#EE667766",width=0.8)+
        theme_bw()+labs(x="Mutation type",y="Average difference between observed & expected MF")+
        theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
        geom_hline(yintercept=0, size=0.4, color="gray10")
ggsave("Output1A/GLM/Mean.diff.Expected.vs.observed.pdf", width=6, height=5)


#####  Look at the results by gene
genes<-read.csv("Data/HCV_annotations2.csv", stringsAsFactors = F)

ByGenes<-list()
for (i in 2:12){
        df<-data[data$gene==i,]
        ByGenes[[(i-1)]]<-df
        names(ByGenes)[(i-1)]<-paste0(genes$Gene[i])
}

pdf(file="Output1A/GLM/Compare_ObservedvsEstimatedMF.pdf", width=8,height=11.5)
par(mfrow = c(6,2),mar=c(2,2,1,1))
for (i in 1:length(ByGenes)){
        df<-ByGenes[[i]]
        plot(df$pos,df$over0.008, pch='', col=cols2[1], cex=.8, main = names(ByGenes)[i], xlab='Genome position', ylab= '(Observed - Estimated) mutation Frequency')  
        points(df$pos[df$over0.008>0],df$over0.008[df$over0.008>0], pch=16, col=cols2[1], cex=.8)  
        points(df$pos[df$over0.008<0],df$over0.008[df$over0.008<0], pch=16, col=cols2[2], cex=.8)  
        abline(h=0, col="gray")
}
dev.off()


#Find the sites with highly conserved (observed-expected to be bottom 5 %) 
l<-nrow(data)
data2<-head(data[order(data$diff, decreasing= F),], n=round(l*0.05))


ByGenes2<-list()
for (i in 2:12){
        df<-data2[data2$gene==i,]
        df<-df[(order(df$pos)),]
        ByGenes2[[(i-1)]]<-df
        names(ByGenes2)[(i-1)]<-paste0(genes$Gene[i])
}

pdf(file="Output1A/GLM/Compare_ObservedvsEstimatedMF.bottom5per.pdf", width=8,height=11.5)
par(mfrow = c(6,2),mar=c(2,2,1,1))
for (i in 1:length(ByGenes2)){
        df<-ByGenes2[[i]]
        plot(df$pos, df$diff, pch=16, col=cols2[2], cex=.8, main = names(ByGenes)[i], xlab='Genome position', ylab= '(Observed - Estimated) mutation Frequency')  
        abline(h=0, col="gray")
}
dev.off()







#################
#Attach the reference sequence/AA info
Transitions<-read.csv("Output1A/Mut.freq.filtered/Summary_Ts.Q35.csv")
data<-merge(data, Transitions[,c("pos","ref","WTAA")], by ="pos")
#most conserved 5% sites (lower than expected mf sites)
data2<-head(data[order(data$diff, decreasing= F),], n=round(l*0.05))

#which nucleotides are over/under- represented in the conserved sites.
Sum<-list()
    Sum[[1]]<-table(data2$ref)
    Sum[[2]]<-table(data$ref)
    Sum[[3]]<-Sum[[1]]/sum(Sum[[1]])
    Sum[[4]]<-Sum[[2]]/sum(Sum[[2]])
    Sum[[5]]<-Sum[[3]]-Sum[[4]]
NT<-as.data.frame(Sum[[5]]/Sum[[4]]*100)


ggplot(NT,aes(x=Var1,y=Freq))+
        geom_bar(stat="identity", color="#EE6677",fill="#EE667766",width=0.8)+
        theme_bw()+labs(x="")+
        theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
        geom_hline(yintercept=0, size=0.4, color="gray10")+
        scale_x_discrete(labels=c("A","C","G","T"))+
        theme(axis.text.x = element_text(size=14),
              axis.text.y = element_text(size=14),axis.title.y = element_text(size = 14))+
        ylab(expression(paste("Over/under-represented nt (%) \n at conserved sites")))+
        theme(plot.margin=unit(c(5,5,5,20),"points"))
ggsave("Output1A/GLM/Overrepresented.Nucleotides.in.Bottom5per.pdf", width=5, height=4)



table(data2$CpG, data2$ref)
#      a   c   g   t
#0  33 175  70  62
#1  20   0   0  38

table(data2$codon)


#proportion of over/under represented sites by gene:
#Adjust the length of the last gene

genes2<-genes[-13,]
genes2$end[12]<-data$pos[nrow(data)]
PercentDeviates<-data.frame("Gene"= paste0(genes$Gene[2:12]), "counts"=rep('',times=11), stringsAsFactors = FALSE)
codonsummary<-data.frame("Gene"= paste0(genes$Gene[2:12]),"first"=rep(0,times=11),stringsAsFactors = FALSE)
AAsummary<-list()
for (i in 2:12){
        vec<-data.frame('pos'=c(genes$start[i]:((genes$end[i]-genes$start[i])+genes$start[i])))
        df<-data2[data2$gene==i,]
        dfm<-merge(df, vec, by="pos", all.y=TRUE)
        dfm$codon<-rep(1:3, times=nrow(dfm)/3) 
        #df<-dfm[!is.na(dfm$mean),]
        PercentDeviates$counts[(i-1)]<-nrow(df)
        PercentDeviates$length[i-1]<-nrow(dfm)
        PercentDeviates$Proportion[(i-1)] <- ((nrow(df)/nrow(data2))/(nrow(dfm)/nrow(data))-1 )*100
        codonsummary$first[i-1] <-sum(dfm$codon[!is.na(dfm$mean)& dfm$codon==1])
        codonsummary$second[i-1]<-sum(dfm$codon[!is.na(dfm$mean)& dfm$codon==2])
        codonsummary$third [i-1]<-sum(dfm$codon[!is.na(dfm$mean)& dfm$codon==3])
        aa<-as.data.frame(table(df$WTAA))
        AAsummary[[i-1]]<-aa
        names(AAsummary)[i-1]<-genes$Gene[i]
        
}

write.csv(PercentDeviates, "Output1A/GLM/Overrepresented.site.proportion.by.gene.csv")

for (i in 1:length(AAsummary)) colnames(AAsummary[[i]])<-c("Gene",paste0(names(AAsummary[i])))
AA<-AAsummary%>% purrr::reduce(full_join, by='Gene')
write.csv(AA,"Output1A/GLM/Overrepresented.site.AA.by.gene.csv")

codonsummary<-codonsummary %>% adorn_totals("row")
write.csv(codonsummary, "Output1A/GLM/Overrepresented.site.codons.csv")


ggplot(PercentDeviates,aes(x=Gene,y=Proportion))+
        geom_bar(stat="identity", color="#EE6677",fill="#EE667766",width=0.8)+
        theme_bw()+labs(x="")+
        theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
        geom_hline(yintercept=0, size=0.4, color="gray10")+
        theme(axis.text.x = element_text(size=12),
              axis.text.y = element_text(size=12),axis.title.y = element_text(size = 14))+
        ylab(expression(paste("Over/under-represented % of conserved sites")))
ggsave("Output1A/GLM/Overrepresented.sites.perGene.Bottom5per.pdf", width=8, height=6)


##### Over or under represented amino acids

AA$Total<-rowSums(AA[2:12])
AAwhole<-as.data.frame(table(data$WTAA))
AAwhole$Prop.total<-AAwhole$Freq/sum(AAwhole$Freq)
AAwhole$Prop.conserved<-AA$Total/sum(AA$Total)
AAwhole$proportion<-(AAwhole$Prop.conserved-AAwhole$Prop.total)/AAwhole$Prop.total*100

ggplot(AAwhole,aes(x=Var1,y=proportion))+
        geom_bar(stat="identity",width=0.8,color="#EE6677",fill="#EE667766")+
        theme_bw()+labs(x="")+
        theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
        geom_hline(yintercept=0, size=0.4, color="gray10")+
        theme(axis.text.x = element_text(size=12),
              axis.text.y = element_text(size=12),axis.title.y = element_text(size = 14))+
        ylab(expression(paste("Over/under-represented AA \n at conserved sites (%)")))+
        theme(plot.margin=unit(c(5,5,5,20),"points"))
ggsave("Output1A/GLM/Overrepresented.AA.perGene.Bottom5per.pdf", width=10, height=6)












########################
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

write.csv(Summary2,"Output1A/GLM/Summary_lowest10_byGene_all_noStop.csv")



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

write.csv(Summary1,"Output1A/GLM/Summary_DeviatedMF_sites_3genes.csv")



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

