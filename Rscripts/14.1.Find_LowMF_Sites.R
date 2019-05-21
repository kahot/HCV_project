library(MASS)
library(betareg)
library(tidyverse)
library(purrr)
library(zoo)


############

Glmdata<-read.csv("Output/GLM/GlmdataFull.Ts.Q35.csv")
GlmData<-Glmdata[,-1]
s<-length(list.files("Output/Overview3/",pattern="overview3.csv"))
data<-GlmData


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


# find the differences betweeen [observed - estimated]
data$diff<- data$EstiamtedMF -data$mean
range(data$diff)  #-0.02639719  0.01081907
plot(data$diff)
#data$over0.008<-sapply(data$diff,function(x){ if (x>0.008) x=x
#                       else x<-NA})
    
write.csv(data,"Output/GLM/Ts_Obs.vs.Estimated.csv")

                        
#data_all<-data
data_s<-data[data$Stop==0,]

# the differece is small for lower end, so focus on the lowest 5%

#####Data Prep
#Attach the reference sequence info
Transitions<-read.csv("Output/Mut.freq.filtered/Summary_Ts.Q35.csv")
Transitions<-Transitions[,-1]
data<-merge(data, Transitions[,c("pos","ref")], by ="pos")

#add the codon position
vec<-data.frame('pos'= 342:8609)
dfm<-merge(data, vec, by="pos", all.y=TRUE)
dfm$codon<-rep(1:3, times=nrow(dfm)/3) 
#remove the sites with NA
dfm<-dfm[!is.na(dfm$mean),]

k=1
MFdata<-list()
for (i in c("a","t","c","g")){
        df<- dfm[dfm$ref==i,]
        for (type in c(0,1)){
                df1<-df[df$Nonsyn==type,]
                
                for (cpg in c(0,1)){
                        df2<-df1[df1$CpG==cpg,]
                        MFdata[[k]]<-df2
                        names(MFdata)[k]<-paste0(i,"_",type,cpg)
                        k=k+1
                }
        }
}
                        
                        
# data files in the order of a-t-c-g with syn-> nonsyn-> noncpg->cpg: a_00 means "syn, noncpg a"
names(MFdata)
v<-c(1:8,9,11,13,15)

##########
### Bottom5%
pdf(file="Output/GLM/Bottom5percentMF_bytype_plot.pdf", width=8,height=11.5)
par(mfrow = c(6,2),mar=c(2,2,1,1))

MFbottom5<-list()
k=1
for (i in v){
        df<- MFdata[[i]]
        no<-as.integer(nrow(df)*0.05)
        Bottom5<-head(df[order(df$mean, decreasing= F),], n=no)
        MFbottom5[[k]]<-Bottom5
        names(MFbottom5)[k]<-names(MFdata)[i]
        k=k+1
        
        plot(Bottom5$pos,Bottom5$mean, pch=16, col="#EE6677", cex=.8, main = names(MFdata)[i], 
             xlab="Genome position",
            ylab="Mutation Frequency")  
        
        
}
dev.off()        
Sum_Bottom5<-do.call(rbind,MFbottom5)      
table(Sum_Bottom5$codon)        
#  1   2   3 
#174 149  68


table(Sum_Bottom5$gene)
#  2  3  4  5  6  7  8  9 10 11 12 
# 46 21  3 33  4 21 85 11 38 61 68
# 45 21  1 24  2 21 94  8 36 61 78
LowMF_genes<-data.frame(table(Sum_Bottom5$gene))
LowMF_genes$Var1<-as.integer(as.character(LowMF_genes$Var1))
genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
LowMF_genes$gene<-genes$Gene[LowMF_genes$Var1]


table(Sum_Bottom5$ref)
#  a   c   g   t 
# 76 120 114  81


# Length of each gene portion
LowMF_genes$length<-''
LowMF_genes$percent<-''
for (i in 1:nrow(LowMF_genes)){
        n<-LowMF_genes$Var1[i]
        df<-dfm[dfm$gene==n,]
        len<-nrow(df)
        LowMF_genes$length[i]<-len
        LowMF_genes$percent[i]<-(LowMF_genes$Freq[i])/len*100       
}

# plot proportion
dat<-LowMF_genes[,c("gene","percent")]
rownames(dat)<-LowMF_genes$gene
dat$percent<-as.numeric(dat$percent)

#z=c(0.7,0.3,0.7,0.3)
plot1<-ggplot(dat,aes(x=gene,y=percent))+geom_bar(stat='identity', fill='#66CCEEE6',width=0.8)+
        labs(x="Genes",y="Proportion of highy conserved sites (%)")+
        theme_classic()
ggsave(filename=paste0("Output/SummaryFigures/Bottom5percent_MF-by-genes.pdf"),plot=plot1, width = 4.7, height = 3.5)




# look at the codon


codon<-data.frame(table(Sum_Bottom5$codon))
colnames(codon)<-c("Codon_position","Count")

plot2<-ggplot(codon,aes(x=Codon_position,y=Count))+geom_bar(stat='identity', fill='#44AA99CC',width=0.8)+
        labs(x='Codon position',y="Counts")+
        theme_classic()
ggsave(filename=paste0("Output/SummaryFigures/ConservedSites-CodonPosition.pdf"),plot=plot2, width = 3, height = 3.5)


#CpG making or not in the bottom 10 with A or T
cpg<-Sum_Bottom5[Sum_Bottom5$ref=='a'|Sum_Bottom5$ref=='t',]
table(cpg$CpG) 
# 0  1 
#101  56 

mean(cpg$mean[cpg$CpG==0])  # 0.002537268     
mean(cpg$mean[cpg$CpG==1])  # 0.002748962    


########################################
#Sel Coeff
mutrates<-read.csv("Data/Geller.mutation.rates.csv")
Sum_Bottom5$TSmutrate[Sum_Bottom5$ref=="a"]<-mutrates$mut.rate[mutrates$mutations=="AG"]
Sum_Bottom5$TSmutrate[Sum_Bottom5$ref=="c"]<-mutrates$mut.rate[mutrates$mutations=="CU"]
Sum_Bottom5$TSmutrate[Sum_Bottom5$ref=="g"]<-mutrates$mut.rate[mutrates$mutations=="GA"]
Sum_Bottom5$TSmutrate[Sum_Bottom5$ref=="t"]<-mutrates$mut.rate[mutrates$mutations=="UC"]

Sum_Bottom5$EstSC<-""
for (i in 1:nrow(Sum_Bottom5)){
        Sum_Bottom5$EstSC[i] <- EstimatedS(Sum_Bottom5$TSmutrate[i],Sum_Bottom5[i,'mean'])
}

Sum_Bottom5$EstSC<-as.numeric(Sum_Bottom5$EstSC)
plot(Sum_Bottom5$EstSC)
hist(Sum_Bottom5$EstSC)
mean(Sum_Bottom5$EstSC) #0.005626306

Sum_Bottom5[which.max(Sum_Bottom5$diff),]
#          pos a t c g Syn Nonsyn Stop CpG bigAAChange        mean gene Core E1 HVR1 E2 NS1 NS2 NS3 NS4A NS4B
#t_00.7962 8303 0 1 0 0   1      0    0   0           0 0.003255265   12    0  0    0  0   0   0   0    0    0
#          NS5A NS5B EstiamtedMF       diff ref codon TSmutrate       EstSC
#t_00.7962    0    1  0.01407433 0.01081907   t     3   1.7e-05 0.005222309




##############################################################################
#############
# The lowest MF sites are not necessarily the most costly sites?
#sc.t from 12.3 .SC site
sc$EstSC
#look at the top 3% in SC
(7958-341)*0.03   #228.5 sotes
sc250<-head(sc[order(sc$EstSC, decreasing= T),], n=250)
table(sc250$gene)
table(sc250$ref)
hist(sc250$EstSC)
mean(sc250$EstSC) #0.007648826



#1. plot by genes


# calculate % based on the length of each gene portion
SCs<-data.frame(table(sc250$gene))
SCs$Var1<-as.character(SCs$Var1)
colnames(SCs)[1]<-"Genes"
SCs$length<-''
SCs$percent<-''
for (i in 1:nrow(SCs)){
        n<-SCs$Genes[i]
        df<-dfm[dfm$gene==n,]
        len<-nrow(df)
        SCs$length[i]<-len
        SCs$percent[i]<-(SCs$Freq[i])/len*100       
}

SCs$percent<-as.numeric(SCs$percent)
plot1<-ggplot(SCs,aes(x=Var1,y=percent))+geom_bar(stat='identity', fill='#66CCEEE6',width=0.8)+
        labs(x="Genes",y="Proportion (%) of most costly sites (top 3%) ")+
        theme_classic()
ggsave(filename=paste0("Output/SummaryFigures/HighestSC-by-genes.pdf"),plot=plot1, width = 4.7, height = 3.5)


# plot by nucleotide bases
all<-data.frame(table(sc.t$ref))
highsc<-data.frame(table(sc$ref))
colnames(all)<-c("nucleotide",'count_all')
colnames(highsc)<-c("nucleotide",'count')
NT<-merge(all, highsc, by ="nucleotide")

NT$percent<-NT$count/ NT$count_all*100
NT$nucleotide<-toupper(NT$nucleotide)
plot2<-ggplot(NT,aes(x=nucleotide,y=percent))+geom_bar(stat='identity', fill='#44AA99CC',width=0.8)+
        labs(x='Nucleotide',y="Percent of most costly sites")+
        theme_classic()+ scale_y_continuous(limits=c(0,8),oob = rescale_none)
ggsave(filename=paste0("Output/SummaryFigures/HighestSC-NucPercent.pdf"),plot=plot2, width = 3, height = 3.5)




####################
#######  Check the total number of reads per site to see if low read numbers cause lower mut frequency ####
reads<-list()
for (i in 1:length(FilteredOverview)){ 
        df<-FilteredOverview[[i]]
        reads[[i]]<-df[,c("pos","TotalReads")]
        }

#assign column names for the list
for (i in 1:length(reads)) {
        colnames(reads[[i]])<-c("pos",paste0(names(FilteredOverview[i])))
}

Reads<-reads %>% purrr::reduce(full_join, by='pos') 
s<-length(reads)
Reads$Ave<-rowMeans(Reads[2:(s+1)], na.rm=T)
         
plot(Ave~pos, data=Reads)   

readave<-Reads[,c('pos','Ave')]
hist(readave$Ave, breaks=50)
lowreads.pos<-readave$pos[readave[,2]<=1000]
mean(data$mean[data$pos %in% lowreads.pos ]) #0.005523429 (<2000),  0.007405742(<1000)
mean(data$mean) # 0.005577214

#the menas for the low reads are not lowere than the overall mean.

##################



###
bottom<-data.frame(table(Sum_Bottom10$ref))
colnames(all)<-c("nucleotide",'count_all')
colnames(bottom)<-c("nucleotide",'count_bottom')
NT<-merge(all, bottom, by ="nucleotide")

NT$percent<-NT$count_bottom/ NT$count_all*100
NT$nucleotide<-toupper(NT$nucleotide)
#66CCEE","#228833","#CCBB44","#EE6677","#AA3377","#4477AA","#BBBBBB")
colors=c("#44AA99","#0077BB","#CC6677" )
plot2<-ggplot(NT,aes(x=nucleotide,y=percent))+geom_bar(stat='identity', fill='#44AA99CC',width=0.8)+
        labs(x='Nucleotide',y="% of highly conserved sites")+
        theme_classic()+ scale_y_continuous(limits=c(0,3),oob = rescale_none)
ggsave(filename=paste0("Output/SummaryFigures/ConservedSites-NucPercent.pdf"),plot=plot2, width = 3, height = 3.5)





############ Bottom 10 #########  No need - delete later

pdf(file="Output/GLM/Bottom10MF_plot.pdf", width=8,height=11.5)
par(mfrow = c(6,2),mar=c(2,2,1,1))

MFbottom10<-list()
k=1
for (i in v){
        df<- MFdata[[i]]
        Bottom10<-head(df[order(df$mean, decreasing= F),], n=10)
        MFbottom10[[k]]<-Bottom10
        names(MFbottom10)[k]<-names(MFdata)[i]
        k=k+1
        
        plot(Bottom10$pos,Bottom10$diff, pch=16, col="#EE6677", cex=.8, main = names(MFdata)[i], 
             xlab='Genome position', 
             ylab= '[Observed - Estimated] mutation Frequency')  
        
        
}
dev.off()        
Sum_Bottom10<-do.call(rbind,MFbottom10)      
table(Sum_Bottom10$codon)        
# 1  2  3 
#37 33 50

table(Sum_Bottom10$gene)
# 2  3  5  6  7  8  9 10 11 12 
#21  1 12  2  9 17  4 10 16 28 

LowMF_genes<-data.frame(table(Sum_Bottom10$gene))
LowMF_genes$Var1<-as.integer(as.character(LowMF_genes$Var1))
genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
LowMF_genes$gene<-genes$Gene[LowMF_genes$Var1]###############

