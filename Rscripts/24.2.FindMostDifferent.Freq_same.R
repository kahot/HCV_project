library(stringr)
library(tidyverse)
library(zoo)
library(reshape2)

cols2<-c("#66CCEE","#EE6677" ,"#228833")
colors<-c("#66CCEE99","#EE667799" ,"#22883399")



#############
# find the sites that are most different in mutation frequencies between the genomes

mf<-read.csv("Output_all/Unfiltered/Ts.Same_Mean_3genotypes.csv", stringsAsFactors = F )
mf<-mf[-c(1:3),-1]
mf<-mf[mf$merged.pos<=8619,]
mf$diff.1A.3A<-mf$mean.1A-mf$mean.3A

hist(mf$diff.1A.3A)

#find the sites with 3A having higher mut greq (mvf.H) and 1A having the higher mut. freq (mvf.L)
mf.H<-mf[order(mf$diff.1A.3A, decreasing = T,na.last = T),]
mf.L<-mf[order(mf$diff.1A.3A,decreasing=F,na.last=T),]

#most and leat different 3% sites
s<-as.integer(0.05*nrow(mf))

mf.H<-mf.H[c(1:s),] #3A has higher mf
mf.L<-mf.L[c(1:s),]  #1A has higher mf


table(mf.H$gene)
#Core      E1      E2    HVR1 NS1(P7)     NS2     NS3    NS4A    NS4B    NS5A    NS5B 
#  10      30      62       5       7      47      98      10      42      67      39  

table(mf.L$gene)
#Core      E1      E2    HVR1 NS1(P7)     NS2     NS3    NS4A    NS4B    NS5A    NS5B 
#  16      35      52       3      12      31      91      10      40      86      41  

table(mf.H$ref.3A)
# a   c   g   t 
#78 155 113  71

higher3A<-table(mf.H$ref.3A,mf.H$Type.3A)
#   nonsyn stop syn
#a     28    0  50
#c     27    1 127
#g     27    2  84
#t     13    0  58
table(mf.H$codon,mf.H$Type.3A)
# 1   2   3 
#67  27 323

table(mf.H$Type.3A)/(95+3+319)
#     nonsyn        stop         syn 
#0.227817746 0.007194245 0.764988010


t1<-table(mf$codon, mf$Type.3A)
margin.table(t1,2)/margin.table(t1)
#    nonsyn       stop        syn 
#0.63938159 0.02684564 0.33377277

###############
## calculate ratio of nonsyn mut freq to syn mut.freq. 

genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)

geno<-c("1A","1B","3A")
ratio<-data.frame(genotype=geno)
mfsummary<-list()
for (g in 1:3){
        dat<-aggregate(mf[,c(paste0("mean.",geno[g]))], list(mf[,paste0("Type.",geno[g])]), mean, na.rm=T)
        colnames(dat)[2]<-geno[g]
        mfsummary[[g]]<-dat
        ratio$ratio[g]<-dat[1,2]/dat[3,2]
        
}
#overall SN Ratio
write.csv(ratio,"Output_all/Unfiltered/NStoSynRatio.csv")

SNratio<-data.frame(Gene=c(genes$Gene[2:12]), Geno.1A=NA,Geno.1B=NA,Geno.3A=NA)
for (g in 1:3){
        for (i in 2:12 ){
                df<-mf[mf$gene==genes$Gene[i], ]
                dat<-aggregate(df[,c(paste0("mean.",geno[g]))], list(df[,paste0("Type.",geno[g])]), mean, na.rm=T)
                #colnames(dat)[2]<-geno[g]
                #mfsummary[[g]]<-dat
                SNratio[i-1,g+1]<-dat[1,2]/dat[3,2]
        }
}

write.csv(SNratio,"Output_all/Unfiltered/NStoSynRatio.by.gene.csv")

colnames(SNratio)[2:4]<-c("1A","1B","3A")
SNratio2<-melt(SNratio)
colnames(SNratio2)[2:3]<-c("Genotype","ratio")

R1<-ggplot(SNratio2,aes(x=Gene,y=ratio,group=Genotype, color=Genotype))+geom_point(position=position_dodge(width=0.3))+scale_color_manual(values=cols2)+
        theme_light()+
        theme(axis.title.x=element_blank())+ylab("Ratio of nonsyn to syn mutation frequency")
R1
ggsave(filename="Output_all/Unfiltered/NSratio.by.gene.by.genotype.pdf",plot=R1, width = 9, height = 8)




t2<-aggregate(mf[,c("mean.3A")],by=list(mf$Type.3A), FUN=mean, na.rm=T)
t2$x[1]/t2$x[3]


margin.table(t1,2)/margin.table(t1)


prop.table(t1,2)


table(mf.H$ref.1A)

table(mf.L$ref.3A)


higher3A<-table(mf.H$ref.3A,mf.H$Type.3A)


#absolute difference top 5%
mf$diff.1a.3a.abs<-abs(mf$diff.1A.3A)
mf.top<-mf[order(mf$diff.1a.3a.abs, decreasing = T,na.last = T),]
mf.top<-mf.top[c(1:s),]

write.csv(mf.top,"Output_all/Unfiltered/Top.Difference.sameasRef.1A-3A.csv")

table(mf.top$gene)
#   Core      E1      E2    HVR1 NS1(P7)     NS2     NS3    NS4A    NS4B    NS5A    NS5B 
#     12      30      59       5      12      34      92      10      41      77      45 

table(mf.top$Type.3A)
#nonsyn   stop    syn 
#    94      3    320

table(mf.top$Type.1A)
#nonsyn   stop    syn 
#    86      2    329 

table(mf.top$codon)
#1   2   3 
#75  30 312


#absolute difference top 1%
s1<-as.integer(0.01*nrow(mf))
mf.top2<-mf.top[c(1:s1),]

write.csv(mf.top,"Output_all/Unfiltered/Top1percent.Difference.sameasRef.1A-3A.csv")


###############
#proportion of most different sites in mut freq
genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)

PercentHighDiff2<-data.frame("Gene"= paste0(genes$Gene[2:12]), "Proportion"=rep('',times=11), stringsAsFactors = FALSE)
DF<-mf[,c("merged.pos", "gene","diff.1a.3a.abs")]
DF<-DF[!is.na(DF$diff.1a.3a.abs),]

for (i in 2:12){
                n<-DF[DF$gene==genes$Gene[i],]
                df<-mf.top[mf.top$gene==genes$Gene[i],]
                
                PercentHighDiff2$Proportion[i-1] <- nrow(df)/(nrow(n))*100
                PercentHighDiff2$Counts[i-1]<-nrow(df)
                PercentHighDiff2$No.Sites[i-1]<-nrow(n)
}

PercentHighDiff2$Proportion<-as.numeric( PercentHighDiff2$Proportion)
p2<-ggplot(PercentHighDiff2,aes(x=Gene,y=Proportion))+geom_bar(stat='identity', fill='#EE667799',width=0.8)+
        labs(x="Genes",y="% of sites highly different in MF betw/ 1A & 3A")+
        theme_classic()
p2
ggsave(filename=paste0("Output_all/Unfiltered/Top5per.MostDifferent.in.MF.by.genes.same.1A-3A.pdf"),plot=p2, width = 5.7, height = 4.5)



### top 1%
PercentHighDiff3<-data.frame("Gene"= paste0(genes$Gene[2:12]), "Proportion"=rep('',times=11), stringsAsFactors = FALSE)
DF<-mf[,c("merged.pos", "gene","diff.1a.3a.abs")]
DF<-DF[!is.na(DF$diff.1a.3a.abs),]

for (i in 2:12){
        n<-DF[DF$gene==genes$Gene[i],]
        df<-mf.top2[mf.top2$gene==genes$Gene[i],]
        
        PercentHighDiff3$Proportion[i-1] <- nrow(df)/(nrow(n))*100
        PercentHighDiff3$Counts[i-1]<-nrow(df)
        PercentHighDiff3$No.Sites[i-1]<-nrow(n)
}



PercentHighDiff3$Proportion<-as.numeric( PercentHighDiff3$Proportion)
p3<-ggplot(PercentHighDiff3,aes(x=Gene,y=Proportion))+geom_bar(stat='identity', fill='#EE667799',width=0.8)+
        labs(x="Genes",y="% of sites highly different in MF betw/ 1A & 3A")+
        theme_classic()
p3
ggsave(filename=paste0("Output_all/Unfiltered/Top1per.MostDifferent.in.MF.by.genes.same.1A-3A.pdf"),plot=p3, width = 5.7, height = 4.5)



table(mf.top2$ref.1A)
# a  c  g  t 
#16 32 12 23 
table(mf.top2$ref.3A)
# a  c  g  t 
#20 26 17 20 



