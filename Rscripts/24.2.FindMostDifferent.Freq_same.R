# explore the sites that differ most in mutation frequency between the two genotypes

library(stringr)
library(tidyverse)
library(zoo)
library(reshape2)
library(ggplot2)
source("Rscripts/baseRscript.R")


cols2.60<-paste0(cols2,"99")
genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)

cols3<-c("#66CCEE4D","#66CCEE","#EE66774D","#EE6677")


#############
# find the sites that are most different in mutation frequencies between the genotypes

mf<-read.csv("Output_all/Unfiltered/Ts.Same_Mean_3genotypes.csv", stringsAsFactors = F )
mf<-mf[-c(1:3),-1]
mf<-mf[mf$merged.pos<=8619,]

geno<-c("1A","1B","3A")
CombG<-t(combn(geno,2))
s<-as.integer(0.05*nrow(mf)) #top 5%

for (g in 1:3){
        g1<-CombG[g,1]
        g2<-CombG[g,2]
        compare<-paste0(g1," vs. ",g2)
        print(compare)
        Freq.by.gene<-data.frame(Gene=genes$Gene[1:12])
        mf1<-mf
        mf1$diff<-mf1[,paste0("mean.",g1)]-mf1[,paste0("mean.",g2)]
        pdf(paste0("Output_all/Difference/Same/Hist.diff.in.mutfreq.",compare,".pdf"),width=6,height=5)
        hist(mf1$diff, breaks=40, xlab="Difference in mutation frequency", main=compare)
        dev.off()
        
        #Order the sites with Genotyepe g1 having higher/lower mut freq than g2 (mvf.H / mvf.L)
        mf.H<-mf1[order(mf1$diff, decreasing = T,na.last = T),]
        mf.L<-mf1[order(mf1$diff,decreasing=F,na.last=T),]
        #select top 5% of the sites
        mf.H<-mf.H[c(1:s),] 
        mf.L<-mf.L[c(1:s),] 
        
        h1<-as.data.frame(table(mf.H$gene))
        h2<-as.data.frame(table(mf.L$gene))
        colnames(h1)<-c("Gene",paste0("higher in ",g1))
        colnames(h2)<-c("Gene",paste0("higher in ",g2))
        Freq.by.gene<-merge(Freq.by.gene,h1, by="Gene",all.x=T)
        Freq.by.gene<-merge(Freq.by.gene,h2, by="Gene",all.x=T)
        write.csv(Freq.by.gene, paste0("Output_all/Difference/Same/Diff.mf.by.gene.",compare,".csv"))
        
        Freq2<-melt(Freq.by.gene, id.vars="Gene")
        colnames(Freq2)[2:3]<-c("Genotype","Counts")
        
        ggplot(Freq2,aes(x=Gene, y=Counts, fill=Genotype))+
                geom_bar(stat='identity', position='dodge')+
                ggtitle("Sites with the largest difference in mutation frequency between the genotypes")+
                theme(legend.title=element_blank())+
                theme(plot.title = element_text(size =11))+
                ylab("Number of sites in top 5%")
        
        ggsave(filename=paste0("Output_all/Difference/Same/Diff.mf.by.gene.",compare,".pdf"),width = 8, height = 6)
        
        #summary tables for g1
        #removed non protein coding sections
        mf.H<-mf.H[mf.H$gene!="5' UTR",]
        mf.L<-mf.L[mf.L$gene!="5' UTR",]
        
        t1<-table(mf.H[,paste0("ref.",g1)],mf.H[,paste0("Type.",g1)])
        t2<-table(mf.H$codon,mf.H[,paste0("Type.",g1)])  
        sum1<-as.data.frame(rbind(t1,t2))
        write.csv(sum1, paste0("Output_all/Difference/Same/Summary.stats.higher.in.",g1,".csv"))
        
        t11<-table(mf.L[,paste0("ref.",g2)],mf.L[,paste0("Type.",g2)])
        t22<-table(mf.L$codon,mf.L[,paste0("Type.",g2)])  
        sum2<-as.data.frame(rbind(t11,t22))
        write.csv(sum2, paste0("Output_all/Difference/Same/Summary.stats.higher.in.",g2,".csv"))
        
        #syn vs. nonsyn proportion
        sink(paste0("Output_all/Difference/Same/syn-nonsyn-ratio.",compare,".txt"))
        print(paste0("Sites with MF higher in ",g1))
        print(margin.table(t2,2)/margin.table(t2))
        print(paste0("Sites with MF higher in ",g2))
        print(margin.table(t22,2)/margin.table(t22))
        sink(NULL)
        
        ## absolute top 5% different sites
        mf1$diff.abs<-abs(mf1$diff)
        mf.top<-mf1[order(mf1$diff.abs, decreasing = T,na.last = T),]
        mf.top<-mf.top[mf.top$gene!="5' UTR",]
        mf.top<-mf.top[c(1:s),]
        write.csv(mf.top,paste0("Output_all/Difference/Same/Top5percent.MostDifferent.MF.", compare,".csv"))
        
        #ta1<-table(mf.top$gene)
        #ta2<-table(mf.top[,paste0("Type.",g1)])
        #ta3<-table(mf.top[,paste0("Type.",g2)])
        #ta4<-table(mf.top$codon)
        
        #proportion of most different sites in mut freq
        PercentHighDiff<-data.frame("Gene"= paste0(genes$Gene[2:12]), "Proportion"=rep('',times=11), stringsAsFactors = FALSE)
        Dat<-mf1[,c("merged.pos", "gene","diff.abs")]
        Dat<-Dat[!is.na(Dat$diff.abs),]
        
        for (i in 2:12){
                n<-Dat[Dat$gene==genes$Gene[i],]
                df<-mf.top[mf.top$gene==genes$Gene[i],]
                
                PercentHighDiff$Proportion[i-1] <- nrow(df)/(nrow(n))*100
                PercentHighDiff$Counts[i-1]<-nrow(df)
                PercentHighDiff$No.Sites[i-1]<-nrow(n)
        }
        
        PercentHighDiff$Proportion<-as.numeric(PercentHighDiff$Proportion)
        write.csv(PercentHighDiff,paste0("Output_all/Difference/Same/Top5percent.MostDifferent.SummarybyGene.",compare,".csv"))
        
        ggplot(PercentHighDiff,aes(x=Gene,y=Proportion))+geom_bar(stat='identity', fill='#66CCEECC',width=0.8)+
                labs(x="Genes",y=paste0("% of sites highly different in MF betw/ ",g1," & ",g2))+
                theme_classic()
        ggsave(filename=paste0("Output_all/Difference/Same/Top5percent.MostDifferent.MF.by.gene.",compare,".pdf"), width = 5.7, height = 4.5)
        
        
        
        #absolute difference top 1%
        s1<-as.integer(0.01*nrow(mf1))
        mf.top2<-mf.top[c(1:s1),]
        
        write.csv(mf.top2,paste0("Output_all/Difference/Same/Top1percent.MostDifferent.MF.",compare,".csv"))
        
        PercentHighDiff2<-data.frame("Gene"= paste0(genes$Gene[2:12]), "Proportion"=rep('',times=11), stringsAsFactors = FALSE)
        for (i in 2:12){
                n<-Dat[Dat$gene==genes$Gene[i],]
                df<-mf.top2[mf.top2$gene==genes$Gene[i],]
                
                PercentHighDiff2$Proportion[i-1] <- nrow(df)/(nrow(n))*100
                PercentHighDiff2$Counts[i-1]<-nrow(df)
                PercentHighDiff2$No.Sites[i-1]<-nrow(n)
        }
        PercentHighDiff2$Proportion<-as.numeric( PercentHighDiff2$Proportion)
        write.csv(PercentHighDiff2,paste0("Output_all/Difference/Same/Top1percent.MostDifferent.SummarybyGene.",compare,".csv"))
        
        ggplot(PercentHighDiff2,aes(x=Gene,y=Proportion))+geom_bar(stat='identity', fill='#EE667799',width=0.8)+
                labs(x="Genes",y=paste0("% of sites highly different in MF betw/ ",g1," & ",g2))+
                theme_classic()
        ggsave(filename=paste0("Output_all/Difference/Same/Top1per.MostDifferent.MF.by.gene.",compare,".pdf"),width = 5.7, height = 4.5)

        #Which Amino acids?
        for (j in 1:2){
                if (j==1) G<-g1; DF<-mf.H
                if (j==2) G<-g2;  DF<-mf.L
                aa1<-table(DF[,paste0("WTAA.",G)],DF[,paste0("Type.",G)])
                AAtotal<-table(mf1[mf1$gene!="5' UTR",paste0("WTAA.",G)],mf1[mf1$gene!="5' UTR",paste0("Type.",G)] )
                AA1<-data.frame(AA=rownames(AAtotal))
                if (length(aa1[,1])!=20){
                        AA.df<-data.frame(AA=rownames(aa1))
                        AA.df$nonsyn<-aa1[,1]
                        AA.df$syn<-aa1[,2]
                        AA1<-merge(AA1, AA.df, by="AA",all.x=T)
                        AA1[is.na(AA1)]<-0
                }
                
                else {  AA1$nonsyn<-aa1[,1]
                        AA1$syn<-aa1[,2]}
                AA1$total.per.nonsyn<-AAtotal[,1]/sum(AAtotal[,1])*100
                AA1$total.per.syn<-AAtotal[,3]/sum(AAtotal[,3])*100
                AA1$per.nonsyn<-AA1$nonsyn/sum(AA1$nonsyn)*100
                AA1$per.syn<-AA1$syn/sum(AA1$syn)*100
                AA1.1<-AA1[,c("AA","total.per.nonsyn","per.nonsyn","total.per.syn","per.syn")]
                colnames(AA1.1)<-c("AA","% total nonsyn", "% nonsyn", "% total syn", "% syn")
                AA1.1<-melt(AA1.1)
                colnames(AA1.1)[2:3]<-c("Group","Proportion")
                
                #ggplot(AA1.1, aes(x=AA, y=Proportion, fill=Group))+geom_bar(stat='identity', position='dodge')+
                #        labs(x="Amino acids", y="% of represented AA")+
                #        ggtitle(paste0("% of AA represented in highly different MF sites of ",G," in ", compare))+
                #        scale_fill_manual("", values = cols3)+theme_classic()+
                #        theme(plot.title = element_text(hjust = 0.5, size=12))
                #ggsave(filename=paste0("Output_all/Difference/Same/Which.AA.in.HigherMFof",G,".in.", compare,".pdf"),width =9, height = 4.5)
                #nonsyn only
                #ggplot(AA1.1[AA1.1$Group=="% total nonsyn"|AA1.1$Group=="% nonsyn",], aes(x=AA, y=Proportion, fill=Group))+geom_bar(stat='identity', position='dodge')+
                #        labs(x="Amino acids", y="% of represented AA")+
                #        ggtitle(paste0("% of AA represented in highly different MF sites of ",G," in ", compare))+
                #        scale_fill_manual("", values = cols3)+theme_classic()+
                #        theme(plot.title = element_text(hjust = 0.5, size=12))
                #ggsave(filename=paste0("Output_all/Difference/Same/Which.AA.in.HigherMFof",G,".in.", compare,".Nonsyn.pdf"),width = 7, height = 4.5)
#
                #show overrepresentation only
                AA1$diff.nonsyn<-(AA1$per.nonsyn-AA1$total.per.nonsyn)/AA1$total.per.nonsyn*100
                ggplot(AA1,aes(x=AA, y=diff.nonsyn))+geom_bar(stat='identity',fill='#66CCEE')+
                        theme_classic()+labs(x="Amino acids", y="% difference")+
                        ggtitle(paste0("Over-/Under-represented AAs of highly different MF sites of ", G, " (",compare,")"))+
                        theme(plot.title = element_text(size=12))
                ggsave(filename=paste0("Output_all/Difference/Same/Overrepresented.AAs.",G,".in.",compare,".pdf"), width = 7, height = 4.5)
                }      

}




###############
## calculate the ratio of nonsyn mut freq to syn mut.freq. 

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

ggplot(SNratio2,aes(x=Gene,y=ratio,group=Genotype, color=Genotype))+geom_point(position=position_dodge(width=0.4))+scale_color_manual(values=cols2)+
        theme_light()+
        theme(axis.title.x=element_blank())+ylab("Ratio of nonsyn to syn mutation frequency")
ggsave(filename="Output_all/Unfiltered/NSratio.by.gene.by.genotype.pdf",width = 9, height = 8)


############################

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



