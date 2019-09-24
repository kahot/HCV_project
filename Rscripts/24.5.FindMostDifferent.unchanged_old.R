# explore the sites that differ most in mutation frequency between the two genotypes

library(stringr)
library(tidyverse)
library(zoo)
library(reshape2)
library(ggplot2)
source("Rscripts/baseRscript.R")

#dir.create("Output_all/Difference/Conserved/")

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

for (g in 1:3){
        g1<-CombG[g,1]
        g2<-CombG[g,2]
        compare<-paste0(g1," vs. ",g2)
        print(compare)
        
        Freq.by.gene<-data.frame(Gene=genes$Gene[1:12])
        mf1<-mf
        #filter the sites that have the same nucleotide between the genotypes (refs):
        n1<-which(colnames(mf1)==paste0("ref.",g1))
        n2<-which(colnames(mf1)==paste0("ref.",g2))
        mf1<-mf1[mf1[n1]==mf[n2],]
        s<-as.integer(0.05*nrow(mf1)) #top 5%
        
        mf1$diff<-mf1[,paste0("mean.",g1)]-mf1[,paste0("mean.",g2)]
        pdf(paste0("Output_all/Difference/Conserved/Hist.diff.in.mutfreq.",compare,".pdf"),width=6,height=5)
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
        write.csv(Freq.by.gene, paste0("Output_all/Difference/Conserved/Diff.mf.by.gene.",compare,".csv"))
        
        Freq2<-melt(Freq.by.gene, id.vars="Gene")
        colnames(Freq2)[2:3]<-c("Genotype","Counts")
        
        ggplot(Freq2,aes(x=Gene, y=Counts, fill=Genotype))+
                geom_bar(stat='identity', position='dodge')+
                ggtitle("Sites with the largest difference in mutation frequency between the genotypes")+
                theme(legend.title=element_blank())+
                theme(plot.title = element_text(size =11))+
                ylab("Number of sites in top 5%")
        
        ggsave(filename=paste0("Output_all/Difference/Conserved/Diff.mf.by.gene.",compare,".pdf"),width = 8, height = 6)
        
        
        
        
        #summary tables for g1
        #removed non protein coding sections
        mf.H<-mf.H[mf.H$gene!="5' UTR",]
        mf.L<-mf.L[mf.L$gene!="5' UTR",]
        
        t1<-table(mf.H[,paste0("ref.",g1)],mf.H[,paste0("Type.",g1)])
        t2<-table(mf.H$codon,mf.H[,paste0("Type.",g1)])  
        sum1<-as.data.frame(rbind(t1,t2))
        write.csv(sum1, paste0("Output_all/Difference/Conserved/Summary.stats.higher.in.",g1,".csv"))
        
        t11<-table(mf.L[,paste0("ref.",g2)],mf.L[,paste0("Type.",g2)])
        t22<-table(mf.L$codon,mf.L[,paste0("Type.",g2)])  
        sum2<-as.data.frame(rbind(t11,t22))
        write.csv(sum2, paste0("Output_all/Difference/Conserved/Summary.stats.higher.in.",g2,".csv"))
        
        #syn vs. nonsyn proportion
        sink(paste0("Output_all/Difference/Conserved/syn-nonsyn-ratio.",compare,".txt"))
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
        write.csv(mf.top,paste0("Output_all/Difference/Conserved/Top5percent.MostDifferent.MF.", compare,".csv"))
        
        
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
        write.csv(PercentHighDiff,paste0("Output_all/Difference/Conserved/Top5percent.MostDifferent.SummarybyGene.",compare,".csv"))
        
        ggplot(PercentHighDiff,aes(x=Gene,y=Proportion))+geom_bar(stat='identity', fill='#66CCEECC',width=0.8)+
                labs(x="Genes",y=paste0("% highly different sites between ",g1," & ",g2))+
                theme_classic()
        ggsave(filename=paste0("Output_all/Difference/Conserved/Top5percent.MostDifferent.MF.by.gene.",compare,".pdf"), width = 5.7, height = 4.5)
        
        
        
        #absolute difference top 1%
        s1<-as.integer(0.01*nrow(mf1))
        mf.top2<-mf.top[c(1:s1),]
        
        write.csv(mf.top2,paste0("Output_all/Difference/Conserved/Top1percent.MostDifferent.MF.",compare,".csv"))
        
        PercentHighDiff2<-data.frame("Gene"= paste0(genes$Gene[2:12]), "Proportion"=rep('',times=11), stringsAsFactors = FALSE)
        for (i in 2:12){
                n<-Dat[Dat$gene==genes$Gene[i],]
                df<-mf.top2[mf.top2$gene==genes$Gene[i],]
                
                PercentHighDiff2$Proportion[i-1] <- nrow(df)/(nrow(n))*100
                PercentHighDiff2$Counts[i-1]<-nrow(df)
                PercentHighDiff2$No.Sites[i-1]<-nrow(n)
        }
        PercentHighDiff2$Proportion<-as.numeric( PercentHighDiff2$Proportion)
        write.csv(PercentHighDiff2,paste0("Output_all/Difference/Conserved/Top1percent.MostDifferent.SummarybyGene.",compare,".csv"))
        
        ggplot(PercentHighDiff2,aes(x=Gene,y=Proportion))+geom_bar(stat='identity', fill='#EE667799',width=0.8)+
                labs(x="Genes",y=paste0("% of sites highly different in MF betw/ ",g1," & ",g2))+
                theme_classic()
        ggsave(filename=paste0("Output_all/Difference/Conserved/Top1per.MostDifferent.MF.by.gene.",compare,".pdf"),width = 5.7, height = 4.5)

        
        
        #Which Amino acids?
        for (j in 1:2){
                if (j==1) G<-g1; DF<-mf.H; G2<-g2
                if (j==2) G<-g2;  DF<-mf.L; G2<-g1
                aa1<-table(DF[,paste0("WTAA.",G)],DF[,paste0("Type.",G)])
                
                #Across the entire genome
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
                AA1$genome.percent.nonsyn<-AAtotal[,1]/sum(AAtotal[,1])*100
                AA1$percent.nonsyn<-AA1$nonsyn/sum(AA1$nonsyn)*100
                AA1$genome.percent.syn<-AAtotal[,3]/sum(AAtotal[,3])*100
                AA1$percent.syn<-AA1$syn/sum(AA1$syn)*100
                AA1.1<-AA1[,c(1,4:7)]
                colnames(AA1.1)<-c("AA","% nonsyn genome-wide", "% nonsyn", "% syn genome-wide", "% syn")
                AA1.1<-melt(AA1.1)
                colnames(AA1.1)[2:3]<-c("Group","Proportion")
                
                ggplot(AA1.1, aes(x=AA, y=Proportion, fill=Group))+geom_bar(stat='identity', position='dodge')+
                        labs(x="Amino acids", y="% of represented AA")+
                        ggtitle(paste0("% of AA in highly different MF sites: higher in ",G,"/lower in ",G2, " (", compare,")"))+
                        scale_fill_manual("", values = cols3)+theme_classic()+
                        theme(plot.title = element_text(hjust = 0.5, size=12))
                ggsave(filename=paste0("Output_all/Difference/Conserved/Which.AA.in.HigherMFof",G,".in.", compare,".pdf"),width =9, height = 4.5)
                #nonsyn only
                ggplot(AA1.1[AA1.1$Group=="% nonsyn genome-wide"|AA1.1$Group=="% nonsyn",], aes(x=AA, y=Proportion, fill=Group))+geom_bar(stat='identity', position='dodge')+
                        labs(x="Amino acids", y="% of represented AA")+
                        ggtitle(paste0("% of AA in highly different MF sites: higher in ",G,"/lower in ",G2, " (", compare,")"))+
                        scale_fill_manual("", values = cols3[1:2])+theme_classic()+
                        theme(plot.title = element_text(hjust = 0.5, size=12))
                ggsave(filename=paste0("Output_all/Difference/Conserved/Which.AA.in.HigherMFof",G,".in.", compare,".Nonsyn.pdf"),width = 9, height = 4.5)
#
                #show overrepresentation only
                AA1$diff.nonsyn<-(AA1$percent.nonsyn-AA1$genome.percent.nonsyn)/AA1$genome.percent.nonsyn*100
                ggplot(AA1,aes(x=AA, y=diff.nonsyn))+geom_bar(stat='identity',fill='#66CCEE')+
                        theme_classic()+labs(x="Amino acids", y="% difference")+
                        ggtitle(paste0("Over/Under-represented AAs of highly different MF sites: higher in ", G, " (",compare,")"))+
                        theme(plot.title = element_text(size=12))
                ggsave(filename=paste0("Output_all/Difference/Conserved/Overrepresented.AAs: higher in ",G,".in.",compare,".pdf"), width = 7, height = 4.5)
        }      
        
        #look at the 5% most different sites (abs, top5)
        
        DF<-mf.top
        aalist<-DF[,c("WTAA.1A","WTAA.1B","WTAA.3A")]
        
        aa1<-table(DF[,paste0("WTAA.",g1)],DF[,paste0("Type.",g1)])
        aa2<-table(DF[,paste0("WTAA.",g2)],DF[,paste0("Type.",g2)])
        
        #Across the entire genome
        AAgenome1<-table(mf1[mf1$gene!="5' UTR",paste0("WTAA.",g1)],mf1[mf1$gene!="5' UTR",paste0("Type.",g1)] )
        AAgenome2<-table(mf1[mf1$gene!="5' UTR",paste0("WTAA.",g2)],mf1[mf1$gene!="5' UTR",paste0("Type.",g2)] )
        AAgenome.comb<-cbind(AAgenome1,AAgenome2)
        colnames(AAgenome.comb)<-c(paste0("nonsyn.",g1),paste0("stop.",g1),paste0("syn.",g1),paste0("nonsyn.",g2),paste0("stop.",g2),paste0("syn.",g2))        
        
        Genome.comb.nonsyn<-(AAgenome.comb[,1]+AAgenome.comb[,4])/(sum(AAgenome.comb[,1])+sum(AAgenome.comb[,4]))
        Genome.comb.syn   <-(AAgenome.comb[,3]+AAgenome.comb[,6])/(sum(AAgenome.comb[,3])+sum(AAgenome.comb[,6]))
        genome.AAcom<-as.data.frame(cbind(Genome.comb.nonsyn,Genome.comb.syn))
        AA<-data.frame(AA=rownames(genome.AAcom))
        
        for (a in 1:2) {
                AA1<-AA
                G<-get(paste0("g",a))
                aa<-get(paste0("aa",a))
                AAgenome<-get(paste0("AAgenome",a))
                AA.df<-data.frame(AA=rownames(aa))
                AA.df$nonsyn<-aa[,"nonsyn"]
                AA.df$syn<-aa[,"syn"]
                AA1<-merge(AA1, AA.df, by="AA",all.x=T)
                AA1[is.na(AA1)]<-0
                
                AA1$genome.percent.nonsyn<-AAgenome[,1]/sum(AAgenome[,1])*100
                AA1$genome.percent.syn<-AAgenome[,3]/sum(AAgenome[,3])*100
                AA1$percent.nonsyn<-AA1$nonsyn/sum(AA1$nonsyn)*100
                AA1$percent.syn<-AA1$syn/sum(AA1$syn)*100
        
        
                AA1.1<-AA1[,c(1,4:7)]
                colnames(AA1.1)<-c("AA","% nonsyn genome-wide", "% nonsyn", "% syn genome-wide", "% syn")
                AA1.1<-melt(AA1.1)
                colnames(AA1.1)[2:3]<-c("Group","Proportion")
                
                ggplot(AA1.1, aes(x=AA, y=Proportion, fill=Group))+geom_bar(stat='identity', position='dodge')+
                        labs(x="Amino acids", y="% of represented AA")+
                        ggtitle(paste0("% of AA represented in highly different MF sites (5%) of ",G," in ", compare))+
                        scale_fill_manual("", values = cols3)+theme_classic()+
                        theme(plot.title = element_text(hjust = 0.5, size=12))
                ggsave(filename=paste0("Output_all/Difference/Conserved/Top5_AAsinDiff.",G,".", compare,".pdf"),width =9, height = 4.5)
                #nonsyn only
                ggplot(AA1.1[AA1.1$Group=="% nonsyn genome-wide"|AA1.1$Group=="% nonsyn",], aes(x=AA, y=Proportion, fill=Group))+geom_bar(stat='identity', position='dodge')+
                        labs(x="Amino acids", y="% of represented AA")+
                        ggtitle(paste0("% of AA represented in highly different MF sites (5%) of ",G," in ", compare))+
                        scale_fill_manual("", values = cols3[1:2])+theme_classic()+
                        theme(plot.title = element_text(hjust = 0.5, size=12))
                ggsave(filename=paste0("Output_all/Difference/Conserved/Top5_AAsinDiff.nonsyn.",G,".", compare,".Nonsyn.pdf"),width = 7, height = 4.5)
                #
                #show overrepresentation only
                AA1$diff.nonsyn<-(AA1$percent.nonsyn-AA1$genome.percent.nonsyn)/AA1$genome.percent.nonsyn*100
                ggplot(AA1,aes(x=AA, y=diff.nonsyn))+geom_bar(stat='identity',fill='#66CCEE')+
                        theme_classic()+labs(x="Amino acids", y="% difference")+
                        ggtitle(paste0("Over-/Under-represented AAs of highly different MF sites (5%) of ", G, " (",compare,")"))+
                        theme(plot.title = element_text(size=12))
                ggsave(filename=paste0("Output_all/Difference/Conserved/Top5_Overrepresented.AAs.",G,".",compare,".pdf"), width = 7, height = 4.5)
        }      
        
        
        ### only look at the same AA sites
        sameAA<-which(aalist[,paste0("WTAA.",g1)] == aalist[,paste0("WTAA.",g2)])
        DF2<-DF[sameAA,]
        aa1<-table(DF2[,paste0("WTAA.",g1)],DF2[,paste0("Type.",g1)])
        aa2<-table(DF2[,paste0("WTAA.",g2)],DF2[,paste0("Type.",g2)])

        AA.df<-data.frame(AA=rownames(aa1))
        AA.df$nonsyn<-aa1[,"nonsyn"]
        AA.df$syn<-aa1[,"syn"]
        AA1<-merge(AA, AA.df, by="AA",all.x=T)
        AA1[is.na(AA1)]<-0
                
        AA1$genome.percent.nonsyn<-genome.AAcom[,1]*100
        AA1$percent.nonsyn<-AA1$nonsyn/sum(AA1$nonsyn)*100
        AA1$genome.percent.syn<-genome.AAcom[,2]*100
        AA1$percent.syn<-AA1$syn/sum(AA1$syn)*100
                
        AA1.1<-AA1[,c(1,4:7)]
        colnames(AA1.1)<-c("AA","% nonsyn genome-wide", "% nonsyn", "% syn genome-wide", "% syn")
        AA1.1<-melt(AA1.1)
        colnames(AA1.1)[2:3]<-c("Group","Proportion")
        
        ggplot(AA1.1, aes(x=AA, y=Proportion, fill=Group))+geom_bar(stat='identity', position='dodge')+
                labs(x="Amino acids", y="% of represented AA")+
                ggtitle(paste0("% of AA in highly different MF sites (5%) between ",g1," & ", g2))+
                scale_fill_manual("", values = cols3)+theme_classic()+
                theme(plot.title = element_text(hjust = 0.5, size=12))
        ggsave(filename=paste0("Output_all/Difference/Conserved/SameAAonly_AAsinTop5Diff.", compare,".pdf"),width =9, height = 4.5)
        #nonsyn only
        ggplot(AA1.1[AA1.1$Group=="% nonsyn genome-wide"|AA1.1$Group=="% nonsyn",], aes(x=AA, y=Proportion, fill=Group))+geom_bar(stat='identity', position='dodge')+
                labs(x="Amino acids", y="% of represented AA")+
                ggtitle(paste0("% of AA in highly different MF sites (5%) between ",g1," & ", g2))+
                scale_fill_manual("", values = cols3[1:2])+theme_classic()+
                theme(plot.title = element_text(hjust = 0.5, size=12))
        ggsave(filename=paste0("Output_all/Difference/Conserved/SameAAonly_AAsinTop5Diff.nonsyn.", compare,".pdf"),width = 9, height = 4.5)
        #
        #show overrepresentation only
        AA1$diff.nonsyn<-(AA1$percent.nonsyn-AA1$genome.percent.nonsyn)/AA1$genome.percent.nonsyn*100
        write.csv(AA1, paste0("Output_all/Difference/Conserved/AAinHighlyDiffrentSites_", compare,".csv"))
        
        ggplot(AA1,aes(x=AA, y=diff.nonsyn))+geom_bar(stat='identity',fill='#66CCEE')+
                theme_classic()+labs(x="Amino acids", y="Over/under-represented AA (%)")+
                ggtitle(paste0("Over/Under-represented AAs of highly different MF sites (5%) between ",g1," & ", g2))+
                theme(plot.title = element_text(size=12))+
                scale_fill_manual("", values = cols3)+theme_classic()
        ggsave(filename=paste0("Output_all/Difference/Conserved/SameAAonly_Overrepresented.AAs.top5.",compare,".pdf"), width = 7, height = 4.5)
        
}


aafile<-list.files("Output_all/Difference/Conserved/", pattern = "AAinHighlyDiffrentSites")
aadiff<-list()
for (i in 1:length(aafile)){
        dt<-read.csv(paste0("Output_all/Difference/Conserved/", aafile[i]), stringsAsFactors = F, row.names = 1)
        aadiff[[i]]<-dt[,c("AA","diff.nonsyn")]
        names(aadiff)[i]<-substr(aafile[i], start=25, stop = 33)
} 

for (i in 1:length(aadiff)) {
        colnames(aadiff[[i]])<-c("AA",paste0(names(aadiff[i])))
}

AAdiff<-aadiff%>% purrr::reduce(full_join, by='AA')
AAdiffm<-melt(AAdiff)
ggplot(AAdiffm,aes(x=AA, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_classic()+labs(x="Amino acids", y="Over/under-represented AA (%)")+
        #ggtitle(paste0("Over/Under-represented AAs of highly different MF sites (5%) between ",g1," & ", g2))+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols3[c(2,1,3)])+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/Difference/Conserved/Overrepresented.AAs.All.top5.pdf"), width = 7, height = 4.5)


genefile<-list.files("Output_all/Difference/Conserved/", pattern="Top5percent.MostDifferent.SummarybyGene")
geneprop<-list()
for (i in 1:length(genefile)){
        dt<-read.csv(paste0("Output_all/Difference/Conserved/", genefile[i]), stringsAsFactors = F, row.names = 1)
        geneprop[[i]]<-dt[,c("Gene","Proportion")]
        names(geneprop)[i]<-substr(genefile[i], start=41, stop = 49)
} 
for (i in 1:length(geneprop)) {
        colnames(geneprop[[i]])<-c("Gene",paste0(names(geneprop[i])))
}

GeneProp<-geneprop%>% purrr::reduce(full_join, by='Gene')
GenePropm<-melt(GeneProp)
ggplot(GenePropm,aes(x=Gene, y=value, fill=variable))+geom_bar(stat='identity', position = "dodge", width = 0.7)+
        theme_classic()+labs(x="", y="% highly different sites between genotypes")+
        #ggtitle(paste0("Over/Under-represented AAs of highly different MF sites (5%) between ",g1," & ", g2))+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols3[c(2,1,3)])+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/Difference/Conserved/All.Top5Proportion.per.gene2.pdf"), width = 7, height = 4.5)




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
write.csv(ratio,"Output_all/Difference/Conserved/NStoSynRatio.csv")

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

write.csv(SNratio,"Output_all/Difference/Conserved/NStoSynRatio.by.gene.csv")

colnames(SNratio)[2:4]<-c("1A","1B","3A")
SNratio$Gene<-as.character(SNratio$Gene)
SNratio[5,1]<-"NS1"
SNratio2<-melt(SNratio)
colnames(SNratio2)[2:3]<-c("Genotype","ratio")

ggplot(SNratio2,aes(x=Gene,y=ratio,group=Genotype, color=Genotype))+geom_point(position=position_dodge(width=0.4), size=4)+
        scale_color_manual(values=cols2)+
        theme_light()+
        theme(axis.title.x=element_blank())+ylab("Ratio of nonsyn to syn mutation frequency")
ggsave(filename="Output_all/Difference/Conserved/NSratio.by.gene.by.genotype.pdf",width = 9, height = 8)
ggsave(filename="Output_all/Difference/Conserved/NSratio.by.gene.by.genotype2.pdf",width = 5.5, height = 4)







