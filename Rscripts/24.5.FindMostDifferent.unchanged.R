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

cols3<-c("#66CCEE66","#66CCEE","#EE667766","#EE6677")


#############
# find the sites that are most different in mutation frequencies between the genotypes

mf<-read.csv("Output_all/Unfiltered/Ts.Same_Mean_3genotypes.csv", stringsAsFactors = F , row.names = 1)
mf<-mf[mf$merged.pos>264&mf$merged.pos<=8608,]

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
        mf1<-mf1[mf1[n1]==mf1[n2],]
        mf1<-mf1[!is.na(mf1[,paste0("mean.",g1)]),]
        mf1<-mf1[!is.na(mf1[,paste0("mean.",g2)]),]
        print(nrow(mf1)) #1A-1B:6426    1A-3A:5533  1B-3A:5657
        s<-as.integer(0.05*nrow(mf1)) #top 5%
        
        mf1$diff<-mf1[,paste0("mean.",g1)]-mf1[,paste0("mean.",g2)]
        pdf(paste0("Output_all/Difference/Conserved/Hist.diff.in.mutfreq.",compare,".pdf"),width=6,height=5)
        hist(mf1$diff, breaks=40, xlab="Difference in mutation frequency", main=compare)
        dev.off()
        
        ## absolute top 5% different sites
        mf1$diff.abs<-abs(mf1$diff)
        mf.top<-mf1[order(mf1$diff.abs, decreasing = T,na.last = T),]
        #exclude "5'UTR"
        mf.top<-mf.top[mf.top$gene!="5' UTR",]
        mf.top<-mf.top[c(1:s),]
        write.csv(mf.top,paste0("Output_all/Difference/Conserved/Top5percent.MostDifferent.MF.", compare,".csv"))
        
        
        #Over/under-representation of highly different sites by gene
        PercentHighDiff<-data.frame("Gene"= paste0(genes$Gene[2:12]), stringsAsFactors = FALSE)
        Dat<-mf1[,c("merged.pos", "gene","diff.abs")]
        Dat<-Dat[!is.na(Dat$diff.abs),]
        
        for (i in 2:12){
                n<-Dat[Dat$gene==genes$Gene[i],]
                df<-mf.top[mf.top$gene==genes$Gene[i],]
                PercentHighDiff$Counts[i-1]<-nrow(df)
                PercentHighDiff$Total[i-1]<-nrow(n)
                PercentHighDiff$Expected[i-1]<-nrow(n)*s/nrow(Dat)
                PercentHighDiff$Difference[i-1] <-(PercentHighDiff$Counts[i-1]-PercentHighDiff$Expected[i-1])/PercentHighDiff$Expected[i-1]*100
        }
        
        write.csv(PercentHighDiff,paste0("Output_all/Difference/Conserved/Top5percent.MostDifferent.SummarybyGene.",compare,".csv"))
        
        ggplot(PercentHighDiff,aes(x=Gene,y=Difference))+geom_bar(stat='identity', fill='#66CCEECC',width=0.8)+
                labs(x="Genes",y=paste0("Over/Under-represented % between ",g1," & ",g2))+
                theme_classic()
        ggsave(filename=paste0("Output_all/Difference/Conserved/Top5percent.MostDifferent.MF.by.gene.",compare,".pdf"), width = 5.7, height = 4.5)
        
        
        
        #Which Amino acids?

        #look at the 5% most different sites (abs, top5)
        DF<-mf.top
        aalist<-DF[,c("WTAA.1A","WTAA.1B","WTAA.3A")]
        # only look at the same AA sites
        sameAA<-which(aalist[,paste0("WTAA.",g1)] == aalist[,paste0("WTAA.",g2)])
        DF<-mf.top[sameAA,]
        # checking: should be aa1 == aa2
        aa1<-table(DF[,paste0("WTAA.",g1)],DF[,paste0("Type.",g1)])
        aa2<-table(DF[,paste0("WTAA.",g2)],DF[,paste0("Type.",g2)])
        
        #genome-wide proportion of AA 
        AAgenome1<-table(mf1[mf1$gene!="5' UTR",paste0("WTAA.",g1)],mf1[mf1$gene!="5' UTR",paste0("Type.",g1)] )
        AAgenome2<-table(mf1[mf1$gene!="5' UTR",paste0("WTAA.",g2)],mf1[mf1$gene!="5' UTR",paste0("Type.",g2)] )
        AAgenome.comb<-cbind(AAgenome1,AAgenome2)
        colnames(AAgenome.comb)<-c(paste0("nonsyn.",g1),paste0("stop.",g1),paste0("syn.",g1),paste0("nonsyn.",g2),paste0("stop.",g2),paste0("syn.",g2))        
        
        Genome.comb.nonsyn<-(AAgenome.comb[,1]+AAgenome.comb[,4])/(sum(AAgenome.comb[,1])+sum(AAgenome.comb[,4]))
        Genome.comb.syn   <-(AAgenome.comb[,3]+AAgenome.comb[,6])/(sum(AAgenome.comb[,3])+sum(AAgenome.comb[,6]))
        genome.AAcom<-as.data.frame(cbind(Genome.comb.nonsyn,Genome.comb.syn))
        
        AA<-data.frame(AA=rownames(genome.AAcom))
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
        colnames(AA1.1)<-c("AA","Nonsyn genome-wide", "Nonsyn highly-different", "Syn genome-wide", "Syn highly-different")
        AA1.1<-melt(AA1.1)
        colnames(AA1.1)[2:3]<-c("Group","Proportion")
        
        ggplot(AA1.1, aes(x=AA, y=Proportion, fill=Group))+geom_bar(stat='identity', position='dodge')+
                labs(x="Amino acids", y="% of represented AA")+
                ggtitle(paste0("% of AA in highly different MF sites (5%) between ",g1," & ", g2))+
                scale_fill_manual("", values = cols3)+theme_classic()+
                theme(plot.title = element_text(hjust = 0.5, size=12))
        ggsave(filename=paste0("Output_all/Difference/Conserved/SameAAonly_AAsinTop5Diff.", compare,".pdf"),width =9, height = 4.5)
        #nonsyn only
        ggplot(AA1.1[AA1.1$Group=="Nonsyn genome-wide"|AA1.1$Group=="Nonsyn highly-different",], aes(x=AA, y=Proportion, fill=Group))+geom_bar(stat='identity', position='dodge')+
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
        
        n<-nrow(mf1)
        
        #type
        type1<-as.data.frame(table(mf.top[,paste0("Type.",g1)]))
        type2<-as.data.frame(table(mf.top[,paste0("Type.",g2)]))
        types<-merge(type1, type2, by="Var1")
        colnames(types)<-c("Type",g1,g2)
        gtype1<-as.data.frame(table(mf1[,paste0("Type.",g1)]))
        gtype2<-as.data.frame(table(mf1[,paste0("Type.",g2)]))
        gtypes<-merge(gtype1, gtype2, by="Var1")
        colnames(gtypes)<-c("Type", paste0("genome.",g1),paste0("genome.",g2) )
        types<-merge(types,gtypes, by="Type")
        types$difference1<-(types[,g1]-(types[,paste0("genome.",g1)]*s/n))/(types[,paste0("genome.",g1)]*s/n)*100
        types$difference2<-(types[,g2]-(types[,paste0("genome.",g2)]*s/n))/(types[,paste0("genome.",g2)]*s/n)*100
        types$difference<-rowMeans(types[,c("difference1", "difference2")])
        write.csv(types,paste0("Output_all/Difference/Conserved/Types_", compare,".csv"))
        
        #by codon positions
        Codons<-as.data.frame(table(mf.top$codon))
        gCodons<-as.data.frame(table(mf1$codon))
        Codons<-merge(Codons,gCodons, by="Var1")
        colnames(Codons)<-c("codon","high.diff", "genome")
        
        Codons$difference<-(Codons$high.diff-(Codons$genome*s/n))/(Codons$genome*s/n)*100
        write.csv(Codons, paste0("Output_all/Difference/Conserved/Codons_", compare,".csv"))
        
        #Nucleotide
        Nt<-as.data.frame(table(mf.top[,paste0("ref.",g1)]))
        gNt<-as.data.frame(table(mf1[,paste0("ref.",g1)]))
        Nt<-merge(Nt,gNt, by="Var1")
        colnames(Nt)<-c("nucleotide","high.diff", "genome")
        Nt$difference<-(Nt$high.diff-(Nt$genome*s/n))/(Nt$genome*s/n)*100
        
        write.csv(Nt,paste0("Output_all/Difference/Conserved/NT_", compare,".csv"))
        
        #CpG Sites
        cpg1<-as.data.frame(table(mf.top[,paste0("makesCpG.",g1)]))
        cpg2<-as.data.frame(table(mf.top[,paste0("makesCpG.",g2)]))
        cpg<-merge(cpg1,cpg2,by="Var1")
        colnames(cpg)<-c("CpG",g1,g2)
        Gcpg1<-as.data.frame(table(mf1[,paste0("makesCpG.",g1)]))
        Gcpg2<-as.data.frame(table(mf1[,paste0("makesCpG.",g2)]))
        Gcpg<-merge(Gcpg1,Gcpg2,by="Var1")
        colnames(Gcpg)<-c("CpG",g1,g2)
        cpg<-merge(cpg,Gcpg, by="CpG")
        cpg$diff1<-(cpg[,paste0(g1,".x")]-(cpg[,paste0(g1,".y")]*s/n))/(cpg[,paste0(g1,".y")]*s/n)*100
        cpg$diff2<-(cpg[,paste0(g2,".x")]-(cpg[,paste0(g2,".y")]*s/n))/(cpg[,paste0(g2,".y")]*s/n)*100
        cpg$difference<-rowMeans(cpg[,c("diff1","diff2")])     
        write.csv(cpg,paste0("Output_all/Difference/Conserved/CpG_", compare,".csv"))
        
        
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
        theme_classic()+labs(x="Amino acid", y="Over/under-represented AA (%)")+
        #ggtitle(paste0("Over/Under-represented AAs of highly different MF sites (5%) between ",g1," & ", g2))+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols4[c(2,1,3)])+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/Difference/Conserved/Overrepresented.AAs.All.top5.pdf"), width = 7, height = 4.5)


genefile<-list.files("Output_all/Difference/Conserved/", pattern="Top5percent.MostDifferent.SummarybyGene")
geneprop<-list()
for (i in 1:length(genefile)){
        dt<-read.csv(paste0("Output_all/Difference/Conserved/", genefile[i]), stringsAsFactors = F, row.names = 1)
        geneprop[[i]]<-dt[,c("Gene","Difference")]
        names(geneprop)[i]<-substr(genefile[i], start=41, stop = 49)
} 
for (i in 1:length(geneprop)) {
        colnames(geneprop[[i]])<-c("Gene",paste0(names(geneprop[i])))
}

GeneProp<-geneprop%>% purrr::reduce(full_join, by='Gene')
GenePropm<-melt(GeneProp)
ggplot(GenePropm,aes(x=Gene, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_classic()+labs(x="", y="Over/under-represented sites per gene (%)")+
        #ggtitle(paste0("Over/Under-represented AAs of highly different MF sites (5%) between ",g1," & ", g2))+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols4[c(2,1,3)])+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/Difference/Conserved/All.Top5Proportion.per.gene.pdf"), width = 7, height = 4.5)

## codon

codonfile<-list.files("Output_all/Difference/Conserved/", pattern = glob2rx("*Codons*.csv*"))
codondiff<-list()
for (i in 1:length(codonfile)){
        dt<-read.csv(paste0("Output_all/Difference/Conserved/", codonfile[i]), stringsAsFactors = F, row.names = 1)
        codondiff[[i]]<-dt[,c("codon","difference")]
        names(codondiff)[i]<-substr(codonfile[i], start=8, stop = 16)
} 

for (i in 1:length(codondiff)) {
        colnames(codondiff[[i]])<-c("codon",paste0(names(codondiff[i])))
}

Codondiff<-codondiff%>% purrr::reduce(full_join, by='codon')
Codondiff$codon<-as.factor(Codondiff$codon)
Codondiffm<-melt(Codondiff)
ggplot(Codondiffm,aes(x=codon, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_classic()+labs(x="Codon position", y="Over/under-represented codon (%)")+
        #ggtitle(paste0("Over/Under-represented codons of highly different MF sites (5%) between ",g1," & ", g2))+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols4[c(2,1,3)])+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/Difference/Conserved/Overrepresented.Codons.all.top5.pdf"), width = 5, height = 3.2)



## NT

ntfile<-list.files("Output_all/Difference/Conserved/", pattern = glob2rx("*NT*.csv*"))
ntdiff<-list()
for (i in 1:length(ntfile)){
        dt<-read.csv(paste0("Output_all/Difference/Conserved/", ntfile[i]), stringsAsFactors = F, row.names = 1)
        ntdiff[[i]]<-dt[,c("nucleotide","difference")]
        names(ntdiff)[i]<-substr(ntfile[i], start=4, stop = 12)
} 

for (i in 1:length(ntdiff)) {
          colnames(ntdiff[[i]])<-c("nt",paste0(names(ntdiff[i])))
}

ntdiff<-ntdiff%>% purrr::reduce(full_join, by='nt')
ntdiffm<-melt(ntdiff)
ggplot(ntdiffm,aes(x=nt, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_classic()+labs(x="", y="Over/under-represented nucleotide (%)")+
        #ggtitle(paste0("Over/Under-represented codons of highly different MF sites (5%) between ",g1," & ", g2))+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols4[c(2,1,3)])+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/Difference/Conserved/Overrepresented.NT.all.top5.pdf"), width = 7, height = 4.5)



## Type

tpfile<-list.files("Output_all/Difference/Conserved/", pattern =glob2rx("*Types*.csv*"))
tpdiff<-list()
for (i in 1:length(tpfile)){
        dt<-read.csv(paste0("Output_all/Difference/Conserved/", tpfile[i]), stringsAsFactors = F, row.names = 1)
        tpdiff[[i]]<-dt[,c("Type","difference")]
        names(tpdiff)[i]<-substr(tpfile[i], start=7, stop = 15)
} 

for (i in 1:length(tpdiff)) {
        colnames(tpdiff[[i]])<-c("Type",paste0(names(tpdiff[i])))
}

tpdiff<-tpdiff%>% purrr::reduce(full_join, by='Type')
tpdiffm<-melt(tpdiff)
ggplot(tpdiffm,aes(x=Type, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_classic()+labs(x="Mutation type", y="Over/under-represented mutation type (%)")+
        #ggtitle(paste0("Over/Under-represented codons of highly different MF sites (5%) between ",g1," & ", g2))+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols4[c(2,1,3)])+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/Difference/Conserved/Overrepresented.Type.all.top5.pdf"), width = 5, height = 3.3)

##CpG
cpgfile<-list.files("Output_all/Difference/Conserved/", pattern = glob2rx("*CpG*.csv*"))
cpgdiff<-list()
for (i in 1:length(cpgfile)){
        dt<-read.csv(paste0("Output_all/Difference/Conserved/", cpgfile[i]), stringsAsFactors = F, row.names = 1)
        cpgdiff[[i]]<-dt[,c("CpG","difference")]
        names(cpgdiff)[i]<-substr(cpgfile[i], start=5, stop = 12)
} 

for (i in 1:length(cpgdiff)) {
        colnames(cpgdiff[[i]])<-c("CpG",paste0(names(cpgdiff[i])))
}

cpgdiff<-cpgdiff%>% purrr::reduce(full_join, by='CpG')
cpgdiff$CpG<-c("Non-CpG","CpG")
cpgdiffm<-melt(cpgdiff)
ggplot(cpgdiffm,aes(x=CpG, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_classic()+labs(x="", y="Over/under-represented CpG sites (%)")+
        #ggtitle(paste0("Over/Under-represented CpG sites of highly different MF sites (5%) between ",g1," & ", g2))+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols4[c(2,1,3)])+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/Difference/Conserved/Overrepresented.CpG.all.top5.pdf"), width = 4.5, height = 3)




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

ggplot(SNratio2,aes(x=Gene,y=ratio,group=Genotype, color=Genotype))+geom_point(position=position_dodge(width=0.4), size=3.5)+
        scale_color_manual(values=cols2)+
        theme_bw()+
        theme(axis.title.x=element_blank())+ylab("Ratio of nonsyn to syn mutation frequency")+
        theme(panel.grid.major.x=element_blank())+
        theme(panel.grid.minor.y=element_blank())+
        theme(panel.grid.major.y=element_line(linetype="dashed"))+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray60", size=.2)
        
#ggsave(filename="Output_all/Difference/Conserved/NSratio.by.gene.by.genotype.pdf",width = 9, height = 8)
ggsave(filename="Output_all/Difference/Conserved/NSratio.by.gene.by.genotype.pdf",width = 5.5, height = 4)







