library(stringr)
library(tidyverse)
library(zoo)
library(reshape2)
library(ggplot2)
library(DescTools)
library(reshape2)
source("Rscripts/baseRscript.R")

#dir.create("Output_all/Fst/")

cols2.60<-paste0(cols2,"99")
cols3<-c("#009988CC" ,"#66CCEECC", "#EE6677CC", "#4477AACC")

genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
geno<-c("1A","1B","3A")

FstList<-list()
fstfiles<-list.files("Output_all/Fst/", pattern=glob2rx("*^Fst*.csv*"))
for (i in 1:3){
        FstList[[i]]<-read.csv(paste0("Output_all/Fst/",fstfiles[i]),stringsAsFactors = F, row.names = 1)
        names(FstList)[i]<-substr(fstfiles[i], start=1, stop = 9)
}

merged.meta<-read.csv("Output_all/merged.metadata.csv", stringsAsFactors = F, row.names = 1)
merged.meta$gene[merged.meta$gene=="NS1(P7)"]<-"NS1"


for (g in 1:3){
        dt<-FstList[[g]]
        print(paste(geno[g],nrow(dt)))
        fname<-substr(names(FstList)[g], start=5, stop = 9)
        g1<-substr(names(FstList)[g], start=5, stop = 6)
        g2<-substr(names(FstList)[g], start=8, stop = 9)
        #top5%
        s<-as.integer(nrow(dt)*0.05)
        #attached the metadata
        dt<-merge(dt, merged.meta, by="merged.pos")
        #dt<-dt[dt$Type.1A!="stop",]
        #dt<-dt[dt$Type.1B!="stop",]
        #dt<-dt[dt$Type.3A!="stop",]
        
        dt.top<-dt[order(dt$Fst,decreasing = T,na.last = T),]
        dt.top<-dt.top[c(1:s),]
        write.csv(dt.top,paste0("Output_all/Fst/Top5percent.", names(FstList)[g],".csv"))
        
        #nonsyn sites in both genotypes only
        dt1<-dt[dt[,paste0("Type.",g1)]=="nonsyn" & dt[,paste0("Type.",g2)]=="nonsyn",]
        dt1.top<-dt1[order(dt1$Fst,decreasing = T,na.last = T),]
        dt1.top<-dt1.top[c(1:s),]
        write.csv(dt1.top,paste0("Output_all/Fst/Top5percent.NS..", fname,".csv"))
        #dt.topNS<-dt.top[dt.top[,paste0("Type.",g1)]=="nonsyn" & dt.top[,paste0("Type.",g2)]=="nonsyn" ,]
        
        #check proportion of sites belonging to each gene
        #1. all
        highFst<-data.frame("Gene"= paste0(genes$Gene[1:12]), stringsAsFactors = FALSE)
        Dat<-dt[!is.na(dt$Fst),]
        for (i in 1:12){
                n<-Dat[Dat$gene==genes$Gene[i],]
                df<-dt.top[dt.top$gene==genes$Gene[i],]
                highFst$Counts[i]<-nrow(df)
                highFst$Total[i]<-nrow(n)
                highFst$Expected[i]<-nrow(n)*s/nrow(Dat)
                highFst$Difference[i] <-(highFst$Counts[i]-highFst$Expected[i])/highFst$Expected[i]*100
        }
        
        write.csv(highFst,paste0("Output_all/Fst/Top5.Fst.SummarybyGene.",fname,".csv"))
        
        ggplot(highFst,aes(x=Gene,y=Difference))+geom_bar(stat='identity', fill='#66CCEECC',width=0.8)+
                labs(x="Genes",y=paste0("Over/Under-represented sites with high-Fst (",g1, "-",g2))+
                theme_classic()
        ggsave(filename=paste0("Output_all/Fst/Top5.Fst.by.gene.",fname,".pdf"), width = 5.7, height = 4.5)
        
        Gtest<-data.frame(test=c("gene","AA-ns","AA-syn","Type", "codon","nt","CpG"))
        k=1
        r1<-GTest(highFst[,2:3])
        Gtest$G[k]<-r1[[1]]
        Gtest$Pvalue[k]<-r1[[3]]
        
        #2. NS only
        highFst1<-data.frame("Gene"= paste0(genes$Gene[1:12]), stringsAsFactors = FALSE)
        Dat1<-dt1[!is.na(dt1$Fst),]
        for (i in 1:12){
                n<-Dat1[Dat1$gene==genes$Gene[i],]
                df<-dt1.top[dt1.top$gene==genes$Gene[i],]
                highFst1$Counts[i]<-nrow(df)
                highFst1$Total[i]<-nrow(n)
                highFst1$Expected[i]<-nrow(n)*s/nrow(Dat1)
                highFst1$Difference[i] <-(highFst1$Counts[i]-highFst1$Expected[i])/highFst1$Expected[i]*100
        }
        
        write.csv(highFst1,paste0("Output_all/Fst/NS/Top5.NS.Fst.SummarybyGene.",fname,".csv"))
        
        ggplot(highFst1,aes(x=Gene,y=Difference))+geom_bar(stat='identity', fill='#66CCEECC',width=0.8)+
                labs(x="Genes",y=paste0("Over/Under-represented sites with high-Fst (",g1, "-",g2))+
                theme_classic()
        ggsave(filename=paste0("Output_all/Fst/NS/Top5.ns.Fst.by.gene.",fname,".pdf"), width = 5.7, height = 4.5)
        
        Gtest1<-data.frame(test=c("gene","AA-ns","codon","nt","CpG"))
        j=1
        r1<-GTest(highFst1[,2:3])
        Gtest1$G[j]<-r1[[1]]
        Gtest1$Pvalue[j]<-r1[[3]]
        
        
    #AA
        #1. all
        DF<-dt.top
        aalist<-DF[,c("WTAA.1A","WTAA.1B","WTAA.3A")]
        # only look at the same AA sites
        sameAA<-which(DF[,paste0("WTAA.",g1)] == DF[,paste0("WTAA.",g2)])
        DF<-dt.top[sameAA,]
        
        #same mutation types
        sametype<-which(DF[,paste0("Type.",g1)] == DF[,paste0("Type.",g2)])
        DF2<-DF[sametype,]
        #different mutation types
        DF3<-DF[-sametype,]        
        
        aa1<-table(DF2[,paste0("WTAA.",g1)],DF2[,paste0("Type.",g1)])
        aa2<-table(DF2[,paste0("WTAA.",g2)],DF2[,paste0("Type.",g2)])
        
        #genome-wide proportion of AA 
        AAgenome1<-table(dt[dt$gene!="5' UTR",paste0("WTAA.",g1)],dt[dt$gene!="5' UTR",paste0("Type.",g1)] )
        AAgenome2<-table(dt[dt$gene!="5' UTR",paste0("WTAA.",g2)],dt[dt$gene!="5' UTR",paste0("Type.",g2)] )
        AAgenome.comb<-cbind(AAgenome1,AAgenome2)
        colnames(AAgenome.comb)<-c(paste0("nonsyn.",g1),paste0("stop.",g1),paste0("syn.",g1),paste0("nonsyn.",g2),paste0("stop.",g2),paste0("syn.",g2))        
        AA.genomes<-as.data.frame(AAgenome.comb)
        AA.genomes$Total<-rowSums(AA.genomes)
        AA.genomes$percent<-AA.genomes$Total/sum(AA.genomes$Total)
        
        #Calculate AA freq of nonsyn & syn sites separately for  both genotypes
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
        colnames(AA1.1)<-c("AA","Nonsyn genome-wide", "Nonsyn high Fst", "Syn genome-wide", "Syn high Fst")
        AA1.1<-melt(AA1.1)
        colnames(AA1.1)[2:3]<-c("Group","Proportion")
        
        ggplot(AA1.1, aes(x=AA, y=Proportion, fill=Group))+geom_bar(stat='identity', position='dodge')+
                labs(x="Amino acids", y="% of represented AA")+
                ggtitle(paste0("% of AA in highly different MF sites (5%) between ",g1," & ", g2))+
                scale_fill_manual("", values = cols3)+theme_classic()+
                theme(plot.title = element_text(hjust = 0.5, size=12))
        ggsave(filename=paste0("Output_all/Fst/AA.same.Top5.", fname,".pdf"),width =9, height = 4.5)
        #nonsyn only
        ggplot(AA1.1[AA1.1$Group=="Nonsyn genome-wide"|AA1.1$Group=="Nonsyn high Fst",], aes(x=AA, y=Proportion, fill=Group))+geom_bar(stat='identity', position='dodge')+
                labs(x="Amino acids", y="% of represented AA")+
                ggtitle(paste0("% of AA in highly different MF sites (5%) between ",g1," & ", g2))+
                scale_fill_manual("", values = cols3[1:2])+theme_classic()+
                theme(plot.title = element_text(hjust = 0.5, size=12))
        ggsave(filename=paste0("Output_all/Fst/AA.same.Top5.nonsyn.", fname,"..pdf"),width = 9, height = 4.5)
        #
        #show overrepresentation only
        AA1$diff.nonsyn<-(AA1$percent.nonsyn-AA1$genome.percent.nonsyn)/AA1$genome.percent.nonsyn*100
        write.csv(AA1, paste0("Output_all/Fst/AAinHighFst.Sites_", fname,".csv"))
        
        ggplot(AA1,aes(x=AA, y=diff.nonsyn))+geom_bar(stat='identity',fill='#66CCEE')+
                theme_classic()+labs(x="Amino acids", y="Over/under-represented AA (%)")+
                ggtitle(paste0("Over/Under-represented AAs of high Fst sites between ",g1," & ", g2))+
                theme(plot.title = element_text(size=12))+
                scale_fill_manual("", values = cols3)+theme_classic()
        ggsave(filename=paste0("Output_all/Fst/AA.Overrepresented.in.highFst.nonsyn.",fname,".pdf"), width = 7, height = 4.5)
        
        
        # all AA
        AA2<-data.frame(AA=rownames(AA.genomes))
        AA2.df<-data.frame(AA=rownames(aa1))
        AA2.df$highFst<-rowSums(aa1) 
        AA2.df$highF.percent<-AA2.df$highFst/sum(AA2.df$highFst)
        AA21<-merge(AA2, AA2.df, by="AA",all.x=T)
        AA21[is.na(AA21)]<-0
        
        AA21$genome.percent<-AA.genomes$percent
        AA21$diff<-(AA21$highF.percent-AA21$genome.percent)/AA21$genome.percent*100
        write.csv(AA21, paste0("Output_all/Fst/All.AAinHighF.Sites_", fname,".csv"))
        
        
        
        
        
        
        #test nonsyn distribution
       
        k=k+1
        r1<-GTest(AA1[,4:5])
        Gtest$G[k]<-r1[[1]]
        Gtest$Pvalue[k]<-r1[[3]]
        #test syn distribution
        k=k+1
        r2<-GTest(AA1[,6:7])  #not significant 
        Gtest$G[k]<-r2[[1]]
        Gtest$Pvalue[k]<-r2[[3]]
        
        n<-nrow(dt)
        n1<-nrow(dt1)
        #2. AA in nonsyn only
        DF1<-dt1.top
        aalist<-DF1[,c("WTAA.1A","WTAA.1B","WTAA.3A")]
        # only look at the same AA sites
        sameAA1<-which(DF1[,paste0("WTAA.",g1)] == DF1[,paste0("WTAA.",g2)]) 
        DF1<-DF1[sameAA1,]
        
        aa1<-table(DF1[,paste0("WTAA.",g1)],DF1[,paste0("Type.",g1)])
        aa2<-table(DF1[,paste0("WTAA.",g2)],DF1[,paste0("Type.",g2)])
        
        #genome-wide proportion of AA (nonsyn only)
        AAgenome11<-table(dt1[dt1$gene!="5' UTR",paste0("WTAA.",g1)])
        AAgenome21<-table(dt1[dt1$gene!="5' UTR",paste0("WTAA.",g2)])
        AAgenome.comb1<-cbind(AAgenome11,AAgenome21)
        colnames(AAgenome.comb1)<-c(paste0("NS.",g1),paste0("NS.",g2))        
        
        #Calculate AA freq of nonsyn sites
        Genome.comb.nonsyn1<-(AAgenome.comb1[,1]+AAgenome.comb1[,2])/(sum(AAgenome.comb1[,1])+sum(AAgenome.comb1[,2]))
        genome.AAcom1<-as.data.frame(Genome.comb.nonsyn1)
        
        #
        AA1<-data.frame(AA=rownames(genome.AAcom1))
        AA.df1<-as.data.frame(aa1)
        AA.df1<-AA.df1[,c(1,3)]
        colnames(AA.df1)<-c("AA","nonsyn")
        
        AA1<-merge(AA1, AA.df1, by="AA",all.x=T)
        AA1[is.na(AA1)]<-0
        
        AA1$genome.percent.nonsyn<-genome.AAcom1[,1]*100
        AA1$percent.nonsyn<-AA1$nonsyn/sum(AA1$nonsyn)*100
        
        AA1.1<-AA1[,c(1,3:4)]
        colnames(AA1.1)<-c("AA","Genome-wide", "High Fst sites")
        AA1.1<-melt(AA1.1)
        colnames(AA1.1)[2:3]<-c("Group","Proportion")
        
        ####
        AA1$diff.nonsyn<-(AA1$percent.nonsyn-AA1$genome.percent.nonsyn)/AA1$genome.percent.nonsyn*100
        write.csv(AA1, paste0("Output_all/Fst/NS/AAinHighFst.ns.Sites_", fname,".csv"))
        
        ggplot(AA1,aes(x=AA, y=diff.nonsyn))+geom_bar(stat='identity',fill='#66CCEE')+
                theme_classic()+labs(x="Amino acids", y="Over/under-represented AA (%)")+
                ggtitle(paste0("Over/Under-represented AAs of high Fst NS sites between ",g1," & ", g2))+
                theme(plot.title = element_text(size=12))+
                scale_fill_manual("", values = cols3)+theme_classic()
        ggsave(filename=paste0("Output_all/Fst/NS/AA.Overrepresented.in.highFst.NS.",fname,".pdf"), width = 7, height = 4.5)
        
        #test nonsyn distribution
        for (row in 1:nrow(AA1)){
                AA1$genome.ns[row]<-sum(AAgenome1[row]+AAgenome2[row])
        }
        
        
        j=j+1
        r1<-GTest(AA1[,c(2,6)])
        Gtest1$G[j]<-r1[[1]]
        Gtest1$Pvalue[j]<-r1[[3]]
        
        
        
        
        
        
        #type
        type1<-as.data.frame(table(dt.top[,paste0("Type.",g1)]))
        type2<-as.data.frame(table(dt.top[,paste0("Type.",g2)]))
        types<-merge(type1, type2, by="Var1")
        colnames(types)<-c("Type",g1,g2)
        gtype1<-as.data.frame(table(dt[,paste0("Type.",g1)]))
        gtype2<-as.data.frame(table(dt[,paste0("Type.",g2)]))
        gtypes<-merge(gtype1, gtype2, by="Var1")
        colnames(gtypes)<-c("Type", paste0("genome.",g1),paste0("genome.",g2) )
        types<-merge(types,gtypes, by="Type")
        types$difference1<-(types[,g1]-(types[,paste0("genome.",g1)]*s/n))/(types[,paste0("genome.",g1)]*s/n)*100
        types$difference2<-(types[,g2]-(types[,paste0("genome.",g2)]*s/n))/(types[,paste0("genome.",g2)]*s/n)*100
        types$difference<-rowMeans(types[,c("difference1", "difference2")])
        write.csv(types,paste0("Output_all/Fst/Types_", fname,".csv"))
        
        
        k=k+1
        types$count<-rowSums(types[,2:3])
        types$total<-rowSums(types[,4:5])
        r1<-GTest(types[,9:10])  
        Gtest$G[k]<-r1[[1]]
        Gtest$Pvalue[k]<-r1[[3]]
        
        
        
        #by codon positions
        #1. all
        Codons<-as.data.frame(table(dt.top$codon))
        gCodons<-as.data.frame(table(dt$codon))
        Codons<-merge(Codons,gCodons, by="Var1")
        colnames(Codons)<-c("codon","highFst", "genome")
        
        Codons$difference<-(Codons$highFst-(Codons$genome*s/n))/(Codons$genome*s/n)*100
        write.csv(Codons, paste0("Output_all/Fst/Codons_", fname,".csv"))
        
        k=k+1
        r1<-GTest(Codons[,2:3])  
        Gtest$G[k]<-r1[[1]]
        Gtest$Pvalue[k]<-r1[[3]]
        
        #2 nonsyn only
        Codons1<-as.data.frame(table(dt1.top$codon))
        gCodons1<-as.data.frame(table(dt1$codon))
        Codons1<-merge(Codons1,gCodons1, by="Var1")
        colnames(Codons1)<-c("codon","highFst", "genome")
        
        Codons1$difference<-(Codons1$highFst-(Codons1$genome*s/n1))/(Codons1$genome*s/n1)*100
        write.csv(Codons1, paste0("Output_all/Fst/NS/Codons.ns_", fname,".csv"))
        
        j=j+1
        r1<-GTest(Codons1[,2:3])  
        Gtest1$G[j]<-r1[[1]]
        Gtest1$Pvalue[j]<-r1[[3]]
        
        
        
        
        #Nucleotide
        Nt<-as.data.frame(table(dt.top[,paste0("ref.",g1)]))
        gNt<-as.data.frame(table(dt[,paste0("ref.",g1)]))
        Nt<-merge(Nt,gNt, by="Var1")
        colnames(Nt)<-c("nucleotide","highFst", "genome")
        Nt$difference<-(Nt$highFst-(Nt$genome*s/n))/(Nt$genome*s/n)*100
        
        write.csv(Nt,paste0("Output_all/Fst/NT_", fname,".csv"))
        
        k=k+1
        r1<-GTest(Nt[,2:3])  
        Gtest$G[k]<-r1[[1]]
        Gtest$Pvalue[k]<-r1[[3]]
        
        #NS only
        Nt1<-as.data.frame(table(dt1.top[,paste0("ref.",g1)]))
        gNt1<-as.data.frame(table(dt1[,paste0("ref.",g1)]))
        Nt1<-merge(Nt1,gNt1, by="Var1")
        colnames(Nt1)<-c("nucleotide","highFst", "genome")
        Nt1$difference<-(Nt1$highFst-(Nt1$genome*s/n1))/(Nt1$genome*s/n1)*100
        
        write.csv(Nt1,paste0("Output_all/Fst/NS/NT.ns_", fname,".csv"))
        
        j=j+1
        r1<-GTest(Nt1[,2:3])  
        Gtest1$G[j]<-r1[[1]]
        Gtest1$Pvalue[j]<-r1[[3]]
        
        
        #CpG Sites
        #1. all
        cpg1<-as.data.frame(table(dt.top[,paste0("makesCpG.",g1)]))
        cpg2<-as.data.frame(table(dt.top[,paste0("makesCpG.",g2)]))
        cpg<-merge(cpg1,cpg2,by="Var1")
        colnames(cpg)<-c("CpG",g1,g2)
        Gcpg1<-as.data.frame(table(dt[,paste0("makesCpG.",g1)]))
        Gcpg2<-as.data.frame(table(dt[,paste0("makesCpG.",g2)]))
        Gcpg<-merge(Gcpg1,Gcpg2,by="Var1")
        colnames(Gcpg)<-c("CpG",g1,g2)
        cpg<-merge(cpg,Gcpg, by="CpG")
        cpg$diff1<-(cpg[,paste0(g1,".x")]-(cpg[,paste0(g1,".y")]*s/n))/(cpg[,paste0(g1,".y")]*s/n)*100
        cpg$diff2<-(cpg[,paste0(g2,".x")]-(cpg[,paste0(g2,".y")]*s/n))/(cpg[,paste0(g2,".y")]*s/n)*100
        cpg$difference<-rowMeans(cpg[,c("diff1","diff2")])     
        write.csv(cpg,paste0("Output_all/Fst/CpG_", fname,".csv"))
        
        k=k+1
        cpg$count<-rowSums(cpg[,2:3])
        cpg$total<-rowSums(cpg[,4:5])
        r1<-GTest(cpg[,9:10])  
        Gtest$G[k]<-r1[[1]]
        Gtest$Pvalue[k]<-r1[[3]]
        
        write.csv(Gtest, paste0("Output_all/Fst/Gtest.results.", fname,".csv"))
        
        #2. NS only
        cpg21<-as.data.frame(table(dt1.top[,paste0("makesCpG.",g2)]))
        cpg11<-as.data.frame(table(dt1.top[,paste0("makesCpG.",g1)]))
        Cpg<-merge(cpg11,cpg21,by="Var1")
        colnames(Cpg)<-c("CpG",g1,g2)
        Gcpg11<-as.data.frame(table(dt1[,paste0("makesCpG.",g1)]))
        Gcpg21<-as.data.frame(table(dt1[,paste0("makesCpG.",g2)]))
        GCpg<-merge(Gcpg11,Gcpg21,by="Var1")
        colnames(GCpg)<-c("CpG",g1,g2)
        Cpg<-merge(Cpg,GCpg, by="CpG")
        Cpg$diff1<-(Cpg[,paste0(g1,".x")]-(Cpg[,paste0(g1,".y")]*s/n1))/(Cpg[,paste0(g1,".y")]*s/n1)*100
        Cpg$diff2<-(Cpg[,paste0(g2,".x")]-(Cpg[,paste0(g2,".y")]*s/n1))/(Cpg[,paste0(g2,".y")]*s/n1)*100
        Cpg$difference<-rowMeans(Cpg[,c("diff1","diff2")])     
        write.csv(Cpg,paste0("Output_all/Fst/NS/CpG.ns_", fname,".csv"))
        
        j=j+1
        Cpg$count<-rowSums(Cpg[,2:3])
        Cpg$total<-rowSums(Cpg[,4:5])
        r1<-GTest(Cpg[,9:10])  
        Gtest1$G[j]<-r1[[1]]
        Gtest1$Pvalue[j]<-r1[[3]]
        
        write.csv(Gtest1, paste0("Output_all/Fst/NS/Gtest.results.NSonly.", fname,".csv"))
        
        
}







#### create figures

aafile<-list.files("Output_all/Fst/", pattern = "AAinHighFst")
aadiff<-list()
for (i in 1:length(aafile)){
        dt<-read.csv(paste0("Output_all/Fst/", aafile[i]), stringsAsFactors = F, row.names = 1)
        aadiff[[i]]<-dt[,c("AA","diff.nonsyn")]
        names(aadiff)[i]<-substr(aafile[i], start=19, stop = 23)
} 

for (i in 1:length(aadiff)) {
        colnames(aadiff[[i]])<-c("AA",paste0(names(aadiff[i])))
}

AAdiff<-aadiff%>% purrr::reduce(full_join, by='AA')
AAdiffm<-melt(AAdiff)
ggplot(AAdiffm,aes(x=AA, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_classic()+labs(x="Amino acid", y="Over/under-represented nonsyn AA (%)")+
        #ggtitle(paste0("Over/Under-represented nonsyn AAs of highly different MF sites (5%) between ",g1," & ", g2))+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols4[c(2,1,3)])+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/Fst/Overrepresented.AAs.NS.pdf"), width = 7, height = 4.5)

## AA plot including syn sites
aaf<-list.files("Output_all/Fst/", pattern = "All.AAinHighF.Sites")
aadiff2<-list()
for (i in 1:length(aafile)){
        dt<-read.csv(paste0("Output_all/Fst/", aaf[i]), stringsAsFactors = F, row.names = 1)
        aadiff2[[i]]<-dt[,c("AA","diff")]
        names(aadiff2)[i]<-substr(aaf[i], start=21, stop = 25)
} 

for (i in 1:length(aadiff2)) {
        colnames(aadiff2[[i]])<-c("AA",paste0(names(aadiff[i])))
}

AAdiff2<-aadiff2%>% purrr::reduce(full_join, by='AA')
AAdiff2m<-melt(AAdiff2)
ggplot(AAdiff2m,aes(x=AA, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_classic()+labs(x="Amino acid", y="Over/under-represented nonsyn AA (%)")+
        #ggtitle(paste0("Over/Under-represented nonsyn AAs of highly different MF sites (5%) between ",g1," & ", g2))+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols4[c(2,1,3)])+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/Fst/Overrepresented.AAs.ALL.pdf"), width = 7, height = 4.5)

## Genes
genefile<-list.files("Output_all/Fst/", pattern="Top5.Fst.SummarybyGene")
geneprop<-list()
for (i in 1:length(genefile)){
        dt<-read.csv(paste0("Output_all/Fst/", genefile[i]), stringsAsFactors = F, row.names = 1)
        geneprop[[i]]<-dt[,c("Gene","Difference")]
        names(geneprop)[i]<-substr(genefile[i], start=24, stop = 28)
} 
for (i in 1:length(geneprop)) {
        colnames(geneprop[[i]])<-c("Gene",paste0(names(geneprop[i])))
}

GeneProp<-geneprop%>% purrr::reduce(full_join, by='Gene')
GeneProp$Gene<-factor(GeneProp$Gene, levels=c("5' UTR","Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

GenePropm<-melt(GeneProp)
ggplot(GenePropm,aes(x=Gene, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_classic()+labs(x="", y="Over/under-represented sites per gene (%)")+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols4[c(2,1,3)])+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/Fst/Overrepresented.Gene.pdf"), width = 7, height = 4.5)

## codon

codonfile<-list.files("Output_all/Fst/", pattern = glob2rx("*Codons*.csv*"))
codondiff<-list()
for (i in 1:length(codonfile)){
        dt<-read.csv(paste0("Output_all/Fst/", codonfile[i]), stringsAsFactors = F, row.names = 1)
        codondiff[[i]]<-dt[,c("codon","difference")]
        names(codondiff)[i]<-substr(codonfile[i], start=8, stop = 12)
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
ggsave(filename=paste0("Output_all/Fst/Overrepresented.Codons.pdf"), width = 5, height = 3.2)



## NT

ntfile<-list.files("Output_all/Fst/", pattern = glob2rx("*NT*.csv*"))
ntdiff<-list()
for (i in 1:length(ntfile)){
        dt<-read.csv(paste0("Output_all/Fst/", ntfile[i]), stringsAsFactors = F, row.names = 1)
        ntdiff[[i]]<-dt[,c("nucleotide","difference")]
        names(ntdiff)[i]<-substr(ntfile[i], start=4, stop = 8)
} 

for (i in 1:length(ntdiff)) {
        colnames(ntdiff[[i]])<-c("nt",paste0(names(ntdiff[i])))
}

ntdiff<-ntdiff%>% purrr::reduce(full_join, by='nt')
ntdiff$nt<-c("A","C","G","T")
ntdiffm<-melt(ntdiff)
ggplot(ntdiffm,aes(x=nt, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_classic()+labs(x="", y="Over/under-represented nucleotide (%)")+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols4[c(2,1,3)])+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/Fst/Overrepresented.NT.pdf"), width = 6, height = 4)



## Type

tpfile<-list.files("Output_all/Fst/", pattern =glob2rx("*Types*.csv*"))
tpdiff<-list()
for (i in 1:length(tpfile)){
        dt<-read.csv(paste0("Output_all/Fst/", tpfile[i]), stringsAsFactors = F, row.names = 1)
        tpdiff[[i]]<-dt[,c("Type","difference")]
        names(tpdiff)[i]<-substr(tpfile[i], start=7, stop = 11)
} 

for (i in 1:length(tpdiff)) {
        colnames(tpdiff[[i]])<-c("Type",paste0(names(tpdiff[i])))
}

tpdiff<-tpdiff%>% purrr::reduce(full_join, by='Type')
tpdiffm<-melt(tpdiff)
ggplot(tpdiffm,aes(x=Type, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_classic()+labs(x="Mutation type", y="Over/under-representation (%)")+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols4[c(2,1,3)])+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/Fst/Overrepresented.Type.pdf"), width = 5, height = 3.3)

##CpG
cpgfile<-list.files("Output_all/Fst/", pattern = glob2rx("*CpG*.csv*"))
cpgdiff<-list()
for (i in 1:length(cpgfile)){
        dt<-read.csv(paste0("Output_all/Fst/", cpgfile[i]), stringsAsFactors = F, row.names = 1)
        cpgdiff[[i]]<-dt[,c("CpG","difference")]
        names(cpgdiff)[i]<-substr(cpgfile[i], start=5, stop = 9)
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
ggsave(filename=paste0("Output_all/Fst/Overrepresented.CpG.pdf"), width = 4.5, height = 3)




########################################
##### Create figures for NS only data



aafile<-list.files("Output_all/Fst/NS/", pattern = "AAinHighFst")
aadiff<-list()
for (i in 1:length(aafile)){
        dt<-read.csv(paste0("Output_all/Fst/NS/", aafile[i]), stringsAsFactors = F, row.names = 1)
        aadiff[[i]]<-dt[,c("AA","diff.nonsyn")]
        names(aadiff)[i]<-substr(aafile[i], start=22, stop = 26)
} 

for (i in 1:length(aadiff)) {
        colnames(aadiff[[i]])<-c("AA",paste0(names(aadiff[i])))
}

AAdiff<-aadiff%>% purrr::reduce(full_join, by='AA')
AAdiffm<-melt(AAdiff)
ggplot(AAdiffm,aes(x=AA, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_classic()+labs(x="Amino acid", y="Over/under-represented nonsyn AA (%)")+
        #ggtitle(paste0("Over/Under-represented nonsyn AAs of highly different MF sites (5%) between ",g1," & ", g2))+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols4[c(2,1,3)])+theme_classic()+
        labs(color="")+
        geom_hline(yintercept = 0,  color = "gray50", size=.5)
ggsave(filename=paste0("Output_all/Fst/NS/Overrepresented.AAs.pdf"), width = 7, height = 4.5)


genefile<-list.files("Output_all/Fst/NS/", pattern="Fst.SummarybyGene")
geneprop<-list()
for (i in 1:length(genefile)){
        dt<-read.csv(paste0("Output_all/Fst/NS/", genefile[i]), stringsAsFactors = F, row.names = 1)
        geneprop[[i]]<-dt[,c("Gene","Difference")]
        names(geneprop)[i]<-substr(genefile[i], start=27, stop = 31)
} 
for (i in 1:length(geneprop)) {
        colnames(geneprop[[i]])<-c("Gene",paste0(names(geneprop[i])))
}

GeneProp<-geneprop%>% purrr::reduce(full_join, by='Gene')
GenePropm<-melt(GeneProp)
ggplot(GenePropm,aes(x=Gene, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_classic()+labs(x="", y="Over/under-represented sites per gene (%)")+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols4[c(2,1,3)])+theme_classic()+
        labs(color="")+
        geom_hline(yintercept = 0,  color = "gray50", size=.5)
ggsave(filename=paste0("Output_all/Fst/NS/Overrepresented.Gene.pdf"), width = 7, height = 4.5)

## codon

codonfile<-list.files("Output_all/Fst/NS/", pattern = glob2rx("*Codons*.csv*"))
codondiff<-list()
for (i in 1:length(codonfile)){
        dt<-read.csv(paste0("Output_all/Fst/NS/", codonfile[i]), stringsAsFactors = F, row.names = 1)
        codondiff[[i]]<-dt[,c("codon","difference")]
        names(codondiff)[i]<-substr(codonfile[i], start=11, stop = 15)
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
        labs(color="")+
        geom_hline(yintercept = 0,  color = "gray50", size=.5)
ggsave(filename=paste0("Output_all/Fst/NS/Overrepresented.Codons.pdf"), width = 5, height = 3.2)



## NT

ntfile<-list.files("Output_all/Fst/NS/", pattern = glob2rx("*NT*.csv*"))
ntdiff<-list()
for (i in 1:length(ntfile)){
        dt<-read.csv(paste0("Output_all/Fst/NS/", ntfile[i]), stringsAsFactors = F, row.names = 1)
        ntdiff[[i]]<-dt[,c("nucleotide","difference")]
        names(ntdiff)[i]<-substr(ntfile[i], start=7, stop = 11)
} 

for (i in 1:length(ntdiff)) {
        colnames(ntdiff[[i]])<-c("nt",paste0(names(ntdiff[i])))
}

ntdiff<-ntdiff%>% purrr::reduce(full_join, by='nt')
ntdiff$nt<-c("A","C","G","T")
ntdiffm<-melt(ntdiff)
ggplot(ntdiffm,aes(x=nt, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_classic()+labs(x="", y="Over/under-represented nucleotide (%)")+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = cols4[c(2,1,3)])+theme_classic()+
        labs(color="")+
        geom_hline(yintercept = 0,  color = "gray50", size=.5)
ggsave(filename=paste0("Output_all/Fst/NS/Overrepresented.NT.pdf"), width = 7, height = 4.5)


##CpG
cpgfile<-list.files("Output_all/Fst/NS/", pattern = glob2rx("*CpG*.csv*"))
cpgdiff<-list()
for (i in 1:length(cpgfile)){
        dt<-read.csv(paste0("Output_all/Fst/NS/", cpgfile[i]), stringsAsFactors = F, row.names = 1)
        cpgdiff[[i]]<-dt[,c("CpG","difference")]
        names(cpgdiff)[i]<-substr(cpgfile[i], start=8, stop = 15)
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
        labs(color="")+
        geom_hline(yintercept = 0,  color = "gray50", size=.5)
ggsave(filename=paste0("Output_all/Fst/NS/Overrepresented.CpG.all.top5.pdf"), width = 4.5, height = 3)




