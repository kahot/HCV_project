library(stringr)
library(tidyverse)
library(zoo)
library(reshape2)
library(ggplot2)
library(DescTools)
library(colorspace)
source("Rscripts/baseRscript.R")

#dir.create("Output_all/Fst/")
colors2<-qualitative_hcl(6, palette="Dark3")
col2_light<-qualitative_hcl(6, palette="Set3")
div.colors<-c(colors2[1],col2_light[1],colors2[3],col2_light[3],colors2[5],col2_light[5])
color.genes<-qualitative_hcl(11, palette="Dark3")

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

#merged.meta<-read.csv("Output_all/merged.metadata.csv", stringsAsFactors = F, row.names = 1)
#merged.meta$gene[merged.meta$gene=="NS1(P7)"]<-"NS1"

Summary<-read.csv("Output_all/Unfiltered/Ts.MinorVariant_Mean_3genotypes.csv",stringsAsFactors = F, row.names = 1)
Summary$gene[Summary$gene=="NS1(P7)"]<-"NS1"

for (g in 1:3){        
        dt<-FstList[[g]]
        fname<-substr(names(FstList)[g], start=5, stop = 9)
        g1<-substr(names(FstList)[g], start=5, stop = 6)
        g2<-substr(names(FstList)[g], start=8, stop = 9)
        
        #attached the metadata
        dt<-merge(dt, Summary, by="merged.pos")
        dt<-dt[dt$Type.1A!="stop",]
        dt<-dt[dt$Type.1B!="stop",]
        dt<-dt[dt$Type.3A!="stop",]
        
        #top5%
        s<-as.integer(nrow(dt)*0.05)
        dt$diff<-dt[,paste0("mean.",g1)]-dt[,paste0("mean.",g2)]
        
        #select sites with MF g1> g2 (costly sites in g2)
        dt1<-dt[dt$diff>0,]
        dt1.top<-dt1[order(dt1$Fst,decreasing = T,na.last = T),]
        dt1.top<-dt1.top[c(1:s),]
        
        #select sites with MF g1 < g2 (costly sites in g1)
        dt2<-dt[dt$diff<0,]
        dt2.top<-dt2[order(dt2$Fst,decreasing = T,na.last = T),]
        dt2.top<-dt2.top[c(1:s),]
        
        #volcano plot like
        plot(dt$diff, dt$Fst, pch=16, cex=.5, col='grey', xlim=c(-0.08,0.08), xlab="Difference in mutation frequency (1A -1B)",
             ylab="Fst")
        points(dt1.top$diff, dt1.top$Fst, pch=16, cex=.5,col=colors2[1])
        points(dt2.top$diff, dt2.top$Fst, pch=16, cex=.5, col=colors2[3])
        abline(v=0, col='grey50', lty=3,lwd=1.3)
        plot(dt$diff, dt$Fst, pch=16, cex=.5, col='grey', xlim=c(-0.08,0.08), xlab="Difference in mutation frequency (1A -1B)",
             ylab="Fst")
        points(dt1.top$diff[dt1.top$gene=="HVR1"],dt1.top$Fst[dt1.top$gene=="HVR1"],pch=16, cex=.5,col=color.genes[1])
        
    #check proportion of sites belonging to each gene
        #1. all
        highFst1<-data.frame("Gene"= paste0(genes$Gene[1:12]), stringsAsFactors = FALSE)
        highFst2<-data.frame("Gene"= paste0(genes$Gene[1:12]), stringsAsFactors = FALSE)
        
        Dat<-dt[!is.na(dt$Fst),]
        for (i in 1:12){
                n<-Dat[Dat$gene==genes$Gene[i],]
                df1<-dt1.top[dt1.top$gene==genes$Gene[i],]
                df2<-dt2.top[dt2.top$gene==genes$Gene[i],]
                highFst1$Counts[i]<-nrow(df1)
                highFst1$Total[i]<-nrow(n)
                highFst1$Expected[i]<-nrow(n)*s/nrow(Dat)
                highFst1$Difference[i] <-(highFst1$Counts[i]-highFst1$Expected[i])/highFst1$Expected[i]*100
        
                highFst2$Counts[i]<-nrow(df2)
                highFst2$Total[i]<-nrow(n)
                highFst2$Expected[i]<-nrow(n)*s/nrow(Dat)
                highFst2$Difference[i] <-(highFst2$Counts[i]-highFst2$Expected[i])/highFst2$Expected[i]*100
        }
        
        write.csv(highFst1,paste0("Output_all/Fst2/HighFst.byGene.",fname,".costly.in.",g2,".csv"))
        write.csv(highFst2,paste0("Output_all/Fst2/HighFst.byGene.",fname,".costly.in.",g1,".csv"))
        
        if (g1=="1A"& g2=="3A") {barcolor<- colors2[1]; barcolor2<- colors2[5]}
        if (g1=="1A"& g2=="1B") {barcolor<- colors2[1]; barcolor2<- colors2[3]}
        if (g1=="1B"& g2=="3A") {barcolor<- colors2[3]; barcolor2<- colors2[5]}
        
        
        ggplot(highFst1,aes(x=Gene,y=Difference))+geom_bar(stat='identity', fill=barcolor2,width=0.8)+
                labs(x="Genes",y=paste0("Over/under-representation (%) \n(costly sites in ",g2, " over ",g1,")"))+
                theme_classic()+
                geom_hline(yintercept = 0,  color = "gray50", size=.5)
        ggsave(filename=paste0("Output_all/Fst2/Top5.Fst.by.gene.",fname,".costly.",g2,".pdf"), width = 5.7, height = 4.5)
        
        ggplot(highFst2,aes(x=Gene,y=Difference))+geom_bar(stat='identity', fill=barcolor,width=0.8)+
                labs(x="Genes",y=paste0("Over/under-representation (%) \n(costly sites in ",g1, " over ",g1,")"))+
                theme_classic()+
                geom_hline(yintercept = 0,  color = "gray50", size=.5)
        ggsave(filename=paste0("Output_all/Fst2/Top5.Fst.by.gene.",fname,".costly.",g1,".pdf"), width = 5.7, height = 4.5)
        
        #Test costly sites in g2
        Gtest1<-data.frame(test=c("gene","AA-ns","AA-syn","Type", "codon","nt","CpG"))
        k=1
        r1<-GTest(highFst1[,2:3])
        Gtest1$G[k]<-r1[[1]]
        Gtest1$Pvalue[k]<-r1[[3]]
        
        #Test costly sites in g1
        Gtest2<-data.frame(test=c("gene","AA-ns","AA-syn","Type", "codon","nt","CpG"))
        k=1
        r2<-GTest(highFst2[,2:3])
        Gtest2$G[k]<-r2[[1]]
        Gtest2$Pvalue[k]<-r2[[3]]
        
    #AA
        #1. all
        DF1<-dt1.top
        DF2<-dt2.top
        aalist1<-DF1[,c("WTAA.1A","WTAA.1B","WTAA.3A")]
        aalist2<-DF2[,c("WTAA.1A","WTAA.1B","WTAA.3A")]
        
        # only look at the same AA sites
        sameAA1<-which(DF1[,paste0("WTAA.",g1)] == DF1[,paste0("WTAA.",g2)])
        sameAA2<-which(DF2[,paste0("WTAA.",g1)] == DF2[,paste0("WTAA.",g2)])
        
        DF1<-dt1.top[sameAA1,]
        DF2<-dt2.top[sameAA2,]
        
        #same mutation types
        sametype1<-which(DF1[,paste0("Type.",g1)] == DF1[,paste0("Type.",g2)])
        DF1<-DF1[sametype1,]
        
        sametype2<-which(DF2[,paste0("Type.",g1)] == DF2[,paste0("Type.",g2)])
        DF2<-DF2[sametype2,]
        
        aa1<-table(DF1[,paste0("WTAA.",g1)],DF1[,paste0("Type.",g1)])
        aa2<-table(DF2[,paste0("WTAA.",g2)],DF2[,paste0("Type.",g2)])
        
        #genome-wide proportion of AA 
        AAgenome1<-table(dt[dt$gene!="5' UTR",paste0("WTAA.",g1)],dt[dt$gene!="5' UTR",paste0("Type.",g1)] )
        AAgenome2<-table(dt[dt$gene!="5' UTR",paste0("WTAA.",g2)],dt[dt$gene!="5' UTR",paste0("Type.",g2)] )
        AAgenome.comb<-cbind(AAgenome1,AAgenome2)
        colnames(AAgenome.comb)<-c(paste0("nonsyn.",g1),paste0("syn.",g1),paste0("nonsyn.",g2),paste0("syn.",g2))        
        AA.genomes<-as.data.frame(AAgenome.comb)
        AA.genomes$Total<-rowSums(AA.genomes)
        AA.genomes$percent<-AA.genomes$Total/sum(AA.genomes$Total)
        
        #Calculate AA freq of nonsyn & syn sites separately for  both genotypes
        Genome.comb.nonsyn<-(AAgenome.comb[,1]+AAgenome.comb[,3])/(sum(AAgenome.comb[,1])+sum(AAgenome.comb[,3]))
        Genome.comb.syn   <-(AAgenome.comb[,2]+AAgenome.comb[,4])/(sum(AAgenome.comb[,2])+sum(AAgenome.comb[,4]))
        genome.AAcom<-as.data.frame(cbind(Genome.comb.nonsyn,Genome.comb.syn))
        
        for (i in 1:2){
                AA<-data.frame(AA=rownames(genome.AAcom))
                aa<-get(paste0("aa",i))
                AA.df<-data.frame(AA=rownames(aa))
                AA.df$all<-rowSums(aa)
                AA.df$nonsyn<-aa[,"nonsyn"]
                AA.df$syn<-aa[,"syn"]
                AA<-merge(AA, AA.df, by="AA",all.x=T)
                AA[is.na(AA)]<-0
                
                AA$genome.percent<-AA.genomes$percent*100
                AA$Fst1.percent<-AA$all/sum(AA1$all)*100
                AA$genome.percent.nonsyn<-genome.AAcom[,1]*100
                AA$percent.nonsyn<-AA$nonsyn/sum(AA$nonsyn)*100
                AA$genome.percent.syn<-genome.AAcom[,2]*100
                AA$percent.syn<-AA$syn/sum(AA$syn)*100
                
                AA1.1<-AA[,c(1,5:10)]
                colnames(AA1.1)<-c("AA","Genome-wide", "High Fst", "Nonsyn genome-wide", "Nonsyn high Fst", "Syn genome-wide", "Syn high Fst")
                AA1.1<-melt(AA1.1)
                colnames(AA1.1)[2:3]<-c("Group","Proportion")
                
                ggplot(AA1.1, aes(x=AA, y=Proportion, fill=Group))+geom_bar(stat='identity', position='dodge')+
                        labs(x="Amino acids", y="% of represented AA")+
                        scale_fill_manual("", values = div.colors)+theme_classic()+
                        theme(plot.title = element_text(hjust = 0.5, size=12))+
                        geom_vline(xintercept = c(1:19)+0.5, color = "gray60", size=.4)
                if (i==1) genotype<-g2
                if (i==2) genotype<-g1
                
                ggsave(filename=paste0("Output_all/Fst2/AA.", fname,".costly.in.",genotype,".pdf"),width =11, height = 4.5)

                #show over/under representation of nonsyn AA
                AA$diff.nonsyn<-(AA$percent.nonsyn-AA$genome.percent.nonsyn)/AA$genome.percent.nonsyn*100
                write.csv(AA, paste0("Output_all/Fst2/AAinHighFst.Sites_", fname,".costly.in.",genotype,".csv"))
                
                ggplot(AA,aes(x=AA, y=diff.nonsyn))+geom_bar(stat='identity',fill='#66CCEE')+
                        theme_classic()+labs(x="Amino acids", y="Over/under-represented AA (%)")+
                        ggtitle(paste0("Over/Under-represented nonsyn AAs (costly sites in ",genotype, " over ",genotype1,")") )+
                        theme(plot.title = element_text(size=12))+
                        scale_fill_manual("", values = cols3)+theme_classic()
                
                ggsave(filename=paste0("Output_all/Fst2/AA.Overrepresented.nonsyn.", fname,".costly.in.", genotype,".pdf"), width = 7, height = 4.5)
                
                tname<-paste0("AA",i)
                assign(tname,AA)
        }
        
      
        #test nonsyn distribution
        k=k+1
        r1<-GTest(AA1[,c("genome.percent.nonsyn", "percent.nonsyn")])
        Gtest1$G[k]<-r1[[1]]
        Gtest1$Pvalue[k]<-r1[[3]]
        r2<-GTest(AA2[,c("genome.percent.nonsyn", "percent.nonsyn")])
        Gtest2$G[k]<-r2[[1]]
        Gtest2$Pvalue[k]<-r2[[3]]
        
        
        #test syn distribution
        k=k+1
        r1<-GTest(AA1[,c("genome.percent.syn", "percent.syn")])
        Gtest1$G[k]<-r1[[1]]
        Gtest1$Pvalue[k]<-r1[[3]]
        r2<-GTest(AA2[,c("genome.percent.syn", "percent.syn")])
        Gtest2$G[k]<-r2[[1]]
        Gtest2$Pvalue[k]<-r2[[3]]
        
        GTest(AA1[,c("genome.percent", "Fst1.percent")]) #ns
        GTest(AA2[,c("genome.percent", "Fst1.percent")]) #ns
  
        #type
        #costly sites in g2
        type1<-as.data.frame(table(dt1.top[,paste0("Type.",g2)]))
        #costly sites in g1
        type2<-as.data.frame(table(dt2.top[,paste0("Type.",g1)]))
        types<-merge(type1, type2, by="Var1")
        colnames(types)<-c("Type",g2,g1)
        
        gtype1<-as.data.frame(table(dt[,paste0("Type.",g2)]))
        gtype2<-as.data.frame(table(dt[,paste0("Type.",g1)]))
        gtypes<-merge(gtype1, gtype2, by="Var1")
        colnames(gtypes)<-c("Type", paste0("genome.",g2),paste0("genome.",g1) )
        
        n<-nrow(dt)
        types<-merge(types,gtypes, by="Type")
        types[,paste0("diff.costly.",g1)]<-(types[,g1]-(types[,paste0("genome.",g1)]*s/n))/(types[,paste0("genome.",g1)]*s/n)*100
        types[,paste0("diff.costly.",g2)]<-(types[,g2]-(types[,paste0("genome.",g2)]*s/n))/(types[,paste0("genome.",g2)]*s/n)*100
        #types$difference<-rowMeans(types[,c("difference1", "difference2")])
        write.csv(types,paste0("Output_all/Fst2/Types_", fname,".csv"))
        
        
        k=k+1
        r1<-GTest(types[,c(2, 4)])  
        Gtest1$G[k]<-r1[[1]]
        Gtest1$Pvalue[k]<-r1[[3]]
        r2<-GTest(types[,c(3, 5)])  
        Gtest2$G[k]<-r2[[1]]
        Gtest2$Pvalue[k]<-r2[[3]]
        
        
        
        #by codon positions
        #1. all
        Codons1<-as.data.frame(table(dt1.top$codon))
        Codons2<-as.data.frame(table(dt2.top$codon))
        gCodons<-as.data.frame(table(dt$codon))
        
        Codons1<-merge(Codons1,gCodons, by="Var1")
        Codons2<-merge(Codons2,gCodons, by="Var1")
        colnames(Codons1)<-c("codon","highFst", "genome")
        colnames(Codons2)<-c("codon","highFst", "genome")
        
        Codons1$difference<-(Codons1$highFst-(Codons1$genome*s/n))/(Codons1$genome*s/n)*100
        write.csv(Codons1, paste0("Output_all/Fst2/Codons_", fname,".costly.in.",g2,".csv"))
        Codons2$difference<-(Codons2$highFst-(Codons2$genome*s/n))/(Codons2$genome*s/n)*100
        write.csv(Codons2, paste0("Output_all/Fst2/Codons_", fname,".costly.in.",g1,".csv"))
        
        k=k+1
        r1<-GTest(Codons1[,2:3])  
        Gtest1$G[k]<-r1[[1]]
        Gtest1$Pvalue[k]<-r1[[3]]
        r2<-GTest(Codons2[,2:3])  
        Gtest2$G[k]<-r2[[1]]
        Gtest2$Pvalue[k]<-r2[[3]]
        
       
        
        #Nucleotide
        Nt1<-as.data.frame(table(dt1.top[,paste0("ref.",g1)]))
        Nt2<-as.data.frame(table(dt2.top[,paste0("ref.",g1)]))
        gNt<-as.data.frame(table(dt[,paste0("ref.",g1)]))
        Nt1<-merge(Nt1,gNt, by="Var1")
        Nt2<-merge(Nt2,gNt, by="Var1")
        colnames(Nt1)<-c("nucleotide","highFst", "genome")
        colnames(Nt2)<-c("nucleotide","highFst", "genome")
        
        Nt1$difference<-(Nt1$highFst-(Nt1$genome*s/n))/(Nt1$genome*s/n)*100
        Nt2$difference<-(Nt2$highFst-(Nt2$genome*s/n))/(Nt2$genome*s/n)*100
        
        write.csv(Nt1,paste0("Output_all/Fst2/NT_", fname,".costly.in.",g2,".csv"))
        write.csv(Nt2,paste0("Output_all/Fst2/NT_", fname,".costly.in.",g1,".csv"))
        
        k=k+1
        r1<-GTest(Nt1[,2:3])  
        Gtest1$G[k]<-r1[[1]]
        Gtest1$Pvalue[k]<-r1[[3]]
        r2<-GTest(Nt2[,2:3])  
        Gtest2$G[k]<-r2[[1]]
        Gtest2$Pvalue[k]<-r2[[3]]
        
        
        #CpG Sites
        #1. all
        cpg1<-as.data.frame(table(dt1.top[,paste0("makesCpG.",g2)]))
        cpg2<-as.data.frame(table(dt2.top[,paste0("makesCpG.",g1)]))
        Gcpg1<-as.data.frame(table(dt[,paste0("makesCpG.",g2)]))
        Gcpg2<-as.data.frame(table(dt[,paste0("makesCpG.",g1)]))
        
        cpg1<-merge(cpg1,Gcpg1,by="Var1")
        colnames(cpg1)<-c("CpG",g2, "genome")
        cpg2<-merge(cpg2,Gcpg2,by="Var1")
        colnames(cpg2)<-c("CpG",g1, "genome")
        
        cpg1$diff<-(cpg1[,g2]-(cpg1[,"genome"]*s/n))/(cpg1[,"genome"]*s/n)*100
        cpg2$diff<-(cpg2[,g1]-(cpg2[,"genome"]*s/n))/(cpg2[,"genome"]*s/n)*100
        write.csv(cpg1,paste0("Output_all/Fst2/CpG_", fname,".costly.in.",g2,".csv"))
        write.csv(cpg2,paste0("Output_all/Fst2/CpG_", fname,".costly.in.",g1,".csv"))
        
        
        k=k+1
        r1<-GTest(cpg1[,2:3])  
        Gtest1$G[k]<-r1[[1]]
        Gtest1$Pvalue[k]<-r1[[3]]
        r2<-GTest(cpg2[,2:3])
        Gtest2$G[k]<-r2[[1]]
        Gtest2$Pvalue[k]<-r2[[3]]
        
        write.csv(Gtest1, paste0("Output_all/Fst2/Gtest.results.", fname,".costly.in.",g2,".csv"))
        write.csv(Gtest2, paste0("Output_all/Fst2/Gtest.results.", fname,".costly.in.",g1,".csv"))
        
}







#### create figures

aafile<-list.files("Output_all/Fst2/", pattern = "AAinHighFst")

#nonsyn
aadiff<-list()
genonames<-c("1A-1B","1A-3A","1B-3A")
for (i in 1:length(aafile)){
        dt<-read.csv(paste0("Output_all/Fst2/", aafile[i]), stringsAsFactors = F, row.names = 1)
        aadiff[[i]]<-dt[,c("AA","diff.nonsyn")]
        names(aadiff)[i]<-substr(aafile[i], start=19, stop = 36)
} 

for (i in 1:length(aadiff)) {
        colnames(aadiff[[i]])<-c("AA",paste0(names(aadiff[i])))
}

AAdiff<-aadiff%>% purrr::reduce(full_join, by='AA')
AAdiffm<-melt(AAdiff)
ggplot(AAdiffm,aes(x=AA, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_classic()+labs(x="Amino acid", y="Over/under-represented nonsyn AA (%)")+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = colors2)+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/Fst2/Overrepresented.AAs.NS.pdf"), width = 7, height = 4.5)


## AA plot including syn sites
aadiff2<-list()
for (i in 1:length(aafile)){
        dt<-read.csv(paste0("Output_all/Fst2/", aafile[i]), stringsAsFactors = F, row.names = 1)
        dt$diff2<-(dt$genome.percent-dt$Fst1.percent)/dt$genome.percent*100
        aadiff2[[i]]<-dt[,c("AA","diff2")]
        names(aadiff)[i]<-substr(aafile[i], start=19, stop = 36)
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
        scale_fill_manual("", values = div.colors)+theme_classic()+
        labs(color="")
ggsave(filename=paste0("Output_all/Fst/Overrepresented.AAs.ALL.pdf"), width = 7, height = 4.5)

comp1<-AAdiff[,c(1,2,3)]
comp1m<-melt(comp1)
ggplot(comp1m,aes(x=AA, y=value, fill=variable))+geom_bar(stat='identity',position='dodge', width = 0.8)+
        theme_light()+labs(x="Amino acid", y="Over/under-represented nonsyn AA (%)")+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = colors2[c(1,3)],labels=c("Costly in 1A", "Costly in 1B"))+
        geom_vline(xintercept = c(1:19)+0.5, color = "gray70", size=.1)+
        geom_hline(yintercept = 0, color = "gray50", size=.4)+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
              panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank(),
              panel.grid.major.y = element_line(size = 0.2, colour = "grey80", linetype=2))
ggsave(filename=paste0("Output_all/Fst2/Overrepresented.AAs.1A-1B.pdf"), width = 9, height = 4.5)

ggplot(comp1m,aes(x=AA, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_light()+labs(x="Amino acid", y="Over/under-represented nonsyn AA (%)")+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = colors2[c(1,3)],labels=c("Costly in 1A", "Costly in 1B"))+
        geom_vline(xintercept = c(1:19)+0.5, color = "gray70", size=.2)+
        geom_hline(yintercept = 0, color = "gray50", size=.4)+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
              panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank(),
              panel.grid.major.y = element_line(size = 0.2, colour = "grey80", linetype=2))
ggsave(filename=paste0("Output_all/Fst2/Overrepresented.AAs.1A-1B_2.pdf"), width = 7, height = 4.5)

comp2<-AAdiff[,c(1,4,5)]
comp2m<-melt(comp2)
ggplot(comp2m,aes(x=AA, y=value, fill=variable))+geom_bar(stat='identity',position='dodge', width = 0.8)+
        theme_light()+labs(x="Amino acid", y="Over/under-represented nonsyn AA (%)")+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = colors2[c(1,5)],labels=c("Costly in 1A", "Costly in 3A"))+
        geom_vline(xintercept = c(1:19)+0.5, color = "gray70", size=.2)+
        geom_hline(yintercept = 0, color = "gray50", size=.4)+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
              panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank(),
              panel.grid.major.y = element_line(size = 0.2, colour = "grey80", linetype=2))
ggsave(filename=paste0("Output_all/Fst2/Overrepresented.AAs.1A-3A.pdf"), width = 9, height = 4.5)
ggplot(comp2m,aes(x=AA, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_light()+labs(x="Amino acid", y="Over/under-represented nonsyn AA (%)")+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = colors2[c(1,5)],labels=c("Costly in 1A", "Costly in 3A"))+
        geom_vline(xintercept = c(1:19)+0.5, color = "gray70", size=.4)+
        geom_hline(yintercept = 0, color = "gray50", size=.4)+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
              panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank())
ggsave(filename=paste0("Output_all/Fst2/Overrepresented.AAs.1A-3A_2.pdf"), width = 7, height = 4.5)

comp3<-AAdiff[,c(1,6,7)]
comp3m<-melt(comp3)
ggplot(comp3m,aes(x=AA, y=value, fill=variable))+geom_bar(stat='identity',position='dodge', width = 0.8)+
        theme_light()+labs(x="Amino acid", y="Over/under-represented nonsyn AA (%)")+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = colors2[c(3,5)],labels=c("Costly in 1B", "Costly in 3A"))+
        geom_vline(xintercept = c(1:19)+0.5, color = "gray70", size=.2)+
        geom_hline(yintercept = 0, color = "gray50", size=.4)+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
              panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank(),
              panel.grid.major.y = element_line(size = 0.2, colour = "grey80", linetype=2))
ggsave(filename=paste0("Output_all/Fst2/Overrepresented.AAs.1B-3A.pdf"), width = 9, height = 4.5)
ggplot(comp3m,aes(x=AA, y=value, fill=variable))+geom_bar(stat='identity')+
        theme_light()+labs(x="Amino acid", y="Over/under-represented nonsyn AA (%)")+
        theme(plot.title = element_text(size=12))+
        scale_fill_manual("", values = colors2[c(1,5)],labels=c("Costly in 1A", "Costly in 3A"))+
        geom_vline(xintercept = c(1:19)+0.5, color = "gray70", size=.2)+
        geom_hline(yintercept = 0, color = "gray50", size=.4)+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
              panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank(),
              panel.grid.major.y = element_line(size = 0.2, colour = "grey80", linetype=2))
ggsave(filename=paste0("Output_all/Fst2/Overrepresented.AAs.1A-3A_2.pdf"), width = 7, height = 4.5)




## Genes
genefile<-list.files("Output_all/Fst2/", pattern="HighFst.byGene")
geneprop<-list()
for (i in 1:length(genefile)){
        dt<-read.csv(paste0("Output_all/Fst2/", genefile[i]), stringsAsFactors = F, row.names = 1)
        geneprop[[i]]<-dt[,c("Gene","Difference")]
        names(geneprop)[i]<-substr(genefile[i], start=16, stop = 33)
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
        scale_fill_manual("", values = colors2)+theme_classic()+
        labs(color="")+
        geom_hline(yintercept = 0, color = "gray50", size=.4)
ggsave(filename=paste0("Output_all/Fst2/Overrepresented.Gene.pdf"), width = 7, height = 4.5)

## codon

codonfile<-list.files("Output_all/Fst2/", pattern = glob2rx("*Codons*.csv*"))
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
ggsave(filename=paste0("Output_all/Fst/Overrepresented.NT.pdf"), width = 7, height = 4.5)



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
ggsave(filename=paste0("Output_all/Fst/Overrepresented.CpG.all.top5.pdf"), width = 4.5, height = 3)




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
        geom_hline(yintercept = 0,  color = "gray50", size=.5)+
        
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




