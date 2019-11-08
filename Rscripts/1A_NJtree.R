library(seqinr)
library(ape)
library(ggplot2)
library(ggtree)

sequences<-ape::read.dna("Output/alignments/HCV1aseq.fasta", format="fasta")
sequences

D<-dist.dna(sequences)
NJ<-nj(D)
ggtree(NJ)+geom_tiplab(aes(angle=angle),size=1, color="blue", vjust=-0.3)

ggtree(NJ,layout="daylight")+geom_tiplab(aes(angle=angle),size=2, color="blue", vjust=-0.3)


NAseq<-ape::read.dna("~/programs/Influenza/NA_alignment.fasta", format="fasta")
D<-dist.dna(NAseq)
NJ<-nj(D)
ggtree(NJ)


NAseqs<-read.fasta("~/programs/Influenza/NA_alignment.fasta")
###########
sequencenames<-row.names(summary(NAseqs))
Data<-data.frame(sequencenames)
Data$Has274Y<-0
Data$R222Q<-0
Data$V234M<-0
for (i in 1:length(NAseqs)){
        AA274<-translate(NAseqs[[i]])[275]
        if (AA274=="Y")
                Data$Has274Y[i]<-1
      R222Q<-translate(NAseqs[[i]])[222]
        if (R222Q=="Q")
                Data$R222Q[i]<-1
      V234M<-translate(NAseqs[[i]])[234]
      if (V234M=="M")
              Data$V234M[i]<-1
      
}


PhyMLAATree<-ape::read.tree("~/programs/Influenza/aaseqs_phy_phyml_tree.txt")

ggtree(PhyMLAATree)
which(PhyMLAATree$tip.label=="gb|CY00931810-1422|")

Data$ColorTree<-"grey"
Data$ColorTree[which(Data$Has274Y==1& Data$R222Q==1 & Data$V234M==1)]<-"red"
Data$ColorTree[which(Data$Has274Y==1& Data$R222Q==0 & Data$V234M==0)]<-"orange"
Data$ColorTree[which(Data$Has274Y==0 & Data$R222Q==1 & Data$V234M==1)]<-"blue"
Data$ColorTree[which(Data$Has274Y==0 & Data$R222Q==0 & Data$V234M==1)]<-"lightblue"
Data$ColorTree[which(Data$Has274Y==0 & Data$R222Q==1 & Data$V234M==0)]<-"lightblue"

ggtree(root(PhyMLAATree,52),lwd=.2)+
        geom_tippoint(color=Data$ColorTree,  size=1)

PhyMLAATree$tip.label[1:5]
Data$sequencenames[1:5]
