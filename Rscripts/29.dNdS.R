library(ggplot2)
library(reshape2)
library(colorspace)
library(ggthemes)
source("Rscripts/baseRscript.R")
colors2<-qualitative_hcl(6, palette="Dark3")
col2_light<-qualitative_hcl(6, palette="Set3")
div.colors<-c(colors2[1],col2_light[1],colors2[3],col2_light[3],colors2[5],col2_light[5])

## 1. Between Genotypes
# Results from SNPGenie Between Group Comparison
BetwResult1<-read.csv("Output_all/dNdS/between_group_codon_results.txt", sep = "\t", stringsAsFactors = F)
BetwResult1$product<-factor(BetwResult1$product, levels=c("Core","E1","HVR1","E2","NS1","NS2","NS3","NS4A","NS4B","NS5A","NS5B") )
colnames(BetwResult1)[2]<-c("gene")

Betwcodons<-list()
Betwcodons[[1]]<-BetwResult1[BetwResult1$group_1=="HCV1A_aligned.cds.fasta"&BetwResult1$group_2=="HCV3A_aligned.cds.fasta",]
names(Betwcodons)[1]<-"Betw1A3A"
Betwcodons[[2]]<-BetwResult1[BetwResult1$group_1=="HCV1B_aligned.cds.fasta"&BetwResult1$group_2=="HCV3A_aligned.cds.fasta",]
names(Betwcodons)[2]<-"Betw1B3A"
Betwcodons[[3]]<-BetwResult1[BetwResult1$group_1=="HCV1B_aligned.cds.fasta"&BetwResult1$group_2=="HCV1A_aligned.cds.fasta",]
names(Betwcodons)[3]<-"Betw1A1B"



#genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
vec.genes<-paste0(unique(BetwResult1$gene))

#count # of conserved sites
cod.sum<-list()
summary<-data.frame(gene=vec.genes)
for (i in 1:3){
        df<-Betwcodons[[i]]
        fname<-names(Betwcodons)[i]
        
        tb<-data.frame(gene=vec.genes)
        
        for (j in 1:length(vec.genes)){
                df1<-df[df$gene==vec.genes[j],]
                tb$conserved[j]<-nrow(df1[df1$variability=="conserved",])
                tb$polymorphic[j]<-nrow(df1[df1$variability=="polymorphic",])
        }
        tb$proportion<-tb$conserved/(tb$polymorphic+tb$conserved)
        plot(tb$proportion~tb$gene, main=fname)
        cod.sum[[i]]<-tb
        names(cod.sum)<-fname
        summary[,fname]<-tb$proportion
        
        
}

summary$gene<-factor(summary$gene, levels=c("Core","E1","HVR1","E2","NS1","NS2","NS3","NS4A","NS4B","NS5A","NS5B") )


write.csv(summary, "Output_all/dNdS/between_codontype_summary.csv")

sum2<-melt(summary)
ggplot(sum2, aes(x=gene, y=value*100, group=variable,fill=variable))+
        geom_bar(position=position_dodge(.9), stat="identity",width=0.8)+
        scale_fill_manual(values=paste0(colors2[c(1,3,5)]), labels=c("1A-3A","1B-3A","1A-1B"))+
        ylab("% conserved codons")+
        theme_bw()+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13),axis.title.x = element_text(size=13))+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray70", size=.5)+
        theme(panel.grid.major.x = element_blank())+
        theme(legend.title = element_blank())

ggsave("Output_all/dNds/ConservedCodonProportion.pdf", width = 7, height = 4.5)


#### plot dN/dS per gene for pairwise comparison between genotypes

BetwResult2<-read.csv("Output_all/dNdS/between_group_product_results.txt", sep = "\t", stringsAsFactors = F)
BetwResult2$product<-factor(BetwResult2$product, levels=c("Core","E1","HVR1","E2","NS1","NS2","NS3","NS4A","NS4B","NS5A","NS5B") )
colnames(BetwResult2)[2]<-c("gene")

BetwResult2$group[BetwResult2$group_1=="HCV1A_aligned.cds.fasta" & BetwResult2$group_2=="HCV3A_aligned.cds.fasta"]<-"1A-3A"
BetwResult2$group[BetwResult2$group_1=="HCV1B_aligned.cds.fasta"&BetwResult2$group_2=="HCV3A_aligned.cds.fasta"]<-"1B-3A"
BetwResult2$group[BetwResult2$group_1=="HCV1B_aligned.cds.fasta"&BetwResult2$group_2=="HCV1A_aligned.cds.fasta"]<-"1A-1B"

ratio<-BetwResult2[,c("gene","group","dN.dS.1")]

ggplot(ratio, aes(x=gene, y=dN.dS.1, group=group,fill=group))+
        geom_bar(position=position_dodge(.9), stat="identity",width=0.8)+
        scale_fill_manual(values=paste0(colors2[c(1,3,5)]), labels=c("1A-3A","1B-3A","1A-1B"))+
        ylab("dN/dS")+
        theme_bw()+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13),axis.title.x = element_text(size=13))+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray70", size=.5)+
        theme(panel.grid.major.x = element_blank())+
        theme(legend.title = element_blank())

ggsave("Output_all/dNds/dNdSratio.betweenGenotypes.SNPGenie.pdf", width = 7, height = 4.5)

## Plot dN & dS separately for each comparison
dn<-BetwResult2[,c("group","gene","dN")]
dn$grouptype<-paste0(dn$group,".dN")
colnames(dn)[3]<-c("divergence")
ds<-BetwResult2[,c("group","gene","dS")]
ds$grouptype<-paste0(ds$group,".dS")
colnames(ds)[3]<-c("divergence")

dNdS<-rbind(dn,ds)
ggplot(dNdS, aes(x=gene, y=divergence, group=grouptype,fill=grouptype))+
        geom_bar(position=position_dodge(.9), stat="identity",width=0.8,colour="gray60",size=.5)+
        scale_fill_manual(values=div.colors)+
        ylab("dN & dS")+
        theme_bw()+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13),axis.title.x = element_text(size=13))+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray70", size=.5)+
        theme(panel.grid.major.x = element_blank())+
        theme(axis.title.x=element_blank())+
        theme(legend.title = element_blank())
        

ggsave("Output_all/dNds/dNdSbetween_genotypes.pdf", width = 9, height = 4.5)





## 2. Within genotype Pi, pN, pS
# 2.1 Average pi (D) per gene by genotype, calculated from the consensus

Genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
Genes$Gene[6]<-"NS1"
gene.vector<-c()
for (i in 1:12){
        gene.vector<-c(gene.vector, rep(Genes$Gene[i],times=Genes$start[i+1]-Genes$start[i]))
}

geno<-c("1A","1B","3A")

for (i in 1:3){
        df<-read.csv(paste0("Output_all/dNdS/",geno[i],"_within_group_site_results.txt"), sep = "\t", stringsAsFactors = F)
        dname<-paste0("r1",geno[i])
        df<-df[1:9410,-1]
        df$Gene<-gene.vector[1:nrow(df)]
        assign(dname, df)
}

pi<-data.frame(Gene=Genes$Gene[1:12])
pi2<-data.frame(Gene=Genes$Gene[1:12])
pi.sum<-data.frame()
Pilist<-list()
for (i in 1:3){
        df<-get(paste0("r1",geno[i]))
        pi2<-data.frame(Gene=Genes$Gene[1:12])
        for (k in 1:nrow(pi)){
                pi[k,paste0("mean.",geno[i])]<-mean(df$pi[df$Gene==pi$Gene[k]], na.rm=T)
                pi[k,paste0("se.",geno[i])]<-std.error(df$pi[df$Gene==pi$Gene[k]], na.rm=T)
                pi2[k,"Mean"]<-mean(df$pi[df$Gene==pi$Gene[k]], na.rm=T)
                pi2[k,"SE"]<-std.error(df$pi[df$Gene==pi$Gene[k]], na.rm=T)
        }
        pi2$Genotype<-geno[i]
        Pilist[[i]]<-pi2
        names(Pilist)[i]<-geno[i]
        pi.sum<-rbind(pi.sum, pi2)
        
}

write.csv(pi, "Output_all/dNdS/Pi_by_geno.by.genotypes.csv")
#write.csv(pi2, "Output_all/dNdS/Pi_by_geno.by.genotypes2.csv")

pi$Gene<-factor(pi$Gene, levels=c("5' UTR","Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))
pi.sum$Gene<-factor(pi$Gene, levels=c("5' UTR","Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))



ggplot(pi.sum,aes(x=Gene,y=Mean,group=Genotype, color=Genotype))+
        geom_point(position=position_dodge(width=0.3),size =1.5)+scale_color_manual(values=colors2[c(1,3,5)])+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, position=position_dodge(width=0.3))+
        theme(axis.title.x=element_blank())+ylab(expression(paste("Average ",pi)))+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:11)+0.5,  
                   color = "gray60", size=.5)

ggsave(filename="Output_all/dNdS/Mean.D(pi).by.gene.by.genotype.pdf",width = 8.5, height = 5)

# 2.2 plot within genotype dN/dS (pN/pS) from SNPGeneie

dnds<-data.frame()
for (i in 1:3){
        df<-read.csv(paste0("Output_all/dNdS/",geno[i],"_within_group_product_results.txt"), sep = "\t", stringsAsFactors = F)
        dname<-paste0("r2.",geno[i])
        df$product<-factor(df$product, levels=c("Core","E1","HVR1","E2","NS1","NS2","NS3","NS4A","NS4B","NS5A","NS5B") )
        df[,1]<-geno[i]
        colnames(df)[1:2]<-c("Genotype","gene")
        dnds<-rbind(dnds,df)
        assign(dname, df)
}


ggplot(dnds, aes(x=gene, y=dN_over_dS, group=Genotype,fill=Genotype))+
        geom_bar(position=position_dodge(.9), stat="identity",width=0.8)+
        scale_fill_manual(values=colors2[c(1,3,5)])+
        ylab("dN/dS")+
        theme_bw()+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13),axis.title.x = element_text(size=13))+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray70", size=.5)+
        theme(panel.grid.major.x = element_blank())

ggsave("Output_all/dNds/dNdSratio.withinGenotype.SNPGenie.pdf", width = 7, height = 4.5)

#2.3 Plot dN & dS within genotype

#formatting the data
dN<-dnds[,c("Genotype","gene","dN","SE_dN")]
dN$group<-paste0(dnds$Genotype,".dN")
colnames(dN)[3:4]<-c("divergence", "SE")
dS<-dnds[,c("Genotype","gene","dS","SE_dS")]
dS$group<-paste0(dnds$Genotype,".dS")
colnames(dS)[3:4]<-c("divergence", "SE")

colour=c(colors2[1],colors2[1] ,colors2[3],colors2[3],colors2[5],colors2[5])

Div<-rbind(dN,dS)
ggplot(Div, aes(x=gene, y=divergence, group=group,fill=group))+
        geom_bar(position=position_dodge(.9), stat="identity",width=0.8,colour="gray60",size=.5)+
        scale_fill_manual(values=div.colors)+
        ylab("dN & dS -within genotype")+
        geom_errorbar(position=position_dodge(.9),aes(ymin=divergence-SE, ymax=divergence+SE), colour="gray60",width=.2)+
        theme_bw()+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13),axis.title.x = element_text(size=13))+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray70", size=.5)+
        theme(panel.grid.major.x = element_blank())

ggsave("Output_all/dNds/dNdSwithin_genotype.pdf", width = 9, height = 4.5)


#########
hyphy1a<-read.csv(paste0("Data/datamonkey-HYPHY_1A.csv"), stringsAsFactors = F)

genes2<-read.csv("Data/HCV_annotations2.csv", stringsAsFactors = F)
genes2$Gene[genes2$Gene=="NS1(P7)"]<-"NS1"
genes2<-genes2[2:13,]
codon.vector2<-c()


for (i in 1:(nrow(genes2)-1)){
        codon.vector2<-c(codon.vector2, rep(genes2$Gene[i],times=(genes2$start[i+1]-genes2$start[i])/3))
}

Hyphy<-list()
hyphy1a$gene<-codon.vector2[1:nrow(hyphy1a)]
Hyphy[[1]]<-hyphy1a
names(Hyphy)[1]<-geno[i]
        




######
#3. Compare HyPhy vs. SNPGenie
# dNdS estimated from HyPhy

geno<-c("1A","1B","3A")

genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
genes<-genes[2:13,]
codon.vector<-c()
for (i in 1:(nrow(genes)-1)){
        codon.vector<-c(codon.vector, rep(genes$Gene[i],times=(genes$start[i+1]-genes$start[i])/3))
}

#Hyphy<-list()
for (i in 2:3){
        hyphy<-read.csv(paste0("Data/datamonkey-HYPHY_", geno[i],".csv"), stringsAsFactors = F)
        hyphy$gene<-codon.vector[1:nrow(hyphy)]
        Hyphy[[i]]<-hyphy
        names(Hyphy)[i]<-geno[i]
        
}

D<-data.frame(gene=genes$Gene[1:11])
D1<-data.frame(gene=genes$Gene[1:11])
Dhy<-data.frame()
for (i in 1:3){
        DF<-Hyphy[[i]]
        for(k in 1:11){
                dt<-DF[DF$gene==genes$Gene[k],]
                D[k,geno[i]]<-mean(dt$X.beta.)/mean(dt$X.alpha.)
                D1$dN_over_dS[k]<-mean(dt$X.beta.)/mean(dt$X.alpha.)
        
        }
        D1$Genotype<-geno[i]
        D1$group<-paste0(geno[i]," HyPhy")
        Dhy<-rbind(Dhy,D1)
}

Dhy$program<-"HyPhy"

write.csv(D, "Output_all/dNdS/dnds.ratio.hyphy.byGene.csv")

#plot hyphy and SNPGenie together
Dsnpg<-data.frame()
for ( i in 1:3){
        df<-dnds[dnds$Genotype==geno[i],]
        df<-df[,c(2,12,1)]
        df$group<-paste0(geno[i]," SNPGenie")
        df$program<-"SNPGenie"
        Dsnpg<-rbind(Dsnpg,df)
}


D2<-rbind(Dsnpg,Dhy)



ggplot(D2, aes(x=gene, y=dN_over_dS, group=group,fill=group))+
        geom_bar(position=position_dodge(.9), stat="identity",width=0.8,colour="gray60",size=.5)+
        scale_fill_manual(values=div.colors)+
        ylab("dN/dS")+
        theme_bw()+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13),axis.title.x = element_text(size=13))+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray70", size=.5)+
        theme(panel.grid.major.x = element_blank())+
        theme(axis.title.x=element_blank())+
        theme(legend.title = element_blank())
ggsave("Output_all/dNds/Compare_dNdSRatio.byPrograms.pdf", width = 9, height = 4.5)






ggplot(D2, aes(x=gene, y=dN_over_dS, group=program,fill=program))+
        geom_bar(position=position_dodge(.9), stat="identity",width=0.8)+
        scale_fill_manual(values=colors2[c(1,4)])+
        ylab("dN/dS")+
        theme_bw()+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13),axis.title.x = element_text(size=13))+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray70", size=.5)+
        theme(panel.grid.major.x = element_blank())

ggsave("Output_all/dNds/Compare_dNdSRatio.within_geno.byPrograms.pdf", width = 9, height = 4.5)




###############
result1A<-read.csv("~/programs/SNPGenie-master/HCV/Within-group/1A_sw_10codons_results.txt", sep = "\t", stringsAsFactors = F)
result1A$product<-factor(result1A$product, levels=c("Core","E1","HVR1","E2","NS1","NS2","NS3","NS4A","NS4B","NS5A","NS5B") )
result1A<-result1A[order(result1A$product),]
result1A$pos<-1:nrow(result1A)
which(result1A$window==100)
result1A<-result1A[1:2430,]
result1A$dN_over_dS<-as.numeric(result1A$dN_over_dS)
plot(result1A$dN_over_dS~result1A$pos, pch=".")

ggplot(result1A, aes(x=pos, y=dN_over_dS))+geom_bar(stat="identity")
ggplot(result1A, aes(x=pos, y=dN_over_dS))+geom_line()
 