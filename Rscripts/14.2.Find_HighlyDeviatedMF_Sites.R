library(ggplot2)
library(miscTools)
library(tidyverse)
library(janitor)
library(reshape2)
library(dplyr)


cols2<-c("#66CCEE","#EE6677" ,"#228833")
cols<-c("#66CCEE", "#228833" ,"#CCBB44", "#EE6677" ,"#AA3377", "#4477AA", "#BBBBBB")
cols3<-c("#009988CC" ,"#66CCEECC", "#EE6677CC", "#4477AACC")

####
data<-read.csv("Output1A/GLM/GlmdataFull.Ts.Q35.csv", stringsAsFactors = F, row.names = 1)

#If run previously, read the EstimatedMF
#data<-read.csv("Output1A/GLM/Ts.with.EstimatedMF.csv",stringsAsFactors = F)


##add the codon position
vec<-data.frame('pos'= 342:8609)
dfm<-merge(data, vec, by="pos", all.y=TRUE)
dfm$codon<-rep(1:3, times=nrow(dfm)/3) 
#remove the sites with NA
data<-dfm[!is.na(dfm$mean),]

#Attach the reference sequence/AA info
Transitions<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F)
data<-merge(data, Transitions[,c("pos","ref","WTAA")], by ="pos")

#remove the stop mutations
data<-data[data$Stop==0,]

s<-length(list.files("Output1A/Overview3/",pattern="overview3.csv"))

modcoef<-read.csv("Output1A/GLM/BetaReg_BestModel.Q35.csv",stringsAsFactors = F)
rownames(modcoef)<-modcoef$X
modcoef<-modcoef[,-1]
coef.vals <- modcoef[,1]


#from GLMPlotFunctions.R
EstimateMF<- function(NsOrNo = 0, CpGorNo = 0, bigAAChangeOrNo = 0,CoreorNo=0,E1orNo=0, HVRorNo=0,E2orNo=0,NS1orNo=0,NS2orNo=0,NS4AorNo=0,NS5BorNo=0){
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


data$EstimatedMF<-0
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
        
        setUpDat <- EstimateMF(NsOrNo, CpGorNo, bigAAChangeOrNo, CoreorNo,E1orNo, HVRorNo,E2orNo,NS1orNo,NS2orNo,NS4AorNo,NS5BorNo)[,rownames(modcoef)]
        
        if (data$a[i]==1) rown =1
        if (data$t[i]==1) rown =2
        if (data$c[i]==1) rown =3
        if (data$g[i]==1) rown =4
        
        data$EstimatedMF[i] <- exp(setUpDat[rown,] %*% coef.vals)
}

# NS3, NS4B and NS5A are not in the Best Fit model # 
#visualize the results
plot(data$mean[500:2000], pch=16, col="#66CCEE", cex=.5, xaxt='n', ylab='Mutation frequency',xlab='Genome position')
points(data$EstimatedMF[500:2000],pch=16, col="blue",cex=0.5)
xticks<-seq(0,1500, by=500)
axis(side=1, at=xticks,labels=seq(500,2000, by = 500))
legend('topright', legend=c('Ovserved', 'Estimated'), col=c("#66CCEE",'blue'), pch= 16,cex=.7)

pdf("Output1A/GLM/Expected.vs.observed.mf.pdf", width=10, height = 6)
plot(data$mean, pch=".", col="#66CCEE", cex=1, ylab='Mutation frequency',xlab='Genome position')
points(data$EstimatedMF,pch=".", col="blue",cex=1)
legend('topright', legend=c('Ovserved', 'Estimated'), col=c("#66CCEE",'blue'), pch= 16,cex=.7)
dev.off()


# find the differences betweeen [estimated - observed]
data$diff<- data$EstimatedMF-data$mean

pdf("Output1A/GLM/Expected.vs.observed.mf-Difference.pdf", width=10, height = 6)
plot(data$diff,pch=".", col="#66CCEE", cex=1, ylab="Expected - Observed mut freq",xlab="Genome position")
dev.off()

# create the summary data based on mutation type
MFdiff<-data.frame("Mutation Type"= c("Overall", "Syn","Nonsyn", "Syn-CpG", "Syn-nonCpG","NS-CpG","NS-nonCPG"))
MFdiff$Mean.diff<-c(mean(data$diff),mean(data$diff[data$Syn==1]),mean(data$diff[data$Nonsyn==1]), 
                    mean(data$diff[data$Syn==1&data$CpG==1]),mean(data$diff[data$Syn==1&data$CpG==0&(data$ref=="a"|data$ref=="t")]),
                    mean(data$diff[data$Nonsyn==1&data$CpG==1]),mean(data$diff[data$Nonsyn==1&data$CpG==0&(data$ref=="a"|data$ref=="t")]))
MFdiff$Mutation.Type<-factor(MFdiff$Mutation.Type, c("Overall", "Syn","Nonsyn","Syn-CpG", "Syn-nonCpG","NS-CpG","NS-nonCPG"))
write.csv(MFdiff, "Output1A/GLM/Expected.vs.observed.mf-Difference2.csv")

label_scientific <- function(l) {
        # turn in to character string in scientific notation
        l <- format(l, scientific = TRUE)
        # quote the part before the exponent to keep all the digits
        l <- gsub("^(.*)e", "'\\1'e", l)
        # turn the 'e+' into plotmath format
        l <- gsub("e", "%*%10^", l)
        # return this as an expression
        parse(text=l)
}
MF<-MFdiff[1:3,]
ggplot(MF,aes(x=Mutation.Type,y=Mean.diff))+
        geom_bar(stat="identity", color="#EE6677",fill="#EE667766",width=0.8)+
        theme_bw()+labs(x="Mutation type",y="Average difference between observed & expected MF")+
        theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
        geom_hline(yintercept=0, size=0.4, color="gray10")+
        scale_y_continuous(labels=label_scientific) 
ggsave("Output1A/GLM/Mean.diff.Exp.vs.obs.pdf", width=4, height=5)

ggplot(MFdiff,aes(x=Mutation.Type,y=Mean.diff))+
        geom_bar(stat="identity", color="#EE6677",fill="#EE667766",width=0.8)+
        theme_bw()+labs(x="Mutation type",y="Ave. [Expected - Observed] MF")+
        theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
        geom_hline(yintercept=0, size=0.4, color="gray10")+
        scale_y_continuous(labels=label_scientific) 
ggsave("Output1A/GLM/Exp-ObsMF_cpg.pdf", width=7, height=5)


#Find the most conserved 5% sites (lower than expected mf sites)
l<-nrow(data)
data2<-head(data[order(data$diff, decreasing= T),], n=round(l*0.05))
range(data2$diff)
#0.00337167 0.01081907

data$top5<-sapply(data$diff,function(x){ if (x>min(data2$diff)) x=x
else x<-NA})
plot(data2$pos, data2$diff,pch=16, col="#66CCEE", cex=.5)

#write.csv(data, "Output1A/GLM/Ts.with.EstimatedMF.csv",row.names = F)



#####  Look at the results by gene
genes<-read.csv("Data/HCV_annotations2.csv", stringsAsFactors = F)

ByGenes<-list()
for (i in 2:12){
        df<-data[data$gene==i,]
        ByGenes[[(i-1)]]<-df
        names(ByGenes)[(i-1)]<-paste0(genes$Gene[i])
}

pdf(file="Output1A/GLM/Compare_Estimated.vsObservedMF_top5.pdf", width=8,height=11.5)
par(mfrow = c(6,2),mar=c(2,2,1,1))
for (i in 1:length(ByGenes)){
        df<-ByGenes[[i]]
        if (length(which(!is.na(df$top5)))==0) next
        else {
                plot(df$pos,df$top5, pch='', col=cols2[1], cex=.8, main = names(ByGenes)[i], xlab='Genome position', ylab= '(Observed - Estimated) mutation Frequency')  
                points(df$pos[df$top5>0],df$top5[df$top5>0], pch=16, col=cols2[1], cex=.8)  
                abline(h=0, col="gray")
        }
}
dev.off()


Counts<-data.frame()
df<-data.frame(gene=names(ByGenes))
for (i in 2:12){
        df$counts[i-1]<-length(which(!is.na(data$top5[data$gene==i])))
}

pdf("Output1A/GLM/ConservedSiteCOunts_byGene.pdf",width = 9, height = 6)
barplot(df$counts, names.arg=df$gene)
dev.off()



#################
#which nucleotides are over/under- represented in the conserved sites.
Sum<-list()
    Sum[[1]]<-table(data2$ref)
    Sum[[2]]<-table(data$ref)
    Sum[[3]]<-Sum[[1]]/sum(Sum[[1]])
    Sum[[4]]<-Sum[[2]]/sum(Sum[[2]])
    Sum[[5]]<-Sum[[3]]-Sum[[4]]
NT<-as.data.frame(Sum[[5]]/Sum[[4]]*100)

NT2<-data.frame(nt=c("A","C","G","T"), Type=rep("all",times=4), Percent=as.vector(t(Sum[[4]])))
NT2$nt<-factor(NT2$nt,c("A","T","C","G"))

NT2.2<-data.frame(nt=c("A","C","G","T"), Type=rep("conserved",times=4), Percent=as.vector(t(Sum[[3]])))
NT2.2$nt<-factor(NT2.2$nt,c("A","T","C","G"))

NT2<-rbind(NT2,NT2.2)
#NT2$pos[1:4]<-cumsum(NT2$Percent[1:4]) - NT2$Percent[1:4]/2
#NT2$pos[5:8]<-cumsum(NT2$Percent[5:8]) - NT2$Percent[5:8]/2

vec<-rev(c(1,4,2,3))
NT2$pos2[1:4]<-cumsum(NT2$Percent[vec])-NT2$Percent[vec]/2
vec<-rev(c(1,4,2,3)+4)
NT2$pos2[5:8]<-cumsum(NT2$Percent[vec])-NT2$Percent[vec]/2


ggplot(NT,aes(x=Var1,y=Freq))+
        geom_bar(stat="identity", color="#EE6677",fill="#EE667766",width=0.8)+
        theme_bw()+labs(x="")+
        theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
        geom_hline(yintercept=0, size=0.4, color="gray10")+
        scale_x_discrete(labels=c("A","C","G","T"))+
        theme(axis.text.x = element_text(size=14),
              axis.text.y = element_text(size=14),axis.title.y = element_text(size = 14))+
        ylab(expression(paste("Over/under-represented NT (%) \n at conserved sites")))+
        theme(plot.margin=unit(c(5,5,5,20),"points"))
ggsave("Output1A/GLM/Overrepresented.Nucleotides.in.Top5.pdf", width=5, height=4)

ggplot(NT,aes(x=Var1,y=Freq))+
        geom_bar(stat="identity", color="#EE6677",fill="#EE667766",width=0.8)+
        theme_bw()+labs(x="")+
        theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
        geom_hline(yintercept=0, size=0.4, color="gray10")+
        scale_x_discrete(labels=c("A","C","G","T"))+
        theme(axis.text.x = element_text(size=14),
              axis.text.y = element_text(size=14),axis.title.y = element_text(size = 14))+
        ylab(expression(paste("Over/under-represented NT (%) \n at conserved sites")))+
        theme(plot.margin=unit(c(5,5,5,20),"points"))
ggsave("Output1A/GLM/Overrepresented.Nucleotides.in.Top5.pdf", width=5, height=4)

### plot side by side

ggplot(NT2,aes(nt,Percent))+
        geom_bar(aes(fill=Type), position="dodge", stat="identity")+
        theme(axis.title.x=element_blank())+ylab("Proportion of each nucleotide")+xlab("")+
        scale_x_discrete(labels=c("A","C","G","T"))+
        scale_fill_manual(values=cols3, labels=c("Entire genome","Conserved sites"))+
        theme_classic()+
        guides(fill=guide_legend(title=NULL))
ggsave(filename="Output1A/GLM/NT_percent_comparison.pdf",width = 5, height = 3.5)

##stacked bar plot
ggplot(NT2, aes(fill=nt, y=Percent, x=Type))+
        geom_bar( stat="identity", position="fill",width = 0.4)+
        ylab("Proportion of each nucleotide")+xlab("")+
        scale_x_discrete(labels=c("Entire genome","Conserved sites"))+
        scale_fill_manual(values=cols3, labels=c("A","T","C","G"))+
        theme_light()+
        theme(legend.position = "none")+
        geom_text(aes(label=c("G","C","T","A","G","C","T","A"),y=pos2),size=3)
ggsave(filename="Output1A/GLM/NT_percent_comparison2.pdf",width = 3.5, height = 5)


table(data2$CpG, data2$ref)
#            a   c   g   t
#0(nonCPG)  44 145  57  72
#1(CpG)     26   0   0  43

#proportion of over/under represented AA by gene:
#Adjust the length of the last gene
genes2<-genes[-13,]
genes2$end[12]<-data$pos[nrow(data)]
PercentDeviates<-data.frame("Gene"= paste0(genes$Gene[2:12]), "counts"=rep('',times=11), stringsAsFactors = FALSE)
codonsummary<-data.frame("Gene"= paste0(genes$Gene[2:12]),"first"=rep(0,times=11),stringsAsFactors = FALSE)
AAsummary<-list()
for (i in 2:12){
        vec<-data.frame('pos'=c(genes$start[i]:((genes$end[i]-genes$start[i])+genes$start[i])))
        df<-data2[data2$gene==i,]
        PercentDeviates$counts[i-1]<-nrow(vec)
        PercentDeviates$length[i-1]<-nrow(df)
        PercentDeviates$Proportion[i-1] <- ((nrow(df)/nrow(data2))/(nrow(vec)/nrow(data))-1 )
        codonsummary$first[i-1] <-sum(df$codon[!is.na(df$mean)& df$codon==1])
        codonsummary$second[i-1]<-sum(df$codon[!is.na(df$mean)& df$codon==2])
        codonsummary$third [i-1]<-sum(df$codon[!is.na(df$mean)& df$codon==3])
        aa<-as.data.frame(table(df$WTAA))
        #aa$Var1<-as.character(aa$Var1)
        AAsummary[[i-1]]<-aa
        names(AAsummary)[i-1]<-genes$Gene[i]
        
}

write.csv(PercentDeviates, "Output1A/GLM/Prop.of.Overrepresented.sites.by.gene.csv")

#NS2 doesn't have any sites so format the data frame
AAsummary[[6]]<-data.frame(Var1='A',Freq=0)
for (i in 1:length(AAsummary)) colnames(AAsummary[[i]])<-c("Gene",paste0(names(AAsummary[i])))
AA<-AAsummary%>% purrr::reduce(full_join, by="Gene")
AA[is.na(AA)]<-0
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
ggsave("Output1A/GLM/Overrepresented.sites.perGene.Top5.pdf", width=8, height=6)


# Over or under represented amino acids
AA$Total<-rowSums(AA[2:12])
colnames(AA)[1]<-"AA"
AA2<-AA[,c("AA","Total","Prop.conserved")]
AA$Prop.conserved<-AA$Total/sum(AA$Total)

AAall<-as.data.frame(table(data$WTAA))
colnames(AAall)<-c("AA","Counts")
AAall$Prop<-AAall$Counts/sum(AAall$Counts)
AAall<-merge(AAall, AA2, by="AA", all.x = T)
AAall[is.na(AAall)]<-0
AAall$proportion<-(AAall$Prop.conserved-AAall$Prop)/AAall$Prop*100

cols4<-substr(cols3, start=1,stop = 7)
ggplot(AAall,aes(x=AA,y=proportion))+
        geom_bar(stat="identity",width=0.8,color=cols4[4],fill=paste0(cols4[4],"CC"))+
        theme_bw()+labs(x="")+
        theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
        geom_hline(yintercept=0, size=0.4, color="gray10")+
        theme(axis.text.x = element_text(size=12),
              axis.text.y = element_text(size=12),axis.title.y = element_text(size = 14))+
        ylab(expression(paste("Over-/Under-represented AA \n at conserved sites (%)")))+
        theme(plot.margin=unit(c(5,5,5,20),"points"))
ggsave("Output1A/GLM/Overrepresented.AA.perGene.Top5_blue.pdf", width=10, height=6)


#Look at the drug resistant sites
dr2<-read.csv("Data/HCV_drugresistance_NTpos.csv")
drpos<-data$pos %in% dr2$NTpos 
DRsites_all<-data[drpos,]
#Ave mut freq of known drug resistant sites
mean(DRsites_all$mean)# 0.004466228

mean(DRsites_all$mean[DRsites_all$Nonsyn==1])  # 0.003669182
mean(DRsites_all$mean[DRsites_all$Nonsyn==0])  # 0.007426685

table(DRsites_all$Nonsyn) # Majority of drug resistant sites are Nonsyn sites
#0  1 
#7 26

######
# Look at the drug resistant sites in the conserved 5% sites: 
DRsites_h<-data2[(data2$pos%in% dr2$NTpos),]
# only 1 site 
# pos a t c g Syn Nonsyn Stop CpG bigAAChange        mean gene Core E1 HVR1 E2 NS1 NS2 NS3 NS4A NS4B NS5A NS5B codon ref WTAA EstimatedMF        diff
#6347 1 0 0 0   1      0    0   0           0 0.007877798   11    0  0    0  0   0   0   0    0    0    1    0     3   a    Q  0.01142332 0.003545527


