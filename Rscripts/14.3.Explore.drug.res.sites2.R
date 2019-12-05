library(plotrix)
library(ggplot2)
library(colorspace)
colors2<-qualitative_hcl(6, palette="Dark3")

#Ts summary with estiamted mutation frequency from the beta regression model
data<-read.csv("Output1A/GLM/Ts.with.EstimatedMF.csv",stringsAsFactors = F)
merge<-read.csv("Output_all/merged.metadata.csv",stringsAsFactors = F,  row.names = 1)
merge<-merge[,c("org.pos.1A","merged.pos")]
colnames(merge)[1]<-"pos"
#SC info
sc<-read.csv("Output1A/SelCoeff/SC.csv", row.names = 1, stringsAsFactors = F)
sc<-sc[,c("pos","EstSC")]
data<-merge(data, sc, by="pos")

#Look at the drug resistant sites
dr<-read.csv("Data/HCV_drugresistance_Geno.csv", stringsAsFactors = F)
dr<-dr[dr$genotype=="1A",]
dr<-merge(dr, merge, by="merged.pos")
drpos<-data$pos %in% dr$pos 
DRsites_all<-data[drpos,]
NonDRsites<-data[!drpos,]


DRsites$gene[DRsites$gene==8]<-"NS3"
DRsites$gene[DRsites$gene==11]<-"NS5A"
DRsites$gene[DRsites$gene==12]<-"NS5B"

#Average sc of dr sites:
mean(DRsites$EstSC)
#0.003064057
std.error(DRsites$EstSC)
#0.0002298942
#average of all EstSC # 0.002681807  

#Test the difference
wilcox.test(DRsites$EstSC, NonDRsites$EstSC[NonDRsites$gene==8|NonDRsites$gene==11|NonDRsites$gene==12], alternative = "greater", paired = FALSE)
# p-value = 0.02441
#SC at DR sites is significantly higher than nonDR sites

wilcox.test(DRsites$EstSC[DRsites$gene=="NS3"], NonDRsites$EstSC[NonDRsites$gene==8], alternative = "greater", paired = FALSE)
#p-value = 0.1736

wilcox.test(DRsites$EstSC[DRsites$gene=="NS5A"], NonDRsites$EstSC[NonDRsites$gene==11], alternative = "greater", paired = FALSE)
#p-value = 0.1555
wilcox.test(DRsites$EstSC[DRsites$gene=="NS5B"], NonDRsites$EstSC[NonDRsites$gene==12], alternative = "greater", paired = FALSE)
#p-value = 0.01898


# Average by genes
Genes<-c('NS3','NS5A','NS5B')
summary<-data.frame(Gene=Genes)
for (i in 1:3){
        dat<-DRsites[DRsites$gene==Genes[i],]
        summary$Total[i]<-nrow(dat)
        summary$SC[i]<-mean(dat$EstSC)
        summary$SC.se[i]<-std.error(dat$EstSC)
        summary$MF[i]<-mean(dat$mean)
        summary$MF.se[i]<-std.error(dat$mean)
        summary$meanMF_syn[i]<-mean(dat$mean[dat$Nonsyn==0])
        summary$menaMF_ns[i]<-mean(dat$mean[dat$Nonsyn==1])
        
   
        nt<-data.frame(table(dat$ref))
        summary$a[i]<-nt$Freq[1]
        summary$t[i]<-nt$Freq[4]
        summary$c[i]<-nt$Freq[2]
        summary$g[i]<-nt$Freq[3]
        
        summary$CpG[i]<-nrow(dat[dat$CpG==1,])
        summary$noCpG[i]<-nrow(dat[dat$CpG==0&(dat$ref=="a"|dat$ref=="t"),])

        
        #codon positions
        summary$codon1[i]<-nrow(dat[dat$codon==1,])
        summary$codon2[i]<-nrow(dat[dat$codon==2,])
        summary$codon3[i]<-nrow(dat[dat$codon==3,])
}

write.csv(summary,"Output1A/DrugRes/Summary_DRsites.csv")

#the highest SC sites
DRsites[which.max(DRsites$EstSC), ]
#   pos a t c g Syn Nonsyn Stop CpG bigAAChange        mean gene Core E1 HVR1 E2 NS1 NS2 NS3 NS4A NS4B NS5A NS5B
#   3813 1 0 0 0   0      1    0   1           0 0.002794247  NS3    0  0    0  0   0   0   1    0    0    0    0
#   codon ref WTAA EstimatedMF        diff top5       EstSC
#       1   a    I 0.004976229 0.002181981   NA 0.006334443

# nonsyn A to G transition mutations in NS3 at codon 1. 

ggplot(summary, aes(x=Gene, y=SC))+
        geom_point(size =2)+scale_color_manual(values=colors2[c(5)])+
        geom_errorbar(aes(ymin=SC-SC.se, ymax=SC+SC.se), width=.2, size=.2)+
        theme(axis.title.x=element_blank())+ylab("SC")+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray60", size=.4)+
        guides(
                shape = guide_legend(order = 2),
                color=guide_legend(order=1)
        )


