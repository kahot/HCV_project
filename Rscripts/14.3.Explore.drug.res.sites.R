## Old list of RAVs: Use New one. 


#Ts summary with estiamted mutation frequency from the beta regression model
data<-read.csv("Output1A/GLM/Ts.with.EstimatedMF.csv",stringsAsFactors = F)


#Look at the drug resistant sites
dr2<-read.csv("Data/HCV_drugresistance_NTpos.csv")
drpos<-data$pos %in% dr2$NTpos 
DRsites_all<-data[drpos,]


## Look at the 3 genes with known drug resistance
genes<-read.csv("Data/HCV_annotations2.csv", stringsAsFactors = F)

ByGenes<-list()
for (i in 2:12){
        df<-data[data$gene==i,]
        ByGenes[[(i-1)]]<-df
        names(ByGenes)[(i-1)]<-paste0(genes$Gene[i])
}

ByGene3<-list()
k=1
for (i in c(7,10,11)){
        ByGene3[[k]]<-ByGenes[[i]]
        names(ByGene3)[k]<-names(ByGenes)[i]
        k=k+1
}

Summary1<-data.frame("Gene"=c('NS3','NS5A','NS5B'))
for (i in 1:3){
        dat<-ByGene3[[i]]
        dat$ref<-as.character(dat$ref)
        #summary of highly deviated sites form expected mut freq
        dat1<-dat[!is.na(dat$top5),]
        
        Summary1$No_sites[i]<-nrow(dat)
        Summary1$No_conserved[i]<-nrow(dat1)
        Summary1$mean_diff[i]<-mean(dat$diff)
        Summary1$meanMF_syn[i]<-mean(dat$mean[dat$Nonsyn==0])
        Summary1$menaMF_ns[i]<-mean(dat$mean[dat$Nonsyn==1])
        
        nt<-data.frame(table(dat$ref))
        Summary1$a[i]<-nt$Freq[1]
        Summary1$t[i]<-nt$Freq[4]
        Summary1$c[i]<-nt$Freq[2]
        Summary1$g[i]<-nt$Freq[3]
        
        Summary1$CpG[i]<-nrow(dat[dat$CpG==1,])
        Summary1$noCpG[i]<-nrow(dat[dat$CpG==0&(dat$ref=="a"|dat$ref=="t"),])
        
        
        Summary1$Max.pos[i]<- dat1$pos[which.max(dat1$diff)]
        Summary1$Max_base[i]<-dat1$ref[which.max(dat1$diff)]
        Summary1$Max_codon[i]<-dat1$codon[which.max(dat1$diff)]
        
        #codon positions
        Summary1$No_codon1[i]<-nrow(dat1[dat1$codon==1,])
        Summary1$No_codon2[i]<-nrow(dat1[dat1$codon==2,])
        Summary1$No_codon3[i]<-nrow(dat1[dat1$codon==3,])
}

write.csv(Summary1,"Output1A/GLM/Summary_DeviatedMF_sites_3DRgenes.csv")


# most conserved drug resistant sites? all 't'
data[data$pos%in%Summary1$Max.pos,]
#      pos a t c g Syn Nonsyn Stop CpG bigAAChange        mean gene Core E1 HVR1 E2 NS1 NS2
#4446 4919 0 1 0 0   1      0    0   0           0 0.003309432    8    0  0    0  0   0   0
#6217 6759 0 1 0 0   1      0    0   0           0 0.005388616   11    0  0    0  0   0   0
#7695 8303 0 1 0 0   1      0    0   0           0 0.003255265   12    0  0    0  0   0   0
#     NS3 NS4A NS4B NS5A NS5B codon ref WTAA EstimatedMF        diff        top5
#4446   1    0    0    0    0     3   t    A  0.01282502 0.009515593 0.009515593
#6217   0    0    0    1    0     1   t    L  0.01282502 0.007436408 0.007436408
#7695   0    0    0    0    1     3   t    R  0.01407433 0.010819067 0.010819067




NS3s_low<-NS3s[NS3s$diff<0,] #65 sites
t1<-data.frame(table(NS3s_low$codon))   # Most are 1st and 2nd positions in the codon
#1  2  3 
#24 40  1 
table(NS3s_low$Nonsyn)  # all are nonsyn mutations
#1 
#65

table(NS3s_low$CpG)  # most are non-CpG sites
# 0  1 
#57  8 
table(NS3s_low$bigAAChange)  # most are Big AA cahnge sites
# 0  1 
#12 53 



###
#2. NS5A

vec2<-data.frame('pos'=c(genes$start[11]:((genes$end[11]-genes$start[11])+genes$start[11])))

NS5A_com<-merge(NS5A, vec2, by="pos", all.y=TRUE)
NS5A_com$codon <- rep(1:3, times=nrow(NS5A_com)/3)

mean(NS5A_com$diff,na.rm=T)  #-0.003766768
NS5A_com[which.max(NS5A_com$diff),]
#     pos a t c g Syn Nonsyn Stop CpG bigAAChange       mean gene hvr EstimatedMF       diff   over0.01 ref codon
#927 7184 0 0 1 0   1      0    0   0           0 0.02834871   11   0 0.007553138 0.02079558 0.02079558   c     3
NS5A_com[which.min(NS5A_com$diff),]
#    pos a t c g Syn Nonsyn Stop CpG bigAAChange        mean gene hvr EstimatedMF        diff    over0.01 ref codon
#35 6292 0 1 0 0   0      1    0   0           1 0.001878365   11   0  0.01412755 -0.01224919 -0.01224919   t     2

####################

data

