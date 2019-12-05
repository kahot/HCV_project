library(plotrix)
library(reshape)
library(tidyverse)
library(zoo)
library(purrr)
library(ape)
library(seqinr)
library(DescTools)
library(ggplot2)
library(colorspace)
library(emojifont)
source("Rscripts/baseRscript.R")
colors2<-qualitative_hcl(6, palette="Dark3")


align<-read.fasta("Output1A/HCV1A_Consenssu_Alignment.fasta",as.string=TRUE)
consensus<-read.fasta("Data/HCVref.fasta",as.string=TRUE)
consen<-paste(consensus)
consen<-unlist(strsplit(consen, "") )
consen<-consen[289:8613]

al.ta<-data.frame(pos=289:8613)
for (i in 1: length(align)){
        seq<-paste(align[i])
        sname<-substr(names(align)[i], start=1,stop =6)
        seq<-unlist(strsplit(seq, "") )
        al.ta[,sname]<-seq
}


dt<-al.ta[,2:196]
dt2<-dt
dt2[,]<-0
for (i in 1:nrow(dt)){
        for (k in 1: ncol(dt)){
                if (dt[i,k]==consen[i]) dt2[i,k]<-1
        }
}

dt2$Sum<-rowSums(dt2[,1:195])      
write.csv(dt2, "Output1A/Identical_sites_matrix.csv")

###
dt2<-read.csv("Output1A/Identical_sites_matrix.csv", row.names = 1, stringsAsFactors = F)

same.sites<-which(dt2$Sum==195)


### Use filtered data first  ###

# separate identical sites and non-identical sites
Ts<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",row.names = 1,stringsAsFactors = F)
Ts<-Ts[Ts$pos>=289& Ts$pos<=8613,]
T_same<-Ts[same.sites,]
T_diff<-Ts[-same.sites,]

mean(T_same$mean, na.rm=T)  #0.004690898
mean(T_diff$mean, na.rm=T) #0.004929429

wilcox.test(T_same$mean, T_diff$mean,  alternative = "less", paired = FALSE)
#W = 7665700, p-value = 0.0003269

nt<-data.frame(table(T_same$ref))
nt2<-data.frame(table(T_diff$ref))

table(T_same$makesCpG)
table(T_diff$makesCpG)

nt$Diff<-nt2$Freq
colnames(nt)[1:2]<-c("nt", "Same")
GTest(nt[,2:3])
#G = 1.8794, X-squared df = 3, p-value = 0.5978


table(T_same$Type)
#nonsyn   stop    syn 
#2543    108   1269 
table(T_diff$Type)
#nonsyn   stop    syn 
#2685    112   1294 

##ad SC

SC<-read.csv("Output1A/SelCoeff/SC.csv", row.names = 1, stringsAsFactors = F)

SC<-SC[SC$pos>=289& SC$pos<=8613,]
Ts_same<-SC[same.sites,]
Ts_diff<-SC[-same.sites,]

mean(Ts_same$EstSC, na.rm = T) #0.002680869
mean(Ts_diff$EstSC, na.rm = T) #0.002578303

wilcox.test(Ts_same$EstSC, Ts_diff$EstSC,  alternative = "greater", paired = FALSE)
#W = 8210700, p-value = 0.03153

## Coding region only
Ts2<-Ts[Ts$pos>=342& Ts$pos<=8613,]
T_same<-Ts2[same.sites,]
T_diff<-Ts2[-same.sites,]

mean(T_same$mean, na.rm=T)  #0.004779994
mean(T_diff$mean, na.rm=T) #0.00487539
wilcox.test(T_same$mean, T_diff$mean,  alternative = "less", paired = FALSE)
#W = 7837100, p-value = 0.2332

SC2<-SC[SC$pos>=342& SC$pos<=8613,]
Ts_same<-SC2[same.sites,]
Ts_diff<-SC2[-same.sites,]
mean(Ts_same$EstSC, na.rm = T) #0.002621477
mean(Ts_diff$EstSC, na.rm = T) #0.002583721

wilcox.test(Ts_same$EstSC, Ts_diff$EstSC,  alternative = "greater", paired = FALSE)
#W = 8210700, p-value = 0.03736



#####
#dinucleotides
fas1<-read.alignment("Output1A/HCV1A_Consensus195_CDS.fasta",format = "fasta")
fas2<-read.alignment("Output1A/HCV1A_comobined_CDS.fasta",format = "fasta")

filen<-c(195,618)
for (f in 1:2){
        fas<-get(paste0("fas",f))
        freq<-list()
        Rho<-list()
        Zscores<-list()
        for (j in 1: length(fas[[3]])){
                seq<-fas[[3]][[j]]
                seq<-substring(seq, 1:nchar(seq),1:nchar(seq))
                
                #record dinucleotide freq
                freq[[j]]<-as.data.frame(seqinr::count(seq,2))
                names(freq)[j]<-fas[[2]][j]
                Rho[[j]]<-as.data.frame(rho(seq,2))
                names(Rho)[j]<-fas[[2]][j]
                Zscores[[j]]<-as.data.frame(zscore(seq, modele="base"))
                names(Zscores)[j]<-fas[[2]][j]
        }

        dint<-data.frame(diNT=freq[[1]][1])
        Rho.values<-data.frame(diNT=freq[[1]][1])
        Z.scores<-data.frame(diNT=freq[[1]][1])
        for (j in 1:length(freq)) {
                dt1<-freq[[j]]
                dt2<-Rho[[j]]
                dt3<-Zscores[[j]]
                dint<-cbind(dint,dt1[,2])
                Rho.values<-cbind(Rho.values, dt2[,2])
                Z.scores<-cbind(Z.scores,dt3[,2])
        }

        tb1<-data.frame(diNT=dint[,1])
        tb1$Freq<-rowMeans(dint[,2:ncol(dint)])
        tb1$Rho<-rowMeans(Rho.values[,2:ncol(Rho.values)], na.rm=1)
        tb1$Z<-rowMeans(Z.scores[,2:ncol(Z.scores)], na.rm=1)
        
        write.csv(tb1, paste0("Output1A/SummaryStats/Dinucleotides.sum.", filen[f],".csv"))
}

# 195 or 618
tb1<-read.csv(paste0("Output1A/SummaryStats/Dinucleotides.sum.", filen[f],".csv"), stringsAsFactors = F, row.names = 1)

tb1$diNT<-toupper(tb1$diNT)
ggplot(tb1, aes(x=diNT, y=Rho), ylim)+geom_point(color=colors2[4], size=3)+
        ylab("Rho statistic")+theme_bw()+ylim(0.7, 1.3)+
        theme(axis.text.x = element_text(size = 12, angle = 90, hjust=1, color="gray20"))+
        xlab("")+annotate("rect", xmin=0, xmax=9.5,ymax=Inf,ymin=-Inf, fill="gray90", alpha=0.1 , color=NA)+
        geom_hline(yintercept = 1.23, color = "gray30", size=.5,linetype=1)+
        geom_hline(yintercept = 0.79, color = "gray30", size=.5,linetype=1)+
        geom_vline(xintercept = c(1:(nrow(tb1)-1))+0.5, color = "gray70", size=.2)+
        theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank(),
              panel.grid.major.y = element_line(size = 0.2, colour = "grey30", linetype=2))+  
        geom_text(label=paste0("\u2193 under-represented"), x= nrow(tb1)+1, y=0.75, size=2.5, color="gray20",hjust=0)+
        geom_text(label=paste0("\u2191 over-represented"), x= nrow(tb1)+1, y=1.25, size=2.5, color="gray20",hjust=0)+
        coord_cartesian(clip = 'off') +
        theme(plot.margin = unit(c(.5, 4, .5, .5), "cm"))
ggsave("Output1A/SummaryFig.Filtered/Dinuc.rho.plot.618.pdf", height = 4,width = 5.5)
#

####  Calculate among hosts mut freq  ###

align2<-read.fasta("Data/HCV1A_CDS_Alignment.fasta",as.string=TRUE)
align2<-read.fasta("Output1A/HCV1A_Consensus195_CDS.fasta",as.string=TRUE)
align2<-read.fasta("Output1A/HCV1A_comobined_CDS.fasta",as.string=TRUE)

consensus<-read.fasta("Data/HCVref.fasta",as.string=TRUE)
al.ta2<-data.frame(pos=342:(9036+341))
al.ta2<-data.frame(pos=342:(8272+341))

for (i in 1: length(align2)){
        seq<-paste(align2[i])
        sname<-substr(names(align2)[i], start=1,stop =10)
        seq<-unlist(strsplit(seq, "") )
        
        al.ta2[,sname]<-seq
}

df<-al.ta2[,2:ncol(al.ta2)]
ta<-data.frame(al.ta2[,1])
for (i in 1:nrow(ta)){
        vec<-paste0(df[i,])
        ta$a[i]<-length(which(vec=="a"))
        ta$t[i]<-length(which(vec=="t"))
        ta$c[i]<-length(which(vec=="c"))
        ta$g[i]<-length(which(vec=="g"))
}
colnames(ta)[1]<-"pos"
###

ta$MajNt<-apply(ta[,2:5],1,function(x) c("a","t","c","g")[which.max(x)])
reference<-read.dna("Data/HCVref.fasta", format = "fasta",as.character=TRUE)
ref.code<-reference[342:ta$pos[nrow(ta)]]
ta$ref<-ref.code
ta$transition.maj<-NA
ta$transition.ref<-NA
for (j in 1:nrow(ta)) ta$transition.maj[j]<-transition(ta$MajNt[j])        
for (j in 1:nrow(ta)) ta$transition.ref[j]<-transition(ta$ref[j])
total<-ncol(al.ta2)-1
#determine Transition mutation freq of every site.
for (k in 1:nrow(ta)){
        if (is.na(ta$MajNt[k])) {
                ta$freq.Ts[k]<-NA #transition mutations
                ta$freq.Ts.ref[k]<-NA
                
                ta$freq.transv[k]<-NA #transversion mutations
                ta$freq.transv.ref[k]<-NA
                ta$freq.transv1[k]<-NA
                ta$freq.transv2[k]<-NA
                ta$freq.transv1.ref[k]<-NA
                ta$freq.transv2.ref[k]<-NA
                
                ta$freq.mutations.ref[k]<-NA #all mutations
                ta$freq.mutations[k]<-NA
                ta$freq.mutations.indels[k]<-NA
                ta$freq.mutations.indels.ref[k]<-NA
        }
        else {
                MajNum <- ta [k,paste0(ta$MajNt[k])]
                MutNum1<- ta [k,paste0(ta$transition.maj[k])]
                WTNum <- ta [k,paste0(ta$ref[k])]
                MutNum2<- ta [k,paste0(ta$transition.ref[k])]
                
                ta$freq.Ts[k]<-MutNum1/total
                ta$freq.Ts.ref[k]<-MutNum2/total
                
                
                #mutation frequencies of all transversion mutataions
                if (ta$MajNt[k]=="a"|ta$MajNt[k]=='g'){
                        TrvMutNum<-ta[k,"c"]+ta[k,"t"]}
                if (ta$MajNt[k]=="c"|ta$MajNt[k]=="t"){
                        TrvMutNum<-ta[k,"a"]+ta[k,"g"]}
                ta$freq.transv[k]<-TrvMutNum/total
                if (ta$ref[k]=="a"|ta$ref[k]=='g'){
                        TrvMutNum2<-ta[k,"c"]+ta[k,"t"]}
                if (ta$ref[k]=="c"|ta$ref[k]=="t"){
                        TrvMutNum2<-ta[k,"a"]+ta[k,"g"]}
                ta$freq.transv.ref[k]<-TrvMutNum2/total
                
                #Frequenceis for specific transversion mutations (1 & 2)
                Tvs1Num<-ta[k,paste0(transv1(ta$MajNt[k]))]
                Tvs2Num<-ta[k,paste0(transv2(ta$MajNt[k]))]
                ta$freq.transv1[k]<-Tvs1Num/total
                ta$freq.transv2[k]<-Tvs2Num/total
                Tvs1rNum<-ta[k,paste0(transv1(ta$ref[k]))]
                Tvs2rNum<-ta[k,paste0(transv2(ta$ref[k]))]
                ta$freq.transv1.ref[k]<-Tvs1rNum/total
                ta$freq.transv2.ref[k]<-Tvs2rNum/total
                
                
                #Frequencies of all SNPs (no indels)
                AllMutNum<-total-MajNum
                AllMutNum2<-total-WTNum
                
                ta$freq.mutations[k]<-AllMutNum/total
                ta$freq.mutations.ref[k]<-AllMutNum2/total        }
}        


TypeOfSite<-c()
TypeOfSite.tv1<-c()
TypeOfSite.tv2<-c()
for (codon in 1:(nrow(ta)/3)) {#for each codon in the sequence
        positions <- c(codon*3-2,codon*3-1, codon*3)
        WTcodon <- ta$ref[positions]  
        if (is.na(WTcodon[1])|is.na(WTcodon[2])|is.na(WTcodon[3])){ 
                WTcodon<-c('n','n','n')
                mutant1codon<-c('n','n','n')
                mutant2codon<-c('n','n','n')
                mutant3codon<-c('n','n','n')}
        else{                        
                mutant1codon <- c(transition(WTcodon[1]), WTcodon[2:3])  #If the first position has transistion mutation, it's labeld as mutatnt1codon.
                mutant2codon <- c(WTcodon[1],transition(WTcodon[2]), WTcodon[3])
                mutant3codon <- c(WTcodon[1:2], transition(WTcodon[3]))
                
                #transversion mutation to 'a' or 'c'
                mutant1codon.tv1 <- c(transv1(WTcodon[1]), WTcodon[2:3]) 
                mutant2codon.tv1 <- c(WTcodon[1],transv1(WTcodon[2]), WTcodon[3])
                mutant3codon.tv1 <- c(WTcodon[1:2], transv1(WTcodon[3]))
                #transversion mutation to 'g' or 't'
                mutant1codon.tv2 <- c(transv2(WTcodon[1]), WTcodon[2:3])  
                mutant2codon.tv2 <- c(WTcodon[1],transv2(WTcodon[2]), WTcodon[3])
                mutant3codon.tv2 <- c(WTcodon[1:2], transv2(WTcodon[3]))
        }
        
        
        TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant1codon))
        TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant2codon))
        TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant3codon))
        
        TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant1codon.tv1))
        TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant2codon.tv1))
        TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant3codon.tv1))
        
        TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant1codon.tv2))
        TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant2codon.tv2))
        TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant3codon.tv2))                
} # This creates a vector showing if there is a transition/transversion mutation at a particular codon, &

ta$Type<-TypeOfSite[1:length(ta$pos)]
ta$Type.tv1<-TypeOfSite.tv1[1:length(ta$pos)]
ta$Type.tv2<-TypeOfSite.tv2[1:length(ta$pos)]

mutrates<-read.csv("Data/Geller.mutation.rates.csv",stringsAsFactors=FALSE)

ta$TSmutrate[ta$ref=="a"]<-mutrates$mut.rate[mutrates$mutations=="AG"]
ta$TSmutrate[ta$ref=="c"]<-mutrates$mut.rate[mutrates$mutations=="CU"]
ta$TSmutrate[ta$ref=="g"]<-mutrates$mut.rate[mutrates$mutations=="GA"]
ta$TSmutrate[ta$ref=="t"]<-mutrates$mut.rate[mutrates$mutations=="UC"]

#ta$TSmutrate.hiv[ta$ref=="a"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="AG"]
#ta$TSmutrate.hiv[ta$ref=="c"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="CU"]
#ta$TSmutrate.hiv[ta$ref=="g"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="GA"]
#ta$TSmutrate.hiv[ta$ref=="t"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="UC"]

ta$TVSmutrate.tv1[ta$ref=="a"]<-mean(mutrates$mut.rate[mutrates$mutations=="AC"])
ta$TVSmutrate.tv1[ta$ref=="c"]<-mean(mutrates$mut.rate[mutrates$mutations=="CA"])
ta$TVSmutrate.tv1[ta$ref=="g"]<-mean(mutrates$mut.rate[mutrates$mutations=="GC"])
ta$TVSmutrate.tv1[ta$ref=="t"]<-mean(mutrates$mut.rate[mutrates$mutations=="UA"])

ta$TVSmutrate.tv2[ta$ref=="a"]<-mean(mutrates$mut.rate[mutrates$mutations=="AU"])
ta$TVSmutrate.tv2[ta$ref=="c"]<-mean(mutrates$mut.rate[mutrates$mutations=="CG"])
ta$TVSmutrate.tv2[ta$ref=="g"]<-mean(mutrates$mut.rate[mutrates$mutations=="GU"])
ta$TVSmutrate.tv2[ta$ref=="t"]<-mean(mutrates$mut.rate[mutrates$mutations=="UG"])

ta$TVSmutrate.tvs[ta$ref=="a"]<-mean(mutrates$mut.rate[mutrates$mutations=="AU"], mutrates$mut.rate[mutrates$mutations=="AC"])
ta$TVSmutrate.tvs[ta$ref=="c"]<-mean(mutrates$mut.rate[mutrates$mutations=="CG"],mutrates$mut.rate[mutrates$mutations=="CA"])
ta$TVSmutrate.tvs[ta$ref=="g"]<-mean(mutrates$mut.rate[mutrates$mutations=="GU"],mutrates$mut.rate[mutrates$mutations=="GC"])
ta$TVSmutrate.tvs[ta$ref=="t"]<-mean(mutrates$mut.rate[mutrates$mutations=="UG"],mutrates$mut.rate[mutrates$mutations=="UA"])


for (k in 1:length(ta$pos)){
        ta$EstSelCoeff[k] <- EstimatedS(ta$TSmutrate[k],ta[k,colnames(ta)=='freq.Ts.ref'])
        #ta$EstSelCoeff_hiv[k] <- EstimatedS(ta$TSmutrate.hiv[k],ta[k,colnames(ta)=='freq.Ts.ref'])
        ta$EstSelCoeff_transv[k] <- EstimatedS(ta$TVSmutrate.tvs[k],ta[k,colnames(ta)=='freq.transv.ref'])
        ta$EstSelCoeff_trans1[k] <- EstimatedS(ta$TVSmutrate.tv1[k],ta[k,colnames(ta)=='freq.transv1.ref'])
        ta$EstSelCoeff_trans2[k] <- EstimatedS(ta$TVSmutrate.tv2[k],ta[k,colnames(ta)=='freq.transv2.ref'])
        
        if (k%%3==1){
                if (is.na(ta$MajNt[k])|is.na(ta$MajNt[k+1])|is.na(ta$MajNt[k+2])) { ta$MajAA[k]<-"NA"
                ta$WTAA[k]<-"NA"
                ta$MUTAA[k]<-"NA"
                ta$TVS1_AA[k]<-"NA"
                ta$TVS2_AA[k]<-"NA"}
                else {  ta$MajAA[k] = seqinr::translate(ta$MajNt[c(k,k+1,k+2)])
                ta$WTAA[k] = seqinr::translate(ta$ref[c(k,k+1,k+2)])
                ta$MUTAA[k] = seqinr::translate(c(transition(ta$ref[k]),ta$ref[c(k+1,k+2)]))
                ta$TVS1_AA[k] = seqinr::translate(c(transv1(ta$ref[k]),ta$ref[c(k+1,k+2)]))
                ta$TVS2_AA[k] = seqinr::translate(c(transv2(ta$ref[k]),ta$ref[c(k+1,k+2)]))}
        } 
        if (k%%3==2){
                if (is.na(ta$MajNt[k-1])|is.na(ta$MajNt[k])|is.na(ta$MajNt[k+1]))  {ta$MajAA[k]<-"NA"
                ta$WTAA[k]<-"NA"
                ta$MUTAA[k]<-"NA"
                ta$TVS1_AA[k]<-"NA"
                ta$TVS2_AA[k]<-"NA"}
                else {  ta$MajAA[k] = seqinr::translate(ta$MajNt[c(k-1,k,k+1)])
                ta$WTAA[k] = seqinr::translate(ta$ref[c(k-1,k,k+1)])
                ta$MUTAA[k] = seqinr::translate(c(ta$ref[c(k-1)],transition(ta$ref[k]),ta$ref[c(k+1)]))
                ta$TVS1_AA[k] = seqinr::translate(c(ta$ref[c(k-1)],transv1(ta$ref[k]),ta$ref[c(k+1)]))
                ta$TVS2_AA[k] = seqinr::translate(c(ta$ref[c(k-1)],transv2(ta$ref[k]),ta$ref[c(k+1)]))}
        }
        if (k%%3==0){
                if (is.na(ta$MajNt[k-2])|is.na(ta$MajNt[k-1])|is.na(ta$MajNt[k]))  {  ta$MajAA[k]<-"NA"
                ta$WTAA[k]<-"NA"
                ta$MUTAA[k]<-"NA"
                ta$TVS1_AA[k]<-"NA"
                ta$TVS2_AA[k]<-"NA"}
                else {  ta$MajAA[k] = seqinr::translate(ta$MajNt[c(k-2,k-1,k)])
                ta$WTAA[k] = seqinr::translate(ta$ref[c(k-2,k-1,k)])
                ta$MUTAA[k] = seqinr::translate(c(ta$ref[c(k-2,k-1)],transition(ta$ref[k])))
                ta$TVS1_AA[k] = seqinr::translate(c(ta$ref[c(k-2,k-1)],transv1(ta$ref[k])))
                ta$TVS2_AA[k] = seqinr::translate(c(ta$ref[c(k-2,k-1)],transv2(ta$ref[k])))}
}}
#Add whether AA change is drastic & makes CpG
ta$bigAAChange<-0
ta$bigAAChange.tv1<-0
ta$bigAAChange.tv2<-0
ta$makesCpG <- 0
ta$makesCpG.tvs <- 0
ta$makesCpG.tv1 <- 0
ta$makesCpG.tv2 <- 0
#ta$makesCpG_all <- 0
#ta$color<-""
for(j in 2:nrow(ta)-1){
        WT <- amCat(ta[j,'WTAA'])
        MUT <- amCat(ta[j,'MUTAA'])
        MUT1<-amCat(ta[j,'TVS1_AA'])
        MUT2<-amCat(ta[j,'TVS2_AA'])
        
        if (WT != MUT) ta$bigAAChange[j] <- 1
        if (WT != MUT1) ta$bigAAChange.tv1[j] <- 1
        if (WT != MUT2) ta$bigAAChange.tv2[j] <- 1
        
        trip <- ta$ref[c(j-1, j,j+1)]
        if (is.na(trip[1])|is.na(trip[2])|is.na(trip[3])) 
                next
        else{
                if (trip[1] == "c" & trip[2] == "a" ) ta$makesCpG[j] <- 1 
                if (trip[2] == "t" & trip[3] == "g")  ta$makesCpG[j] <- 1
                if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) ta$makesCpG.tvs[j] <- 1
                if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) ta$makesCpG.tvs[j] <- 1
                
                if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) ta$makesCpG.tv2[j] <- 1                                
                if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) ta$makesCpG.tv1[j] <- 1
                
        }
}


write.csv(ta,"Output1A/NCBI_mutfreq.csv")

write.csv(ta,"Output1A/HCV1A_195_mutfreq.csv")

write.csv(ta,"Output1A/HCV1A_Combined_mutfreq.csv")



#calculate means

#1. NCBI 423 SEQUENCES
mean(ta$freq.Ts)  #0.03995104
mean(ta$freq.Ts.ref)  #0.051109464
mean(ta$freq.transv1)  #0.004104151
mean(ta$freq.transv2)  #0.002803862
mean(ta$freq.transv) #0.006908013

#2. My 195 SEQUENCES
mean(ta$freq.Ts)  #0.04391026
mean(ta$freq.Ts.ref)  #0.05600853
mean(ta$freq.transv1)  #0.004454322
mean(ta$freq.transv2)  #0.003167311
mean(ta$freq.transv)   #0.007621634


#3. Combined 618 SEQUENCES
mean(ta$freq.Ts)  #0.04223532
mean(ta$freq.Ts.ref)  #0.0535479
mean(ta$freq.transv1)  #0.004369441
mean(ta$freq.transv2)  #0.003038871
mean(ta$freq.transv)   #0.007408312

