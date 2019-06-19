#Script to analyse the frequency data and associate with features
library(dplyr)
source("Rscripts/baseRscript.R")

#get the file name
HCVFiles_SeqData<-list.files("Output1A/SeqData/",pattern="SeqData")

Overview<-list()
for (i in 1:length(HCVFiles_SeqData)){   
        id<-substr(paste(HCVFiles_SeqData[i]),start=9,stop=15)
        print(id)
        OverviewDF<-read.csv(paste0("Output1A/SeqData/",HCVFiles_SeqData[i]),stringsAsFactors=FALSE)
        OverviewDF<-OverviewDF[,-1]
                
        TypeOfSite<-c()
        TypeOfSite.tv1<-c()
        TypeOfSite.tv2<-c()
        for (codon in 1:(nrow(OverviewDF)/3)) {#for each codon in the sequence
                positions <- c(codon*3-2,codon*3-1, codon*3)
                WTcodon <- OverviewDF$MajNt[positions]  
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
                  # whehter the mutation will be Syn, nonSyn, or Stop codon.

        OverviewDF$Type<-TypeOfSite[1:length(OverviewDF$pos)]
        OverviewDF$Type.tv1<-TypeOfSite.tv1[1:length(OverviewDF$pos)]
        OverviewDF$Type.tv2<-TypeOfSite.tv2[1:length(OverviewDF$pos)]
        #OverviewDF$consensus<-ref[1:length(OverviewDF$pos)]

        
        Overview[[i]]<-OverviewDF[,-c(1:6)]
        names(Overview)[i]<-id   
}




###############################
#Mut rates from Geller paper 

## 2. Using 12 different Mutation Frequencies baed on from/to nucleotides

mutrates<-read.csv("Data/Geller.mutation.rates.csv",stringsAsFactors=FALSE)
#Mut rates and sel coefficients from Abrams paper (for HIV)
#read.csv("~/programs/BachelerProject_New/Data/HIVMutRates/HIVMutRates.csv")->mutrates1

Overview_sum<-list()
for (i in 1:length(Overview)){

        OverviewDF<-Overview[[i]]
        id<-substr(paste(HCVFiles_SeqData[i]),start=9,stop=15)

        OverviewDF$TSmutrate[OverviewDF$MajNt=="a"]<-mutrates$mut.rate[mutrates$mutations=="AG"]
        OverviewDF$TSmutrate[OverviewDF$MajNt=="c"]<-mutrates$mut.rate[mutrates$mutations=="CU"]
        OverviewDF$TSmutrate[OverviewDF$MajNt=="g"]<-mutrates$mut.rate[mutrates$mutations=="GA"]
        OverviewDF$TSmutrate[OverviewDF$MajNt=="t"]<-mutrates$mut.rate[mutrates$mutations=="UC"]
        
        #OverviewDF$TSmutrate.hiv[OverviewDF$MajNt=="a"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="AG"]
        #OverviewDF$TSmutrate.hiv[OverviewDF$MajNt=="c"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="CU"]
        #OverviewDF$TSmutrate.hiv[OverviewDF$MajNt=="g"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="GA"]
        #OverviewDF$TSmutrate.hiv[OverviewDF$MajNt=="t"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="UC"]

        
        OverviewDF$TVSmutrate1[OverviewDF$MajNt=="a"]<-mutrates$mut.rate[mutrates$mutations=="AC"]
        OverviewDF$TVSmutrate1[OverviewDF$MajNt=="c"]<-mutrates$mut.rate[mutrates$mutations=="CA"]
        OverviewDF$TVSmutrate1[OverviewDF$MajNt=="g"]<-mutrates$mut.rate[mutrates$mutations=="GC"]
        OverviewDF$TVSmutrate1[OverviewDF$MajNt=="t"]<-mutrates$mut.rate[mutrates$mutations=="UA"]
        
        OverviewDF$TVSmutrate2[OverviewDF$MajNt=="a"]<-mutrates$mut.rate[mutrates$mutations=="AU"]
        OverviewDF$TVSmutrate2[OverviewDF$MajNt=="c"]<-mutrates$mut.rate[mutrates$mutations=="CG"]
        OverviewDF$TVSmutrate2[OverviewDF$MajNt=="g"]<-mutrates$mut.rate[mutrates$mutations=="GU"]
        OverviewDF$TVSmutrate2[OverviewDF$MajNt=="t"]<-mutrates$mut.rate[mutrates$mutations=="UG"]
        
        OverviewDF$TVSmutrate.tvs[OverviewDF$MajNt=="a"]<-mean(mutrates$mut.rate[mutrates$mutations=="AU"], mutrates$mut.rate[mutrates$mutations=="AC"])
        OverviewDF$TVSmutrate.tvs[OverviewDF$MajNt=="c"]<-mean(mutrates$mut.rate[mutrates$mutations=="CG"],mutrates$mut.rate[mutrates$mutations=="CA"])
        OverviewDF$TVSmutrate.tvs[OverviewDF$MajNt=="g"]<-mean(mutrates$mut.rate[mutrates$mutations=="GU"],mutrates$mut.rate[mutrates$mutations=="GC"])
        OverviewDF$TVSmutrate.tvs[OverviewDF$MajNt=="t"]<-mean(mutrates$mut.rate[mutrates$mutations=="UG"],mutrates$mut.rate[mutrates$mutations=="UA"])

        OverviewDF$MajAA<-""
        OverviewDF$WTAA<-""
        OverviewDF$MUTAA<-""
        OverviewDF$TVS1_AA<-""
        OverviewDF$TVS2_AA<-""
        
        
        for (k in 1:length(OverviewDF$pos)){
                OverviewDF$EstSelCoeff[k] <- EstimatedS(OverviewDF$TSmutrate[k],OverviewDF[k,colnames(OverviewDF)=='freq.Ts'])
                #OverviewDF$EstSelCoeff_hiv[k] <- EstimatedS(OverviewDF$TSmutrate.hiv[k],OverviewDF[k,colnames(OverviewDF)=='freq.Ts'])
                OverviewDF$EstSelCoeff_T1[k] <- EstimatedS(OverviewDF$TVSmutrate1[k],OverviewDF[k,colnames(OverviewDF)=='freq.transv1'])
                OverviewDF$EstSelCoeff_T2[k] <- EstimatedS(OverviewDF$TVSmutrate2[k],OverviewDF[k,colnames(OverviewDF)=='freq.transv2'])
                OverviewDF$EstSelCoeff_transv[k] <- EstimatedS(OverviewDF$TVSmutrate.tvs[k],OverviewDF[k,colnames(OverviewDF)=='freq.transv'])
                
                if (k%%3==1){
                        if (is.na(OverviewDF$MajNt[k])|is.na(OverviewDF$MajNt[k+1])|is.na(OverviewDF$MajNt[k+2])) { 
                                OverviewDF$MajAA[k]<-"NA"
                                OverviewDF$WTAA[k]<-"NA"
                                OverviewDF$MUTAA[k]<-"NA"
                                OverviewDF$TVS1_AA[k]<-"NA"
                                OverviewDF$TVS2_AA[k]<-"NA"}
                        else {  OverviewDF$MajAA[k] = seqinr::translate(OverviewDF$MajNt[c(k,k+1,k+2)])
                        OverviewDF$WTAA[k] = seqinr::translate(OverviewDF$MajNt[c(k,k+1,k+2)])
                        OverviewDF$MUTAA[k] = seqinr::translate(c(transition(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1,k+2)]))
                        OverviewDF$TVS1_AA[k] = seqinr::translate(c(transv1(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1,k+2)]))
                        OverviewDF$TVS2_AA[k] = seqinr::translate(c(transv2(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1,k+2)]))}
                } 
                if (k%%3==2){
                        if (is.na(OverviewDF$MajNt[k-1])|is.na(OverviewDF$MajNt[k])|is.na(OverviewDF$MajNt[k+1]))  {OverviewDF$WTAA[k]<-"NA"
                        OverviewDF$MajAA[k]<-"NA"
                        OverviewDF$MUTAA[k]<-"NA"
                        OverviewDF$TVS1_AA[k]<-"NA"
                        OverviewDF$TVS2_AA[k]<-"NA"}
                        else {  OverviewDF$MajAA[k] = seqinr::translate(OverviewDF$MajNt[c(k-1,k,k+1)])
                        OverviewDF$WTAA[k] = seqinr::translate(OverviewDF$MajNt[c(k-1,k,k+1)])
                        OverviewDF$MUTAA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-1)],transition(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1)]))
                        OverviewDF$TVS1_AA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-1)],transv1(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1)]))
                        OverviewDF$TVS2_AA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-1)],transv2(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1)]))}
                }
                if (k%%3==0){
                        if (is.na(OverviewDF$MajNt[k-2])|is.na(OverviewDF$MajNt[k-1])|is.na(OverviewDF$MajNt[k]))  {  
                                OverviewDF$MajAA[k]<-"NA"
                                OverviewDF$WTAA[k]<-"NA"
                                OverviewDF$MUTAA[k]<-"NA"
                                OverviewDF$TVS1_AA[k]<-"NA"
                                OverviewDF$TVS2_AA[k]<-"NA"}
                        else {  OverviewDF$MajAA[k] = seqinr::translate(OverviewDF$MajNt[c(k-2,k-1,k)])
                        OverviewDF$WTAA[k] = seqinr::translate(OverviewDF$MajNt[c(k-2,k-1,k)])
                        OverviewDF$MUTAA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-2,k-1)],transition(OverviewDF$MajNt[k])))
                        OverviewDF$TVS1_AA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-2,k-1)],transv1(OverviewDF$MajNt[k])))
                        OverviewDF$TVS2_AA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-2,k-1)],transv2(OverviewDF$MajNt[k])))}
                }
        }
        
        
        #Add whether AA change is drastic & makes CpG
        OverviewDF$bigAAChange<-0
        OverviewDF$bigAAChange.ref<-0
        OverviewDF$bigAAChange.tv1<-0
        OverviewDF$bigAAChange.tv2<-0
        OverviewDF$makesCpG <- 0
        OverviewDF$makesCpG.ref <- 0
        OverviewDF$makesCpG.tvs <- 0
        OverviewDF$makesCpG.tv1 <- 0
        OverviewDF$makesCpG.tv2 <- 0
        OverviewDF$makesCpG_all <- 0
        #OverviewDF$color<-""
        for(j in 2:nrow(OverviewDF)-1){
                MAJ <- amCat(OverviewDF[j,'MajAA'])
                WT <- amCat(OverviewDF[j,'WTAA'])
                MUT <- amCat(OverviewDF[j,'MUTAA'])
                MUT1<-amCat(OverviewDF[j,'TVS1_AA'])
                MUT2<-amCat(OverviewDF[j,'TVS2_AA'])
                
                if (MAJ!= MUT) OverviewDF$bigAAChange[j] <- 1
                if (WT != MUT) OverviewDF$bigAAChange.ref[j] <- 1
                if (WT != MUT1) OverviewDF$bigAAChange.tv1[j] <- 1
                if (WT != MUT2) OverviewDF$bigAAChange.tv2[j] <- 1
                
                trip <- OverviewDF$MajNt[c(j-1, j,j+1)]
                if (is.na(trip[1])|is.na(trip[2])|is.na(trip[3])) 
                        next
                else{
                        if (trip[1] == "c" & trip[2] == "a" ) OverviewDF$makesCpG[j] <- 1 
                        if (trip[2] == "t" & trip[3] == "g")  OverviewDF$makesCpG[j] <- 1
                        if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) OverviewDF$makesCpG.tvs[j] <- 1
                        if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) OverviewDF$makesCpG.tvs[j] <- 1
                        
                        if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) OverviewDF$makesCpG.tv2[j] <- 1                                
                        if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) OverviewDF$makesCpG.tv1[j] <- 1
                        if (trip[1] == "c" & trip[2] != "g") OverviewDF$makesCpG_all[j] <- 1
                        if (trip[2] != "c" & trip[3] == "g") OverviewDF$makesCpG_all[j] <- 1}
                
                trip2 <- OverviewDF$ref[c(j-1, j,j+1)]
                if (is.na(trip2[1])|is.na(trip2[2])|is.na(trip2[3])) 
                        next
                else{
                        if (trip2[1] == "c" & trip2[2] == "a" ) OverviewDF$makesCpG.ref[j] <- 1 
                        if (trip2[2] == "t" & trip2[3] == "g")  OverviewDF$makesCpG.ref[j] <- 1
                }
        }
        
                
        write.csv(OverviewDF,paste0("Output1A/Overview1/",id,"_overview.csv"))
        Overview_sum[[i]]<-OverviewDF
        names(Overview_sum)[i]<-id
        print(id)
}        
