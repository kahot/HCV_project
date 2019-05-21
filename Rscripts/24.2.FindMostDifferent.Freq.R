library(stringr)
library(tidyverse)
library(zoo)

#############
# find the sites that are most different in mutation frequencies between the genomes

mvf<-read.csv("Output_all/Unfiltered/Ts_Mean_3genotypes.csv", stringsAsFactors = F )
mvf<-mvf[,-1]
mvf<-mvf[c(264:8685),]


mvf$diff.1A.3A<-mvf$mean.1A-mvf$mean.3A

hist(mvf$diff.1A.3A)

#find the sites with 3A having higher mut greq (mvf.H) and 1A having the higher mut. freq (mvf.L)
mvf.H<-mvf[order(mvf$diff.1A.3A, decreasing = T,na.last = T),]
mvf.L<-mvf[order(mvf$diff.1A.3A,decreasing=F,na.last=T),]

#top 5% sites
s<-as.integer(0.05*nrow(mvf))

mvf.H<-mvf.H[c(1:s),]
mvf.L<-mvf.L[c(1:s),]


table(mvf.H$gene)
#   Core      E1      E2    HVR1 NS1(P7)     NS2     NS3    NS4A    NS4B    NS5A    NS5B 
#     13      25      46       9       7      35      91      12      50      85      48 

table(mvf.L$gene)
# Core      E1      E2    HVR1 NS1(P7)     NS2     NS3    NS4A    NS4B    NS5A    NS5B 
#   16      27      45      10      13      35      83      11      49      79      53 


#absolute difference top 5%
mvf$diff.1a.3a.abs<-abs(mvf$diff.1A.3A)
mvf.top<-mvf[order(mvf$diff.1a.3a.abs, decreasing = T,na.last = T),]
mvf.top<-mvf.top[c(1:s),]

write.csv(mvf.top,"Output_all/Unfiltered/Top.Difference.1A-3A.csv")

table(mvf.top$gene)
#   Core      E1      E2    HVR1 NS1(P7)     NS2     NS3    NS4A    NS4B    NS5A    NS5B 
#     12      21      41      14      11      28      84      11      46      92      61 



mvfh.core<-mvf.H[mvf.H$gene=="Core",]

mvf.L[mvf.L$gene=="Core",]
