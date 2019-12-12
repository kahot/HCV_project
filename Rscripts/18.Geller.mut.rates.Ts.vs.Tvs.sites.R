# Explore Transv > Trans sites in vitro

Ts<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",row.names = 1,stringsAsFactors = F)
Tv1.MutFreq<-read.csv("Output1A/MutFreq.filtered/Filtered.Tv1.MutFreq.Q35.csv",row.names = 1,stringsAsFactors = F)
Tv2.MutFreq<-read.csv("Output1A/MutFreq.filtered/Filtered.Tv2.MutFreq.Q35.csv",row.names = 1,stringsAsFactors = F)
Tvs.MutFreq<-read.csv("Output1A/MutFreq.filtered/Filtered.Tvs.MutFreq.Q35.csv",row.names = 1,stringsAsFactors = F)
AllMutFreq<-read.csv("Output1A/MutFreq.filtered/Filtered.AllMutFreq.Q35.csv",row.names = 1,stringsAsFactors = F)



S4080Ts<-mfs[mfs$pos==4080,]
S4080All<-AllMutFreq[AllMutFreq$pos==4080,]
S4080Tv1<-Tv1.MutFreq[Tv1.MutFreq$pos==4080,]
S4080Tv2<-Tv2.MutFreq[Tv2.MutFreq$pos==4080,]
S4080Tvs<-Tvs.MutFreq[Tvs.MutFreq$pos==4080,]

S4080Ts$mean    #0.00132958
S4080All$mean   #0.002376884
S4080Tvs$mean   #0.001047295
S4080Tv1$mean   #0.001017858
S4080Tv2$mean   #2.943693e-05

S4080Ts$ref #c

#Transition > Transversions at Pos 4080


geller<-read.csv("Data/Geller.mutation.rates.csv")
geller<-geller[c(2,3,4,5,7,8,9,10,12,17,18,19),]
ggplot(geller, aes(x= mutations, y=mut.rate))+geom_point()

       