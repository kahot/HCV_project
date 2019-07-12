#######
# Calculate the Fst between genotyeps

Overview1A<-Overview1.1A[1:10]
Overview1B<-Overview1.1B[1:10]
Overview3A<-Overview1.3A[1:10]
mergepos<-read.csv("Data/MergedPositionInfo.csv", row.names = 1)
mergepos<-mergepos[265:8610,]
geno<-c("1A","1B","3A")

for ( g in 1:3){
        OverV<-get(paste0("Overview",geno[g]))
        
        OvervDf<-list()
        for ( i in 1:length(OverV)){
                df<-OverV[[i]]
                n<-which(colnames(mergepos)==paste0("org.pos.",geno[g]))
                df<-merge(mergepos,df, by.x=colnames(mergepos)[n], by.y="pos",all.x = T)
                OvervDf[[i]]<-df
                names(OvervDf)[i]<-names(OverV)[i]
        }
        listname<-paste0("OverV.", geno[g])
        assign(listname, OvervDf)
}

O<-do.call(c, list(OverV.1A,OverV.1B,OverV.3A))
#O<-append(OverV.1A,OverV.1B)
        
ids<-names(O)
Comb<-t(combn(ids,2))
FstDF<-data.frame(matrix(ncol=0,nrow=nrow(Comb)))

Pi<-read.csv(paste0("Output_all/Diversity/NucDiversity.per.sample2.1A.csv"), stringsAsFactors = F, row.names = 1)
for(g in 2:3){
        dt<-read.csv(paste0("Output_all/Diversity/NucDiversity.per.sample2.",geno[g],".csv"), stringsAsFactors = F, row.names = 1)
        Pi<-rbind(Pi, dt)
}


for (i in 246:nrow(Comb)){
        samp1<-Comb[i,1]
        samp2<-Comb[i,2]
        df1<-O[[samp1]]
        df2<-O[[samp2]]
        fname<-paste(samp1,samp2)
        FstDF$comb[i]<-fname
        
        df1$DT<-""
        for (k in 1:nrow(df1)){
                if (is.na(df1$freq.Ts[k])|is.na(df2$freq.Ts[k])){
                        df1$DT[k]<-NA
                }
                else{
                        ac<-(df1[k, "a"]+df2[k, "a"])*(df1[k,"c"]+df2[k,"c"])
                        ag<-(df1[k, "a"]+df2[k, "a"])*(df1[k,"g"]+df2[k,"g"])
                        at<-(df1[k, "a"]+df2[k, "a"])*(df1[k,"t"]+df2[k,"t"])
                        cg<-(df1[k, "c"]+df2[k, "c"])*(df1[k,"g"]+df2[k,"g"])
                        ct<-(df1[k, "c"]+df2[k, "c"])*(df1[k,"t"]+df2[k,"t"])
                        gt<-(df1[k, "g"]+df2[k, "g"])*(df1[k,"t"]+df2[k,"t"])
                        m<-((df1[k,"TotalReads"]+df2[k,"TotalReads"])^2-df1[k,"TotalReads"]-df2[k,"TotalReads"])/2
                        df1$DT[k]<-(ac+ag+at+cg+ct+gt)/m
                }
        }
        df1$DT<-as.numeric(df1$DT)
        n<-nrow(df1[!is.na(df1$DT),])
        FstDF$piT[i]<-sum(df1$DT, na.rm=T)/n
        FstDF$piS.ave[i]<-mean(c(Pi$Pi[Pi$SampleID==samp1],Pi$Pi[Pi$SampleID==samp2]),na.rm=T )
        FstDF$Fst[i]<-(FstDF$piT[i]- FstDF$piS.ave[i])/ FstDF$piT[i]
        print(i)
}

write.csv(FstDF,paste0("Output_all/Diversity/AllGenotype_Fst2.csv"))


id1A<-names(Overview1.1A)
id1B<-names(Overview1.1B)
id3A<-names(Overview1.3A)

data<-FstDF
splitted<-strsplit(data$comb, split=" ")
for (i in 1:nrow(data)){
        data$Sample1[i]<-splitted[[i]][1]
        data$Sample2[i]<-splitted[[i]][2]
}

dat<-data[,4:6]
dat1<-dat
for (i in 1:nrow(dat)){
        if (dat$Sample1[i] %in% id1A) dat$Sample1[i]<-paste0("1A.",dat$Sample1[i])
        if (dat$Sample1[i] %in% id1B) dat$Sample1[i]<-paste0("1B.",dat$Sample1[i])
        if (dat$Sample1[i] %in% id3A) dat$Sample1[i]<-paste0("3A.",dat$Sample1[i])
        if (dat$Sample2[i] %in% id1A) dat$Sample2[i]<-paste0("1A.",dat$Sample2[i])
        if (dat$Sample2[i] %in% id1B) dat$Sample2[i]<-paste0("1B.",dat$Sample2[i])
        if (dat$Sample2[i] %in% id3A) dat$Sample2[i]<-paste0("3A.",dat$Sample2[i])
}

colorFst<-c(rep(cols2[1], times=10), rep(cols2[2], times=10), rep(cols2[3], times=10))
colorFst2<-c(rep(cols2[1], times=9), rep(cols2[2], times=10), rep(cols2[3], times=10))

ggplot(data=dat, aes(Sample1,Sample2, fill=Fst)) +geom_tile(color="gray40")+ 
        scale_fill_gradient2(low="white", high="#0000FFCC",mid="#00FDFFCC", midpoint=0.5,limit=c(0,1))+
        theme_light()+
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),panel.grid.major = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.text.x=element_text(size=9, angle=45,hjust=1))+
        theme(axis.text.x=element_text(color=colorFst))+
        theme(axis.text.y=element_text(color=colorFst2))

ggsave(file=paste0("Output_all/Diversity/Fst_3Gentypes.Plot.pdf"),width=9, height=9, units='in',device='pdf')

