library(plotrix)
library(reshape)
library(tidyverse)
library(zoo)
library(purrr)
library(colorspace)
colors2<-qualitative_hcl(6, palette="Dark3")


#Create the filtered mut frequency table of all samples 

HCVFiles_overview3<-list.files("Output1A/Overview3/",pattern="overview3.csv")
FilteredOverview2<-list()
for (i in 1:length(HCVFiles_overview3)){ 
        overviews2<-read.csv(paste0("Output1A/Overview3/",HCVFiles_overview3[i]),stringsAsFactors=FALSE, row.names=1)
        FilteredOverview2[[i]]<-overviews2
        names(FilteredOverview2)[i]<-substr(paste(HCVFiles_overview3[i]),start=1,stop=7)
}

##################################
MutFreq_Ts<-list()
MutFreq_tv1<-list()
MutFreq_tv2<-list()
MutFreq_tvs<-list()
MutFreq_all<-list()

for (i in 1:length(FilteredOverview2)){
        dat<-FilteredOverview2[[i]]
        filename<-names(FilteredOverview2)[i]
        
        MutFreq_Ts[[i]]<-dat[,c("pos","freq.Ts.ref")] 
        MutFreq_tv1[[i]]<-dat[,c("pos","freq.transv1.ref")] 
        MutFreq_tv2[[i]]<-dat[,c("pos","freq.transv2.ref")] 
        MutFreq_tvs[[i]]<-dat[,c("pos","freq.transv.ref")] 
        MutFreq_all[[i]]<-dat[,c("pos","freq.mutations.ref")] 
        
        names(MutFreq_Ts)[i]<-filename
        names(MutFreq_tv1)[i]<-filename
        names(MutFreq_tv2)[i]<-filename
        names(MutFreq_tvs)[i]<-filename
        names(MutFreq_all)[i]<-filename
        
        
}
#assign column names for the list
for (i in 1:length(MutFreq_Ts)) {
        colnames(MutFreq_Ts[[i]])<-c("pos",paste0(names(MutFreq_Ts[i])))
        colnames(MutFreq_tv1[[i]])<-c("pos",paste0(names(MutFreq_tv1[i])))
        colnames(MutFreq_tv2[[i]])<-c("pos",paste0(names(MutFreq_tv2[i])))
        colnames(MutFreq_tvs[[i]])<-c("pos",paste0(names(MutFreq_tvs[i])))
        colnames(MutFreq_all[[i]])<-c("pos",paste0(names(MutFreq_all[i])))
}

Ts<-MutFreq_Ts%>% purrr::reduce(full_join, by='pos')
Tv1.MutFreq<-MutFreq_tv1 %>% purrr::reduce(full_join, by='pos')
Tv2.MutFreq<-MutFreq_tv2 %>% purrr::reduce(full_join, by='pos')
Tvs.MutFreq<-MutFreq_tvs %>% purrr::reduce(full_join, by='pos')
AllMutFreq<-MutFreq_all %>% purrr::reduce(full_join, by='pos')

### all mut. freq with metadata ###

files<-c("Ts","Tv1.MutFreq","Tv2.MutFreq","Tvs.MutFreq","AllMutFreq" )
cnames<-c("",".tv1",".tv2",".tvs","")
mf.files<-list()
s<-length(FilteredOverview2)
M<-FilteredOverview2[[3]]
colnames(M)[33:35]<-c("MutAA","MutAA.tv1","MutAA.tv2")
#muttypes1<-M[,c("pos","ref","Type","Type.tv1","Type.tv2","WTAA","MUTAA","TVS1_AA","TVS2_AA","makesCpG","makesCpG.tv1","makesCpG.tv2","bigAAChange","bigAAChange.tv1","bigAAChange.tv2")]
#muttypes<-M[,c(1,4,5,20,23,27,32,36,39)]

for (i in 1:5) {
        dat<-get(files[i])
        dat$mean<-rowMeans(dat[2:(s+1)],na.rm=T, dims=)
        if (i==1|i==5) muttypes<-M[,c("pos","ref", "Type","WTAA","MutAA","makesCpG","bigAAChange")]
        if(i==2|i==3) muttypes<-M[,c("pos","ref", paste0("Type",cnames[i]),"WTAA",paste0("MutAA",cnames[i]), paste0("makesCpG",cnames[i]), paste0("bigAAChange",cnames[i]))]
        if(i==4) muttypes<-M[,c("pos","ref","Type.tv1","Type.tv2","WTAA","MutAA.tv1","MutAA.tv2","makesCpG.tv1","makesCpG.tv1","bigAAChange.tv1","bigAAChange.tv2")]
        
        dat2<-merge(dat,muttypes,by="pos")
        mf.files[[i]]<-dat2
        names(mf.files)[i]<-files[i]
        #dfname<-paste0(files[i],".F")
        #assign(dfname, dat2)
        write.csv(dat2,paste0("Output1A/MutFreq.filtered/Filtered.",files[i],".Q35.csv"))
}



##############################
##############################
#read in files
files<-c("Ts","Tv1.MutFreq","Tv2.MutFreq","Tvs.MutFreq","AllMutFreq")
mf.files<-list()
for (i in 1:5){
        data<-read.csv(paste0("Output1A/MutFreq.filtered/Filtered.", files[i],".Q35.csv"),row.names = 1,stringsAsFactors = F)
        assign(paste0(files[i]),data)
        mf.files[[i]]<-data
        names(mf.files)[i]<-files[i]
}



##### Chceck the mean mut freq
#Summarize the mean and se for all types of mutations

tb<-data.frame(type=c("Ts","Ts.syn","Ts.ns","Ts.stop", "Tv1","Tv1.syn","Tv1.ns","Tv1.stop","Tv2","Tv2.syn","Tv2.ns","Tv2.stop","Tvs","All" ))
files<-c("Ts","Tv1.MutFreq","Tv2.MutFreq","Tvs.MutFreq","AllMutFreq" )

for (i in 1:length(files)){
        dt<-mf.files[[i]]
        #coding region only
        dt<-dt[dt$pos>341,]
        
        if (i<=3){
                colnames(dt)[199]<-"Type"
                k<-(i-1)*4+1
                tb$mean[k]<-mean(dt$mean,na.rm=T)
                tb$se[k]<-std.error(dt$mean,na.rm=T)
                k<-k+1
                tb$mean[k]<-mean(dt$mean[dt$Type=="syn"],na.rm=T)
                tb$se[k]<-std.error(dt$mean[dt$Type=="syn"],na.rm=T)
                k<-k+1
                tb$mean[k]<-mean(dt$mean[dt$Type=="nonsyn"],na.rm=T)
                tb$se[k]<-std.error(dt$mean[dt$Type=="nonsyn"],na.rm=T)
                k<-k+1
                tb$mean[k]<-mean(dt$mean[dt$Type=="stop"],na.rm=T)
                tb$se[k]<-std.error(dt$mean[dt$Type=="stop"],na.rm=T)
        }
        
        if (i==4) {k=13
                tb$mean[k]<-mean(dt$mean,na.rm=T)
                tb$se[k]<-std.error(dt$mean,na.rm=T)}
        if (i==5) {k=14
                tb$mean[k]<-mean(dt$mean,na.rm=T)
                tb$se[k]<-std.error(dt$mean,na.rm=T)}
}

tb$type<-as.character(tb$type)

addtb<-data.frame(type=c("tvs.syn","tvs.ns","tvs.stop", "all.syn","all.ns","all.stop"),
                    mean=c(mean(tb$mean[6],tb$mean[10]),mean(tb$mean[7],tb$mean[11]),mean(tb$mean[8],tb$mean[12]),
                           mean(tb$mean[2],tb$mean[6],tb$mean[10]),mean(tb$mean[3],tb$mean[7],tb$mean[11]),mean(tb$mean[4],tb$mean[8],tb$mean[12])),
                    se  =c(mean(tb$se[6],tb$se[10]),mean(tb$se[7],tb$se[11]),mean(tb$se[8],tb$se[12]),
                           mean(tb$se[2],tb$se[6],tb$se[10]),mean(tb$se[3],tb$se[7],tb$se[11]),mean(tb$se[4],tb$se[8],tb$se[12])))

tb2<-rbind(tb, addtb)

write.csv(tb2, "Output1A/MutFreq.filtered/MF.Mean.SE.summary.csv")



#####
#compare the Q30 and Q35 mut freq
#Q30 mut freq
#q30<-read.csv("Output1A/mut.freq.filtered/Summary_TsMutFreq.csv",stringsAsFactors = F)
#r1<-wilcox.test(Ts$mean[Ts$Type=="syn"], q30$mean[q30$Type=="syn"], alternative = "less", paired = FALSE) 
#W = 2011800, p-value < 2.2e-16
#r1[3]  # p.value =  5.283758e-24


######

## Create a summary plot
Ts<-Ts[Ts$pos>=342, ]

summary<-data.frame(Mutation=rep("Transition", times=4), Type=c("All", "Syn", "Nonsyn", "Nonsense"), 
                    Mean=c(mean(Ts$mean, na.rm=T), mean(Ts$mean[Ts$Type=="syn"]), mean(Ts$mean[Ts$Type=="nonsyn"]), mean(Ts$mean[Ts$Type=="stop"])),
                    SE= c(std.error(Ts$mean, na.rm=T), std.error(Ts$mean[Ts$Type=="syn"]),std.error(Ts$mean[Ts$Type=="nonsyn"]), std.error(Ts$mean[Ts$Type=="stop"])))

# Ave. mut freq all
al<-mf.files[[5]]
mean(al$mean) #0.005771959

tvs<-mf.files[[4]]
tvs<-tvs[tvs$pos>=342,]
mean(tvs$mean) #0.0009650088

#Transversion with mean
Tv1<-Tv1.MutFreq
Tv1<-Tv1[Tv1$pos>=342, ]
Tv2<-Tv2.MutFreq
Tv2<-Tv2[Tv2$pos>=342, ]

All<-c(Tv1$mean,Tv2$mean) #0.0004825044
Syn<-c(Tv1$mean[Tv1$Type.tv1=="syn"],Tv2$mean[Tv2$Type.tv2=="syn"])
Nonsyn <- c(Tv1$mean[Tv1$Type.tv1=="nonsyn"],Tv2$mean[Tv2$Type.tv2=="nonsyn"])
Stop <- c(Tv1$mean[Tv1$Type.tv1=="stop"],Tv2$mean[Tv2$Type.tv2=="stop"])

summary2<-data.frame(Mutation=rep("Tranversion", times=4), Type=c("All", "Syn", "Nonsyn", "Nonsense"), 
                     Mean=c(mean(All), mean(Syn), mean(Nonsyn), mean(Stop)),
                     SE= c(std.error(All), std.error(Syn), std.error(Nonsyn),std.error(Stop)))

Summary<-rbind(summary, summary2)
Summary$Type<-factor(Summary$Type, levels=c("All", "Syn","Nonsyn","Nonsense"))

ggplot(Summary,aes(x=Type,y=Mean, ymin=Mean-SE, ymax=Mean+SE, fill=Mutation))+
        geom_bar(position=position_dodge(.9), stat="identity",width=0.8)+
        scale_fill_manual(values=paste0(colors2[c(1,4)],"E6"), labels=c("Transition","Transversion"))+
        geom_errorbar(position=position_dodge(.9), width=.2, color="gray30")+
        theme(axis.title.x=element_blank())+ylab("Mean mutation frequency")+
        theme_linedraw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12, angle=45, hjust=1),axis.title.y = element_text(size=12))+
        geom_vline(xintercept = c(1:3)+0.5,  
                   color = "gray60", size=.4)+
        theme(panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(linetype=2, colour="gray60"),
              panel.grid.minor.y = element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(size=12))

ggsave("Output1A/MutFreq.filtered/MeanMF.byTyep.pdf", width = 5, height = 3.4)
ggsave("Output1A/MutFreq.filtered/MeanMF.byTyep_Presentation.pdf", width = 5, height = 2.8)



### Read depth 
Reads<-data.frame(Seq=names(FilteredOverview2))

for (i in 1:length(FilteredOverview2)){
        dat<-FilteredOverview2[[i]]
        dat<-dat[!is.na(dat$freq.Ts),]
        Reads$ave[i]<-mean(dat$TotalReads, na.rm=T) 
        Reads$median[i]<-median(dat$TotalReads, na.rm = T)
        Reads$max[i]<-max(dat$TotalReads, na.rm=T)
        Reads$min[i]<-min(dat$TotalReads, na.rm=T)
        
}

mean(Reads$ave)
mean(Reads$median) #5282.449
mean(Reads$max) #22005.66

sum(Reads$ave>=5000)

### Structural vs. non-structural genes mut. freq
TS<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F,row.names=1)

st<-TS[TS$pos>341& TS$pos<=2579,]
nonst<-TS[TS$pos<=2579,]

mean(st$mean, na.rm=T) #0.004964672
mean(nonst$mean, na.rm = T)  #0.004872566

r1<-wilcox.test(st$mean,nonst$mean, alternative = "greater", paired = FALSE) 
r1[[3]]  #P=0.1054109 Not Significant
