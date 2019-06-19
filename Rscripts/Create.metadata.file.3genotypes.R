source("Rscripts/baseRscript.R")

## Create a metadata file for 3 genotypes with amino acid info to all bases (some in overview files are missing)
geno<-c("1A","1B","3A")
M1<-read.csv("Output_all/Overview1A/D75002-_overview4.csv", stringsAsFactors = F)
M2<-read.csv("Output_all/Overview1B/D75046-_overview4.csv", stringsAsFactors = F)
M3<-read.csv("Output_all/Overview3A/D75003-_overview4.csv", stringsAsFactors = F)


for (f in 1:3){
        m<-get(paste0("M",f))
        muttypes<-m[,c(2,3,4,23,26,30,35,36,39,42)]
        muttypes1<-muttypes[muttypes$nuc!="-",]
        for (k in 3:nrow(muttypes1)){
                if (k%%3==0){
                        muttypes1$WTAA[k] = seqinr::translate(muttypes1$nuc[c(k,k+1,k+2)])
                        muttypes1$MUTAA[k] = seqinr::translate(c(transition(muttypes1$nuc[k]),muttypes1$nuc[c(k+1,k+2)])) }
                if (k%%3==1){
                        muttypes1$WTAA[k] = seqinr::translate(muttypes1$nuc[c(k-1,k,k+1)])
                        muttypes1$MUTAA[k] = seqinr::translate(c(muttypes1$nuc[c(k-1)],transition(muttypes1$nuc[k]),muttypes1$nuc[c(k+1)]))  }
                if (k%%3==2){
                        muttypes1$WTAA[k] = seqinr::translate(muttypes1$nuc[c(k-2,k-1,k)])
                        muttypes1$MUTAA[k] = seqinr::translate(c(muttypes1$nuc[c(k-2,k-1)],transition(muttypes1$nuc[k]))) }
        }
        
        muttypes2<-merge(muttypes[,c(1,3)], muttypes1[,c(1,2,4:10)], by="merged.pos",all.x=T)
        columnnames1<-paste0(colnames(muttypes2),".",geno[f])
        colnames(muttypes2)[2:10]<-columnnames1[2:10]
        filename<-paste0("muttypes",geno[f])
        assign(filename,muttypes2)
        
}
M<-merge(muttypes3A, muttypes1A, by="merged.pos")
M<-merge(M, muttypes1B, by="merged.pos")


genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}

genetable<-data.frame("merged.pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector
M<-merge(genetable,M, by="merged.pos")


write.csv(M,"Output_all/merged.metadata.csv")

