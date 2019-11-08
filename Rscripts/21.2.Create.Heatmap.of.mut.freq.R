# Create heatmap of mutation frequency by gene to compare across genotypes 


#library(colorspace)
library(gplots)

mf<-read.csv("Output_all/mutfreq_3genotypes.csv", stringsAsFactors = F)
ave<-colMeans(mf[2:4],na.rm = T)
row.names(mf)<-mf$merged.pos
mf<-mf[,-1]

#loop doesn't work due to error but can generate figures by running one by one\

genenames<-unique(mf$gene)
par(mar=c(3,1,1,4))
for (i in 2:length(genenames)){
        mname<-genenames[i]
        dat<-mf[mf$gene==mname,]
        dat<-dat[,2:4]
        dat.m<-data.matrix(dat)
        l<-nrow(dat.m)
        pdf(paste0("Output_all/Heatmap.",mname,".pdf"), width = l/20, height = 3)
        par(mar=c(1,1,1,4))
        heatmap.2(t(dat.m), margin=c(2,6),col =colorRampPalette(c("white","red"))(50), 
                  cexCol=.3,cexRow = 1.3, trace="none",density.info="none",Rowv=FALSE,
                  key=F,na.color="grey",Colv = FALSE, lhei=c(.2, 4), lwid=c(0.2,4))
        dev.off()
}


