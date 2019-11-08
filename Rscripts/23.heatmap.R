
#library(colorspace)
library(gplots)

mf<-read.csv("Output_all/mutfreq_3genotypes.csv", stringsAsFactors = F)
ave<-colMeans(mf[2:4],na.rm = T)
row.names(mf)<-mf$merged.pos
mf<-mf[,-1]

#loop doesn't work due to error but can generate figures by running one by one\

genenames<-unique(mf$gene)
par(mar=c(3,4,4,5))
for (i in 6:length(genenames)){
        mname<-genenames[i]
        dat<-mf[mf$gene==mname,]
        dat<-dat[,2:4]
        dat.m<-data.matrix(dat)
        l<-nrow(dat.m)
        pdf(paste0("Output_all/Heatmap.",mname,".pdf"), width = l/20, height = 3)
        heatmap.2(t(dat.m), margin=c(2,6),col =colorRampPalette(c("white","red"))(50), 
                  cexCol=.3,cexRow = 1.3, trace="none",density.info="none",Rowv=FALSE,
                  key=F,na.color="grey",Colv = FALSE,
                  lmat=rbind(c(3,4),c(2,1)),lhei=c(.2,4),lwid = c(0.2,4))
        dev.off()
        dev.off()
        
        
}

margin=c(2,6),
heatmap.2(Core, col =colorRampPalette(c("white","red"))(50), 
          cexCol=1,cexRow=.3, trace="none",density.info="none",Rowv=FALSE,
          keysize=1, key.xlab="Mut freq", key.title="",na.color="grey",Colv = FALSE)
dev.off()


#############
# find the most different
