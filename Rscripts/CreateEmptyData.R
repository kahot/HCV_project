CreateEmptyData<-function(){
        dfnames<-c('m','m1','m2','se','se1','se2')
        dfnames1<-sapply(dfnames,function(x) paste0(x,"_CpG"))
        dfnames2<-sapply(dfnames,function(x) paste0(x,"_nonCpG"))
        dfnames<-c(dfnames,dfnames1,dfnames2)
        df<-sapply(dfnames,function(x) x<-data.frame())
        names(df)<-dfnames
        list2env(df,envir=.GlobalEnv)
        
        lnames<-c('mut','mut1','mut2','SE','SE1','SE2')
        lnames1<-sapply(lnames,function(x) paste0(x,".CpG"))
        lnames2<-sapply(lnames,function(x) paste0(x,".nonCpG"))
        lnames<-c(lnames,lnames1,lnames2)
        
        lists<-sapply(lnames,function(x) x<-list())
        unlist(lists)
        names(lists)<-lnames
        list2env(lists,envir=.GlobalEnv)
        
        vnames<-c("A","T","C","G")
        for (k in vnames){
                vn1<-paste0(vnames,"_syn")
                vn2<-paste0(vnames,"_nonsyn")
                vn3<-paste0(vnames, "_tv1_syn")
                vn4<-paste0(vnames,"_tv2_syn")
                vn5<-paste0(vnames, "_tv1_nonsyn")
                vn6<-paste0(vnames,"_tv2_nonsyn")
                vn7<-paste0(vnames, "_syn_cpg")
                vn8<-paste0(vnames,"_syn_noncpg")
                vn9<-paste0(vnames, "_nonsyn_cpg")
                vn10<-paste0(vnames,"_nonsyn_noncpg")
                vn11<-paste0(vnames, "_tv1_syn_cpg")
                vn12<-paste0(vnames,"_tv1_syn_noncpg")
                vn13<-paste0(vnames, "_tv1_nonsyn_cpg")
                vn14<-paste0(vnames,"_tv1_nonsyn_noncpg")
                vn15<-paste0(vnames, "_tv2_syn_cpg")
                vn16<-paste0(vnames,"_tv2_syn_noncpg")
                vn17<-paste0(vnames, "_tv2_nonsyn_cpg")
                vn18<-paste0(vnames,"_tv2_nonsyn_noncpg")
                vn19<-paste0(vnames,"_stop")
                vn20<-paste0(vnames,"_tv1_stop")
                vn21<-paste0(vnames,"_tv2_stop")
        }
       
        v2names<-c(noquote(paste0("vn",seq(1:21))))
        for (i in (1:length(v2names))){
                v<-get(v2names[i])
                vct<-sapply(v,function(x) x<-c())
                list2env(vct,envir=.GlobalEnv)}
        
}
