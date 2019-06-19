## prepare summary data for mutation frequenceis across all samples

library(plotrix)
source("Rscripts/CreateEmptyData.R")
CreateEmptyData()
for (i in 1:length(Overview_summary)){
        dat<-Overview_summary[[i]]
        filename<-names(Overview_summary)[i]

        for (typeofsite in c("syn", "nonsyn")){
                for (wtnt in c("a", "t", "c", "g")){
                        mutrate<- dat$freq.Ts[dat$Type==typeofsite & dat$MajNt==wtnt]
                        m[typeofsite,wtnt]<-mean(mutrate[!is.na(mutrate)])
                        se[typeofsite,wtnt]<-std.error(mutrate[!is.na(mutrate)])
                        
                        m_NonCpG<-dat$freq.Ts[dat$Type==typeofsite & dat$MajNt==wtnt & dat$makesCpG==0]
                        m_nonCpG[typeofsite,wtnt]<-mean(m_NonCpG[!is.na(m_NonCpG)])
                        se_nonCpG[typeofsite,wtnt]<-std.error(m_NonCpG[!is.na(m_NonCpG)])
                        
                        mu_CpG<-dat$freq.Ts[dat$Type==typeofsite & dat$MajNt==wtnt & dat$makesCpG==1]
                        m_CpG[typeofsite,wtnt]<-mean(mu_CpG[!is.na(mu_CpG)])
                        se_CpG[typeofsite,wtnt]<-std.error(mu_CpG[!is.na(mu_CpG)])
                        
                        #TV1
                        mutrate1<-dat$freq.transv1[dat$Type.tv1==typeofsite & dat$MajNt==wtnt]
                        m1[typeofsite,wtnt]<-mean(mutrate1[!is.na(mutrate1)])
                        se1[typeofsite,wtnt]<-std.error(mutrate1[!is.na(mutrate1)])
                        
                        m1_NonCpG<-dat$freq.transv1[dat$Type.tv1==typeofsite & dat$MajNt==wtnt & dat$makesCpG.tv1==0]
                        m1_nonCpG[typeofsite,wtnt]<-mean(m1_NonCpG[!is.na(m1_NonCpG)])
                        se1_nonCpG[typeofsite,wtnt]<-std.error(m1_NonCpG[!is.na(m1_NonCpG)])
                        
                        mu1_CpG<-dat$freq.transv1[dat$Type.tv1==typeofsite & dat$MajNt==wtnt & dat$makesCpG.tv1==1]
                        m1_CpG[typeofsite,wtnt]<-mean(mu1_CpG[!is.na(mu1_CpG)])
                        se1_CpG[typeofsite,wtnt]<-std.error(mu1_CpG[!is.na(mu1_CpG)])
                        
                        #tv2
                        mutrate2<-dat$freq.transv2[dat$Type.tv2==typeofsite & dat$MajNt==wtnt]
                        m2[typeofsite,wtnt]<-mean(mutrate2[!is.na(mutrate2)])
                        se2[typeofsite,wtnt]<-std.error(mutrate2[!is.na(mutrate2)])
                        
                        m2_NonCpG<-dat$freq.transv2[dat$Type.tv2==typeofsite & dat$MajNt==wtnt & dat$makesCpG.tv2==0]
                        m2_nonCpG[typeofsite,wtnt]<-mean(m2_NonCpG[!is.na(m2_NonCpG)])
                        se2_nonCpG[typeofsite,wtnt]<-std.error(m2_NonCpG[!is.na(m2_NonCpG)])
                        
                        mu2_CpG<-dat$freq.transv2[dat$Type.tv2==typeofsite & dat$MajNt==wtnt & dat$makesCpG.tv2==1]
                        m2_CpG[typeofsite,wtnt]<-mean(mu2_CpG[!is.na(mu2_CpG)])
                        se2_CpG[typeofsite,wtnt]<-std.error(mu2_CpG[!is.na(mu2_CpG)])

                        vectorname<<-paste0(typeofsite,"_",wtnt)
                        assign(vectorname, mutrate)
                        vectorname.tv1<-paste0(typeofsite,"_",wtnt,"_tv1")
                        vectorname.tv2<-paste0(typeofsite,"_",wtnt,"_tv2")
                        assign (vectorname.tv1, mutrate1)
                        assign (vectorname.tv2, mutrate2)
                        
                        vname1<<-paste0(typeofsite,"_",wtnt,"_noncpg")
                        assign(vname1, m_NonCpG)
                        vectorname.tv1<-paste0(typeofsite,"_",wtnt,"_tv1_noncpg")
                        vectorname.tv2<-paste0(typeofsite,"_",wtnt,"_tv2_noncpg")
                        assign (vectorname.tv1, m1_NonCpG)
                        assign (vectorname.tv2, m2_NonCpG)
                        
                        vname2<<-paste0(typeofsite,"_",wtnt,"_cpg")
                        assign(vname2, mu_CpG)
                        vectorname.tv1<-paste0(typeofsite,"_",wtnt,"_tv1_cpg")
                        vectorname.tv2<-paste0(typeofsite,"_",wtnt,"_tv2_cpg")
                        assign (vectorname.tv1, mu1_CpG)
                        assign (vectorname.tv2, mu2_CpG)
                      
                        }
        }
        
        A_syn<<-c(A_syn,syn_a)
        T_syn<<-c(T_syn,syn_t)
        C_syn<<-c(C_syn,syn_c)
        G_syn<<-c(G_syn,syn_g)
        A_nonsyn<<-c(A_nonsyn,nonsyn_a)
        T_nonsyn<<-c(T_nonsyn,nonsyn_t)
        C_nonsyn<<-c(C_nonsyn,nonsyn_c)
        G_nonsyn<<-c(G_nonsyn,nonsyn_g)
        
        A_syn_noncpg<<-c(A_syn_noncpg,syn_a_noncpg)
        T_syn_noncpg<<-c(T_syn_noncpg,syn_t_noncpg)
        C_syn_noncpg<<-c(C_syn_noncpg,syn_c_noncpg)
        G_syn_noncpg<<-c(G_syn_noncpg,syn_g_noncpg)
        A_nonsyn_noncpg<<-c(A_nonsyn_noncpg,nonsyn_a_noncpg)
        T_nonsyn_noncpg<<-c(T_nonsyn_noncpg,nonsyn_t_noncpg)
        C_nonsyn_noncpg<<-c(C_nonsyn_noncpg,nonsyn_c_noncpg)
        G_nonsyn_noncpg<<-c(G_nonsyn_noncpg,nonsyn_g_noncpg)
        A_syn_noncpg<<-c(A_syn_noncpg,syn_a_noncpg)
        T_syn_noncpg<<-c(T_syn_noncpg,syn_t_noncpg)
        C_syn_noncpg<<-c(C_syn_noncpg,syn_c_noncpg)
        G_syn_noncpg<<-c(G_syn_noncpg,syn_g_noncpg)
        A_nonsyn_noncpg<<-c(A_nonsyn_noncpg,nonsyn_a_noncpg)
        T_nonsyn_noncpg<<-c(T_nonsyn_noncpg,nonsyn_t_noncpg)
        C_nonsyn_noncpg<<-c(C_nonsyn_noncpg,nonsyn_c_noncpg)
        G_nonsyn_noncpg<<-c(G_nonsyn_noncpg,nonsyn_g_noncpg)
        
        A_tv1_syn_cpg<<-c(A_tv1_syn_cpg,syn_a_tv1_cpg)
        A_tv1_syn_noncpg<<-c(A_tv1_syn_noncpg,syn_a_tv1_noncpg)
        A_tv1_nonsyn_cpg<<-c(A_tv1_nonsyn_cpg,nonsyn_a_tv1_cpg)
        A_tv1_nonsyn_noncpg<<-c(A_tv1_nonsyn_noncpg,nonsyn_a_tv1_noncpg)
        
        C_tv1_syn_noncpg<<-c(C_tv1_syn_noncpg,syn_c_tv1_noncpg)
        C_tv1_nonsyn_noncpg<<-c(C_tv1_nonsyn_noncpg,nonsyn_c_tv1_noncpg)
        
        G_tv1_syn_cpg<<-c(G_tv1_syn_cpg,syn_g_tv1_cpg)
        G_tv1_syn_noncpg<<-c(G_tv1_syn_noncpg,syn_g_tv1_noncpg)
        G_tv1_nonsyn_cpg<<-c(G_tv1_nonsyn_cpg,nonsyn_g_tv1_cpg)
        G_tv1_nonsyn_noncpg<<-c(G_tv1_nonsyn_noncpg,nonsyn_g_tv1_noncpg)
        
        T_tv1_syn_noncpg<<-c(T_tv1_syn_noncpg,syn_t_tv1_noncpg)
        T_tv1_nonsyn_noncpg<<-c(T_tv1_nonsyn_noncpg,nonsyn_t_tv1_noncpg)
        
        
        A_tv2_syn<<-c(A_tv2_syn,syn_a_tv2)
        T_tv2_syn<<-c(T_tv2_syn,syn_t_tv2)
        C_tv2_syn<<-c(C_tv2_syn,syn_c_tv2)
        G_tv2_syn<<-c(G_tv2_syn,syn_g_tv2)
        A_tv2_nonsyn<<-c(A_tv2_nonsyn,nonsyn_a_tv2)
        T_tv2_nonsyn<<-c(T_tv2_nonsyn,nonsyn_t_tv2)
        C_tv2_nonsyn<<-c(C_tv2_nonsyn,nonsyn_c_tv2)
        G_tv2_nonsyn<<-c(G_tv2_nonsyn,nonsyn_g_tv2)
        
        A_tv2_syn_noncpg<<-c(A_tv2_syn_noncpg,syn_a_tv2_noncpg)
        A_tv2_nonsyn_noncpg<<-c(A_tv2_nonsyn_noncpg,nonsyn_a_tv2_noncpg)
        
        C_tv2_syn_cpg<<-c(C_tv2_syn_cpg,syn_c_tv2_cpg)
        C_tv2_syn_noncpg<<-c(C_tv2_syn_noncpg,syn_c_tv2_noncpg)
        C_tv2_nonsyn_cpg<<-c(C_tv2_nonsyn_cpg,nonsyn_c_tv2_cpg)
        C_tv2_nonsyn_noncpg<<-c(C_tv2_nonsyn_noncpg,nonsyn_c_tv2_noncpg)
        
        G_tv2_syn_noncpg<<-c(G_tv2_syn_noncpg,syn_g_tv2_noncpg)
        G_tv2_nonsyn_noncpg<<-c(G_tv2_nonsyn_noncpg,nonsyn_g_tv2_noncpg)
        
        T_tv2_syn_cpg<<-c(T_tv2_syn_cpg,syn_t_tv2_cpg)
        T_tv2_syn_noncpg<<-c(T_tv2_syn_noncpg,syn_t_tv2_noncpg)
        T_tv2_nonsyn_cpg<<-c(T_tv2_nonsyn_cpg,nonsyn_t_tv2_cpg)
        T_tv2_nonsyn_noncpg<<-c(T_tv2_nonsyn_noncpg,nonsyn_t_tv2_noncpg)
        
        mut[[i]]<-m
        mut.CpG[[i]]<-m_CpG
        mut.nonCpG[[i]]<-m_nonCpG
        names(mut)[i]<-filename
        names(mut.CpG)[i]<-filename
        names(mut.nonCpG)[i]<-filename
        
        mut1[[i]]<-m1
        mut1.CpG[[i]]<-m1_CpG
        mut1.nonCpG[[i]]<-m1_nonCpG
        names(mut1)[i]<-filename
        names(mut1.CpG)[i]<-filename
        names(mut1.nonCpG)[i]<-filename
        
        mut2[[i]]<-m2
        mut2.CpG[[i]]<-m2_CpG
        mut2.nonCpG[[i]]<-m2_nonCpG
        names(mut2)[i]<-filename
        names(mut2.CpG)[i]<-filename
        names(mut2.nonCpG)[i]<-filename
        
        SE[[i]]<-se
        SE.CpG[[i]]<-se_CpG
        SE.nonCpG[[i]]<-se_nonCpG
        names(SE)[i]<-paste0(filename,".se")
        names(SE.CpG)[i]<-paste0(filename,".se")
        names(SE.nonCpG)[i]<-paste0(filename,".se")
        
        SE1[[i]]<-se1
        SE1.CpG[[i]]<-se1_CpG
        SE1.nonCpG[[i]]<-se1_nonCpG
        names(SE1)[i]<-paste0(filename,".se")
        names(SE1.CpG)[i]<-paste0(filename,".se")
        names(SE1.nonCpG)[i]<-paste0(filename,".se")
        
        SE2[[i]]<-se2
        SE2.CpG[[i]]<-se2_CpG
        SE2.nonCpG[[i]]<-se2_nonCpG
        names(SE2)[i]<-filename
        names(SE2.CpG)[i]<-filename
        names(SE2.nonCpG)[i]<-filename


}


Summary<-data.frame(matrix(ncol = 0, nrow = 4),row.names =c("A","T",'C','G') )
dir="./Output1A/MutFreq/Maj/"

for (i in c("mut","mut1","mut2")){
        if (i=="mut") { k=1
        n<-"" 
        error<-SE
        fname<-"Transition"}
        if (i=="mut1") {k=2
        n=1
        error<-SE1
        fname<-"Transversion_1"}
        if (i=="mut2") {k=3; n=2
        error<-SE2
        fname<-"Transversion_2"}
        
        syn<-do.call(rbind, lapply(get(paste0(i)), function(x) x[1,]))
        nonsyn<-do.call(rbind, lapply(get(paste0(i)), function(x) x[2,]))
        syn.se<-do.call(rbind, lapply(error, function(x) x[1,]))
        nonsyn.se<-do.call(rbind, lapply(error, function(x) x[2,]))
        SYN<-rbind(syn, syn.se)
        NonSYN<-rbind(nonsyn, nonsyn.se)
        
        Summary$ave.syn<-colMeans(syn)
        Summary$ave.nonsyn<-colMeans(nonsyn)
        Summary$ave.syn.se<-colMeans(syn.se)
        Summary$ave.nonsyn.se<-colMeans(nonsyn.se)
        
        CpG<-do.call(rbind, lapply(get(paste0(i,".CpG")), function(x) x[1,]))
        CpG.ns<-do.call(rbind, lapply(get(paste0(i,".CpG")), function(x) x[2,]))
        CpG.se<-do.call(rbind, lapply(get(paste0("SE",n,".CpG")), function(x) x[1,]))
        CpG.ns.se<-do.call(rbind, lapply(get(paste0("SE",n,".CpG")), function(x) x[2,]))
        CPG.s<-rbind(CpG,CpG.se)
        CPG.ns<-rbind(CpG.ns,CpG.ns.se)
        
        Summary$CpG.ave<-colMeans(CpG)
        Summary$CpG.ns.ave<-colMeans(CpG.ns)
        Summary$CpG.ave.se<-colMeans(CpG.se)
        Summary$CpG.ns.ave.se<-colMeans(CpG.se)
        
        nonCpG<-do.call(rbind, lapply(get(paste0(i,".nonCpG")), function(x) x[1,]))
        nonCpG.se<-do.call(rbind, lapply(get(paste0("SE",n,".nonCpG")), function(x) x[1,]))
        nonCpG.ns<-do.call(rbind, lapply(get(paste0(i,".nonCpG")), function(x) x[2,]))
        nonCpG.ns.se<-do.call(rbind, lapply(get(paste0("SE",n,".nonCpG")), function(x) x[2,]))
        nonCPG.s<-rbind(nonCpG,nonCpG.se)
        nonCPG.ns<-rbind(nonCpG.ns,nonCpG.ns.se)
        
        Summary$nonCpG.ave<-colMeans(nonCpG)
        Summary$nonCpG.ns.ave<-colMeans(nonCpG.ns)
        Summary$nonCpG.ave.se<-colMeans(nonCpG.se)
        Summary$nonCpG.ns.ave.se<-colMeans(nonCpG.se)
        
        write.csv(Summary,paste0(dir,"Summary_",fname,".csv"))
        write.csv(SYN, paste0(dir,"Synonymous_",fname,".csv"))
        write.csv(NonSYN,paste0(dir,"Nonsynonymous_",fname,".csv"))
        write.csv(CPG.s,paste0(dir, "CpG_Synonymous_",fname,".csv"))
        write.csv(CPG.ns,paste0(dir, "CpG_Nonsynonymous_",fname,".csv"))
        write.csv(nonCPG.s,paste0(dir,"nonCpG_Synonymous_", fname,".csv"))
        write.csv(nonCPG.ns,paste0(dir,"nonCpG_Nonsynonymous_", fname,".csv"))
}


