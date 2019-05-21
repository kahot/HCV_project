library(stringr)
library(R.utils)

setwd("~/programs/HCV/")
#other genotypes

#cm1B<-readLines("Data/template/D7xxxx1_1B.sh")
#cm2.1B<-readLines("Data/template/D7xxxx2_1B.sh")
#cm3.1B<-readLines("Data/template/D7xxxx3_1B.sh")
cm1<-readLines("Data/template/D7xxxx1_3A.sh")
#cm2<-readLines("Data/template/D7xxxx2_3A.sh")
cm3<-readLines("Data/template/D7xxxx3_3A.sh")





genotypes<-c("1B","2A","2B","3A")

for (geno in 1:length(genotypes)){
        dir.create(paste0("~/programs/HCV/Bashscripts/bash1_",genotypes[geno]))
        dir.create(paste0("~/programs/HCV/Bashscripts/bash2_",genotypes[geno]))
        dir.create(paste0("~/programs/HCV/Bashscripts/bash3_",genotypes[geno]))
}


for (geno in 1:length(genotypes)){
        fq<-list.files(paste0("Data/fastq_",genotypes[geno],"/"), pattern="fastq") 
        
        #create vector of odd numbers:
        n<-seq(1, by = 2, len = (length(fq)/2))
        fq2<-fq[n]
        for (i in 1:length(fq2)){
                #choose the paired reads fastq files
                fa1<-fq2[i]
                fa2<-gsub(pattern="R1",replace="R2",x=fa1)
                fname<-substr(fa1,start=1,stop=7)
                new<-gsub(pattern="D75000-HCV_S11_L001_R1_001.fastq", replace=paste0(fa1),x=cm1)
                new<-gsub(pattern="D75000-HCV_S11_L001_R2_001.fastq", replace=paste0(fa2),x=new)
                new<-gsub(pattern="unzipped", replace=paste0("fastq_",genotypes[geno]), x=new)
                new<-gsub(pattern="REFERENCE", replace=paste0("HCV",genotypes[geno]), x=new)
                new<-gsub(pattern="D75000",replace=paste0(fname),x=new)
                writeLines(new, con=paste0("Bashscripts/bash1_",genotypes[geno],"/",fname,".sh"))
                
                #new2<-gsub(pattern="D75000",replace=paste0(fname),x=cm2.1B)
                #writeLines(new2, con=paste0("Bashscripts/bash2_",genotypes[geno],"/",fname,"2.sh"))
                
                new3<-gsub(pattern="D75000",replace=paste0(fname),x=cm3)
                writeLines(new3, con=paste0("Bashscripts/bash3_",genotypes[geno],"/",fname,"3.sh"))
        }
}


### MKDIR ###
cm<-readLines("Data/template/mkdir.sh")

f<-list.files("~/programs/HCV/Output2/Consensus/1B/", pattern="fasta") 

#create vector of odd numbers:
i=1
dir.create("Bashscripts/bash_mkdir/")
for (i in 1:length(f)){
        fname<-substr(f[i],start=1,stop=7)
        new<-gsub(pattern="D75000",replace=paste0(fname),x=cm)
        writeLines(new, con=paste0("Bashscripts/bash_mkdir/",fname,".sh"))
        }
