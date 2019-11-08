#Step 1. Create bash executable files to process (filter, trim and map) the fastq fies
library(stringr)
library(R.utils)


#create the bash files to run bbmap and bwa
# read the template command text file:
cmmd<-readLines("Data/template/D7xxxx1_Q35.sh")
cmmd2<-readLines("Data/template/D7xxxx2.sh")
cmmd3<-readLines("Data/template/D7xxxx3.sh")

cm<-readLines("Data/template/D7xxxx1_Q35_2.sh")

#choose the fastq files to be prrocessed
fq<-list.files("Data/unzipped/",pattern="fastq") 

#create vector of odd numbers:
n<-seq(1, by = 2, len = (length(fq)/2))
fq2<-fq[n]
for (i in 1:length(fq2)){
        #choose the paired reads fastq files
        fa1<-fq2[i]
        fa2<-gsub(pattern="R1",replace="R2",x=fa1)
        fname<-substr(fa1,start=1,stop=7)
        new<-gsub(pattern="D75000-HCV_S11_L001_R1_001.fastq", replace=paste0(fa1),x=cm)
        new<-gsub(pattern="D75000-HCV_S11_L001_R2_001.fastq", replace=paste0(fa2),x=new)
        new<-gsub(pattern="D75000",replace=paste0(fname),x=new)
        writeLines(new, con=paste0("Bashscripts/Bash1_Q35/",fname,".sh"))
        
        #new2<-gsub(pattern="D75000",replace=paste0(fname),x=cmmd2)
        #writeLines(new2, con=paste0("Bashscripts/bash2_Q35/",fname,"2.sh"))
        
        new3<-gsub(pattern="D75000",replace=paste0(fname),x=cmmd3)
        writeLines(new3, con=paste0("Bashscripts/bash3_Q35/",fname,"3.sh"))
        
        
}


########
#other genotypes

cmmd<-readLines("Data/template/D7xxxx1_genotypes.sh")
cmmd2<-readLines("Data/template/D7xxxx2.sh")
cmmd3<-readLines("Data/template/D7xxxx3.sh")

genotypes<-c("1A","2A","2B","3A")
geno=1
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
                new<-gsub(pattern="D75335R-HCV_S11_L001_R1_001.fastq", replace=paste0(fa1),x=cmmd)
                new<-gsub(pattern="D75335R-HCV_S11_L001_R2_001.fastq", replace=paste0(fa2),x=new)
                new<-gsub(pattern="REFERENCE", replace=paste0("HCV",genotypes[geno]), x=new)
                new<-gsub(pattern="REFERENCE", replace=paste0("HCV",genotypes[geno]), x=new)
                new<-gsub(pattern="D75335",replace=paste0(fname),x=new)
                writeLines(new, con=paste0("Bashscripts/bash1/",fname,".sh"))
                
                new2<-gsub(pattern="D75000",replace=paste0(fname),x=cmmd2)
                writeLines(new2, con=paste0("Bashscripts/bash2/",fname,"2.sh"))
                
                new3<-gsub(pattern="D75000",replace=paste0(fname),x=cmmd3)
                writeLines(new3, con=paste0("Bashscripts/bash3/",fname,"3.sh"))
                
                
                
        }
        
        

#### Bash1.2 script is for the sequences I finished up to clean.fq. Starting with clumping ~


cmmd1<-readLines("Data/template/D7xxxx1.2.sh")
cmmd3<-readLines("Data/template/D7xxxx3.sh")

f<-list.files("Output/clean.fq/",pattern=".fq$")

for (i in 1:length(f)){
        fname<-substr(f[i],start=1,stop=7)
        new<-gsub(pattern="D75000",replace=paste0(fname),x=cmmd1)
        writeLines(new, con=paste0("Bashscripts/bash1.2/",fname,".sh"))
        
        new3<-gsub(pattern="D75000",replace=paste0(fname),x=cmmd3)
        writeLines(new3, con=paste0("Bashscripts/bash3.2/",fname,"3.sh"))
        
        
        
}


#### Bash1.3 script is for the sequences with its own consensus made but rerun with adapter removal


cmmd1.3<-readLines("Data/template/D7xxxx1.3.sh")
cmmd3<-readLines("Data/template/D7xxxx3.sh")

#create vector of odd numbers:
fq<-list.files("Data/unzipped/",pattern="fastq") 
n<-seq(1, by = 2, len = (length(fq)/2))
fq2<-fq[n]
for (i in 1:length(fq2)){
        fa1<-fq2[i]
        fa2<-gsub(pattern="R1",replace="R2",x=fa1)
        fname<-substr(fa1,start=1,stop=7)
        new<-gsub(pattern="D75000-HCV_S11_L001_R1_001.fastq", replace=paste0(fa1),x=cmmd1.3)
        new<-gsub(pattern="D75000-HCV_S11_L001_R2_001.fastq", replace=paste0(fa2),x=new)

        new<-gsub(pattern="D75000",replace=paste0(fname),x=new)
        writeLines(new, con=paste0("Bashscripts/bash1.3/",fname,".sh"))
        
        new3<-gsub(pattern="D75000",replace=paste0(fname),x=cmmd3)
        writeLines(new3, con=paste0("Bashscripts/bash3.3/",fname,"3.sh"))
     
}

###11.26.18 remapping to create better consensus sequences
cm1<-readLines("Data/template/D7xxx2_map.sh")
f<-list.files("Output/Consen/",pattern=".fasta")

for (i in 1:length(f)){
        fname<-substr(f[i],start=1,stop=7)
        new<-gsub(pattern="D75000",replace=paste0(fname),x=cm1)
        writeLines(new, con=paste0("Bashscripts/",fname,".sh"))
        
}

cm2<-readLines("Data/template/D7xxxx3_map.sh")

f<-list.files("Output/clumped.fq/",pattern=".fq$")

for (i in 1:length(f)){
        fname<-substr(f[i],start=1,stop=7)
        new<-gsub(pattern="D75000",replace=paste0(fname),x=cm2)
        writeLines(new, con=paste0("Bashscripts/",fname,".sh"))
}



##############

remove<-list.files("Output/Bash/",pattern=".fq", recursive = T)
# in bash 
find . -type f -name '*.fq' -delete




#######################################################


#####Section1 -Skip to Section2 once all fiels are unzipped and moved to fastq ######
#Sec.1
#Find the sample ID fos sepcific strains
samples<-read.csv("Data/CoInfectionStudy_HCV_Genomes.csv")

str(samples)
genotypes<-levels(samples$Sequence.Genotype)
#"1A" "1B" "2A" "2B" "2Q" "3A" "3I" "4D" "4K" "4Q" "4R" "4T" "4V" "6L"

#select samples of each genotype
gtypes<-c()
for (i in 1:length(genotypes)){
        dataname<-paste0("s",genotypes[i])
        gtypes[i]<-dataname
        data<-subset(samples,Sequence.Genotype==genotypes[i])
        assign(dataname,data)
}


#list all compressed fastq files
gz<-list.files("/Volumes/Kaho_Backup/HCV/Data/fastq/",pattern="gz") #916 files

#select the files contatinig 'reads' (R1/R2 fastq files, not the I1/I2 files)
gz2<-sapply("R1|R2", function(y){gz[grepl(pattern =y, x=gz)]}) #551 files


####################################
#### process genotype 1A only #######
#select genotype 1A files
gz1A<-sapply(s1A$Enumber, function(y){gz2[grepl(pattern =y, x=gz2)]})
#eliminate the empty elements
gz1A<-gz1A[lapply(gz1A,length)>0] #196

#move the 1A file 
dir.create("/Volumes/Kaho_Backup/HCV/Data/1A")
for (i in 1:length(gz1A)){
        file.copy(paste0("/Volumes/Kaho_Backup/HCV/Data/fastq/",gz1A[[i]][1]), 
                  paste0("/Volumes/Kaho_Backup/HCV/Data/1A/",gz1A[[i]][1]))
        file.copy(paste0("/Volumes/Kaho_Backup/HCV/Data/fastq/",gz1A[[i]][2]), 
                  paste0("/Volumes/Kaho_Backup/HCV/Data/1A/",gz1A[[i]][2]))
}


#unzip 1A genotype files:
for (i in 1:length(gz1A)){
        gunzip(paste0("/Volumes/Kaho_Backup/HCV/Data/1A/",gz1A[[i]][1]))
        gunzip(paste0("/Volumes/Kaho_Backup/HCV/Data/1A/",gz1A[[i]][2]))
}

#move the fastq files to 'unzipped/' directory
#fastqfiles<-list.files("Data/fastq/", pattern="fastq$")
#for (i in 1:length(fastqfiles)) {
#        file.rename(paste0("Data/fastq/",fastqfiles[i]), paste0("Data/unzipped/", fastqfiles[i]))
#}



#########process all genotypes#######
gz<-list.files("Data/fastq/",pattern="gz")
gz2<-sapply("R1|R2", function(y){gz[grepl(pattern =y, x=gz)]})

#### move each genotype fastq files into a separate separately ###
# no 2Q,3I,4D,4K,4Q,4R,4T,4V file exists (skip k=5,7:13)
for (k in 1:c(1:4,6)){
        data<-get(gtypes[k])
        gzsubset<-sapply(data$Enumber, function(y){gz2[grepl(pattern =y, x=gz2)]})
        gzsubset<-gzsubset[lapply(gzsubset,length)>0]
        
        dir.create(paste0("Data/fastq_",genotypes[k]))
        for (i in 1:length(gzsubset)){
                gunzip(paste0("Data/fastq/",gzsubset[[i]][1]),remove=F)
                gunzip(paste0("Data/fastq/",gzsubset[[i]][2]),remove=F)
        }
        
        fastqfiles<-list.files("Data/fastq/", pattern=".fastq$")
        for (i in 1:length(fastqfiles)) {
                file.copy(paste0("Data/fastq/",fastqfiles[i]), 
                          paste0("Data/fastq_",genotypes[k],"/", fastqfiles[i]))
                file.rename(paste0("Data/fastq/",fastqfiles[i]), paste0("Data/unzipped/",fastqfiles[i]))
        }
        
} 




##################################################
##################################################
