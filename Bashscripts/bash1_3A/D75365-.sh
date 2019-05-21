#!/bin/bash


mkdir Output3A/Bash/D75365-
 
#0 adapter trimming
bbduk.sh in1=~/programs/HCV/Data/fastq_3A/D75365-HCV_S24_L001_R1_001.fastq in2=~/programs/HCV/Data/fastq_3A/D75365-HCV_S24_L001_R2_001.fastq  out=~/programs/HCV/Output3A/Bash/D75365-/D75365-_adp.trimmed.fastq ref=/Users/kahotisthammer/programs/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=~/programs/HCV/Output3A/Bash/D75365-/stats35_0.txt

#1 Trim reads at both ends at Score<15
bbduk.sh in=~/programs/HCV/Output3A/Bash/D75365-/D75365-_adp.trimmed.fastq out=~/programs/HCV/Output3A/Bash/D75365-/D75365-_trimmed.Q35.fastq qtrim=rl trimq=35 stats=~/programs/HCV/Output3A/Bash/D75365-/stats35_1.txt

#2. Kmer filtering
bbduk.sh in=~/programs/HCV/Output3A/Bash/D75365-/D75365-_trimmed.Q35.fastq out=~/programs/HCV/Output3A/Bash/D75365-/D75365-_unmatched.Q35.fq ref=/Users/kahotisthammer/programs/bbmap/resources/phix174_ill.ref.fa k=31 hdist=1 stats=~/programs/HCV/Output3A/Bash/D75365-/stats35_2.txt

#3.Remove reads with ave score <35 (Phred score =20 99% accuracy 1% chance of mistake)
bbduk.sh in=~/programs/HCV/Output3A/Bash/D75365-/D75365-_unmatched.Q35.fq out=~/programs/HCV/Output3A/clean.fq/D75365-_clean.Q35.fq maq=35 stats=~/programs/HCV/Output3A/Bash/D75365-/stats35_3.txt

#4.deduplication
clumpify.sh in=~/programs/HCV/Output3A/clean.fq/D75365-_clean.Q35.fq out=Output3A/clumped.fq/D75365-_clumped.fq dedupe subs=0

#5.merge
bbmerge.sh in=Output3A/clumped.fq/D75365-_clumped.fq out=Output3A/merged/D75365-_merged.fq outu=Output3A/unmerged/D75365-_unmerged.fq ihist=Output3A/Bash/D75365-/D75365-_ihist.txt

#6 Map to the reference to create consensus 


bwa index -p D75365- -a is ~/programs/HCV/Output3A/Consensus/D75365-_consensus.fasta

bwa mem -t 8 -k 15 D75365- ~/programs/HCV/Output3A/merged/D75365-_merged.fq  > ~/programs/HCV/Output3A/Bash/D75365-/D75365-_ConMapped.merge.sam
bwa mem -t 8 -k 15 D75365- ~/programs/HCV/Output3A/unmerged/D75365-_unmerged.fq  > ~/programs/HCV/Output3A/Bash/D75365-/D75365-_ConMapped.un.sam


#2 Hard clipped the aligned reads
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output3A/Bash/D75365-/D75365-_ConMapped.merge.sam > ~/programs/HCV/Output3A/sam/D75365-_me_clipped.sam
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output3A/Bash/D75365-/D75365-_ConMapped.un.sam > ~/programs/HCV/Output3A/sam/D75365-_un_clipped.sam

rm ~/programs/HCV/Output3A/Bash/D75365-/D75365-_ConMapped.merge.sam ~/programs/HCV/Output3A/Bash/D75365-/D75365-_ConMapped.un.sam

