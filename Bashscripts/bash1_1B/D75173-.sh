#!/bin/bash


mkdir Output1B/Bash/D75173-
 
#0 adapter trimming
bbduk.sh in1=~/programs/HCV/Data/fastq_1B/D75173-HCV_S23_L001_R1_001.fastq in2=~/programs/HCV/Data/fastq_1B/D75173-HCV_S23_L001_R2_001.fastq  out=~/programs/HCV/Output1B/Bash/D75173-/D75173-_adp.trimmed.fastq ref=/Users/kahotisthammer/programs/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=~/programs/HCV/Output1B/Bash/D75173-/stats35_0.txt

#1 Trim reads at both ends at Score<15
bbduk.sh in=~/programs/HCV/Output1B/Bash/D75173-/D75173-_adp.trimmed.fastq out=~/programs/HCV/Output1B/Bash/D75173-/D75173-_trimmed.Q35.fastq qtrim=rl trimq=35 stats=~/programs/HCV/Output1B/Bash/D75173-/stats35_1.txt

#2. Kmer filtering
bbduk.sh in=~/programs/HCV/Output1B/Bash/D75173-/D75173-_trimmed.Q35.fastq out=~/programs/HCV/Output1B/Bash/D75173-/D75173-_unmatched.Q35.fq ref=/Users/kahotisthammer/programs/bbmap/resources/phix174_ill.ref.fa k=31 hdist=1 stats=~/programs/HCV/Output1B/Bash/D75173-/stats35_2.txt

#3.Remove reads with ave score <35 (Phred score =20 99% accuracy 1% chance of mistake)
bbduk.sh in=~/programs/HCV/Output1B/Bash/D75173-/D75173-_unmatched.Q35.fq out=~/programs/HCV/Output1B/clean.fq/D75173-_clean.Q35.fq maq=35 stats=~/programs/HCV/Output1B/Bash/D75173-/stats35_3.txt

#4.deduplication
clumpify.sh in=~/programs/HCV/Output1B/clean.fq/D75173-_clean.Q35.fq out=Output1B/clumped.fq/D75173-_clumped.fq dedupe subs=0

#5.merge
bbmerge.sh in=Output1B/clumped.fq/D75173-_clumped.fq out=Output1B/merged/D75173-_merged.fq outu=Output1B/unmerged/D75173-_unmerged.fq ihist=Output1B/Bash/D75173-/D75173-_ihist.txt

bwa index -p D75173- -a is ~/programs/HCV/Output1B/Consensus/D75173-_consensus.fasta

bwa mem -t 8 -k 15 D75173- ~/programs/HCV/Output1B/merged/D75173-_merged.fq  > ~/programs/HCV/Output1B/Bash/D75173-/D75173-_ConMapped.merge.sam
bwa mem -t 8 -k 15 D75173- ~/programs/HCV/Output1B/unmerged/D75173-_unmerged.fq  > ~/programs/HCV/Output1B/Bash/D75173-/D75173-_ConMapped.un.sam


#2 Hard clipped the aligned reads
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75173-/D75173-_ConMapped.merge.sam > ~/programs/HCV/Output1B/sam/D75173-_me_clipped.sam
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75173-/D75173-_ConMapped.un.sam > ~/programs/HCV/Output1B/sam/D75173-_un_clipped.sam

rm ~/programs/HCV/Output1B/Bash/D75173-/D75173-_ConMapped.merge.sam ~/programs/HCV/Output1B/Bash/D75173-/D75173-_ConMapped.un.sam

