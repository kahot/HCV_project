#!/bin/bash


mkdir Output1B/Bash/D75000
 
#0 adapter trimming
bbduk.sh in1=~/programs/HCV/Data/fastq_1B/D75000-HCV_S11_L001_R1_001.fastq in2=~/programs/HCV/Data/fastq_1B/D75000-HCV_S11_L001_R2_001.fastq  out=~/programs/HCV/Output1B/Bash/D75000/D75000_adp.trimmed.fastq ref=/Users/kahotisthammer/programs/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=~/programs/HCV/Output1B/Bash/D75000/stats35_0.txt

#1 Trim reads at both ends at Score<15
bbduk.sh in=~/programs/HCV/Output1B/Bash/D75000/D75000_adp.trimmed.fastq out=~/programs/HCV/Output1B/Bash/D75000/D75000_trimmed.Q35.fastq qtrim=rl trimq=35 stats=~/programs/HCV/Output1B/Bash/D75000/stats35_1.txt

#2. Kmer filtering
bbduk.sh in=~/programs/HCV/Output1B/Bash/D75000/D75000_trimmed.Q35.fastq out=~/programs/HCV/Output1B/Bash/D75000/D75000_unmatched.Q35.fq ref=/Users/kahotisthammer/programs/bbmap/resources/phix174_ill.ref.fa k=31 hdist=1 stats=~/programs/HCV/Output1B/Bash/D75000/stats35_2.txt

#3.Remove reads with ave score <35 (Phred score =20 99% accuracy 1% chance of mistake)
bbduk.sh in=~/programs/HCV/Output1B/Bash/D75000/D75000_unmatched.Q35.fq out=~/programs/HCV/Output1B/clean.fq/D75000_clean.Q35.fq maq=35 stats=~/programs/HCV/Output1B/Bash/D75000/stats35_3.txt

#4.deduplication
clumpify.sh in=~/programs/HCV/Output1B/clean.fq/D75000_clean.Q35.fq out=Output1B/clumped.fq/D75000_clumped.fq dedupe subs=0

#5.merge
bbmerge.sh in=Output1B/clumped.fq/D75000_clumped.fq out=Output1B/merged/D75000_merged.fq outu=Output1B/unmerged/D75000_unmerged.fq ihist=Output1B/Bash/D75000/D75000_ihist.txt

bwa index -p D75000 -a is ~/programs/HCV/Output1B/Consensus/D75000_consensus.fasta

bwa mem -t 8 -k 15 D75000 ~/programs/HCV/Output1B/merged/D75000_merged.fq  > ~/programs/HCV/Output1B/Bash/D75000/D75000_ConMapped.merge.sam
bwa mem -t 8 -k 15 D75000 ~/programs/HCV/Output1B/unmerged/D75000_unmerged.fq  > ~/programs/HCV/Output1B/Bash/D75000/D75000_ConMapped.un.sam


#2 Hard clipped the aligned reads
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75000/D75000_ConMapped.merge.sam > ~/programs/HCV/Output1B/sam/D75000_me_clipped.sam
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75000/D75000_ConMapped.un.sam > ~/programs/HCV/Output1B/sam/D75000_un_clipped.sam

rm ~/programs/HCV/Output1B/Bash/D75000/D75000_ConMapped.merge.sam ~/programs/HCV/Output1B/Bash/D75000/D75000_ConMapped.un.sam

