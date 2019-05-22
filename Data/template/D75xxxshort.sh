#!/bin/bash


mkdir Output/Bash/D75335 
#0 adapter trimming
bbduk.sh in1=~/programs/HCV/Data/unzipped/D75335R-HCV_S11_L001_R1_001.fastq in2=~/programs/HCV/Data/unzipped/D75335R-HCV_S11_L001_R2_001.fastq  out=~/programs/HCV/Output/Bash/D75335/D75335_adp.trimmed.fastq ref=/Users/kahotisthammer/programs/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=~/programs/HCV/Output/Bash/D75335/stats30_0.txt

#1 Trim reads at both ends at Score<15
bbduk.sh in=~/programs/HCV/Output/Bash/D75335/D75335_adp.trimmed.fastq out=~/programs/HCV/Output/Bash/D75335/D75335_trimmed.q30.fastq qtrim=rl trimq=30 stats=~/programs/HCV/Output/Bash/D75335/stats30_1.txt

#2. Kmer filtering
bbduk.sh in=~/programs/HCV/Output/Bash/D75335/D75335_trimmed.q30.fastq out=~/programs/HCV/Output/Bash/D75335/D75335_unmatched.q30.fq ref=/Users/kahotisthammer/programs/bbmap/resources/phix174_ill.ref.fa k=31 hdist=1 stats=~/programs/HCV/Output/Bash/D75335/stats30_2.txt

#3.Remove reads with ave score <30 (Phred score =20 99% accuracy 1% chance of mistake)
bbduk.sh in=~/programs/HCV/Output/Bash/D75335/D75335_unmatched.q30.fq out=~/programs/HCV/Output/clean.fq/D75335_clean.q30.fq maq=30 stats=~/programs/HCV/Output/Bash/D75335/stats30_3.txt

#4. Align the file using bwa to the consensus 
bwa index -p D75335 -a is ~/programs/HCV/Output/Consensus/D75335_consensus.fasta
bwa mem -t 4 -k 15 D75335 ~/programs/HCV/Output/clean.fq/D75335_clean.q30.fq  > ~/programs/HCV/Output/Bash/D75335/D75335_consensus_mapped.sam

#5. Hard clipped the aligned reads
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output/Bash/D75335/D75335_consensus_mapped.sam > ~/programs/HCV/Output/sam/D75335_clipped.sam

rm ~/programs/HCV/Output/Bash/D75335/D75335_adp.trimmed.fastq ~/programs/HCV/Output/Bash/D75335/D75335_trimmed.q30.fastq ~/programs/HCV/Output/Bash/D75335/D75335_unmatched.q30.fq ~/programs/HCV/Output/Bash/D75335/D75335_consensus_mapped.sam




