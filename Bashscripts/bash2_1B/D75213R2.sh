#!/bin/bash


#1  Align the filtered reads using bwa to the consensus sequence. 

bwa index -p D75213R -a is ~/programs/HCV/Output1B/Consensus/D75213R_consensus.fasta

bwa mem -t 8 -k 15 D75213R ~/programs/HCV/Output1B/merged/D75213R_merged.fq  > ~/programs/HCV/Output1B/Bash/D75213R/D75213R_ConMapped.merge.sam
bwa mem -t 8 -k 15 D75213R ~/programs/HCV/Output1B/unmerged/D75213R_unmerged.fq  > ~/programs/HCV/Output1B/Bash/D75213R/D75213R_ConMapped.un.sam


#2 Hard clipped the aligned reads
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75213R/D75213R_ConMapped.merge.sam > ~/programs/HCV/Output1B/sam/D75213R_me_clipped.sam
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75213R/D75213R_ConMapped.un.sam > ~/programs/HCV/Output1B/sam/D75213R_un_clipped.sam

rm ~/programs/HCV/Output1B/Bash/D75213R/D75213R_ConMapped.merge.sam ~/programs/HCV/Output1B/Bash/D75213R/D75213R_ConMapped.un.sam

