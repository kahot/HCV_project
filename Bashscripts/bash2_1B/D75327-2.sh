#!/bin/bash


#1  Align the filtered reads using bwa to the consensus sequence. 

bwa index -p D75327- -a is ~/programs/HCV/Output1B/Consensus/D75327-_consensus.fasta

bwa mem -t 8 -k 15 D75327- ~/programs/HCV/Output1B/merged/D75327-_merged.fq  > ~/programs/HCV/Output1B/Bash/D75327-/D75327-_ConMapped.merge.sam
bwa mem -t 8 -k 15 D75327- ~/programs/HCV/Output1B/unmerged/D75327-_unmerged.fq  > ~/programs/HCV/Output1B/Bash/D75327-/D75327-_ConMapped.un.sam


#2 Hard clipped the aligned reads
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75327-/D75327-_ConMapped.merge.sam > ~/programs/HCV/Output1B/sam/D75327-_me_clipped.sam
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75327-/D75327-_ConMapped.un.sam > ~/programs/HCV/Output1B/sam/D75327-_un_clipped.sam

rm ~/programs/HCV/Output1B/Bash/D75327-/D75327-_ConMapped.merge.sam ~/programs/HCV/Output1B/Bash/D75327-/D75327-_ConMapped.un.sam

