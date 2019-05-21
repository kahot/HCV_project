#!/bin/bash


#1  Align the filtered reads using bwa to the consensus sequence. 

bwa index -p D75224- -a is ~/programs/HCV/Output1B/Consensus/D75224-_consensus.fasta

bwa mem -t 8 -k 15 D75224- ~/programs/HCV/Output1B/merged/D75224-_merged.fq  > ~/programs/HCV/Output1B/Bash/D75224-/D75224-_ConMapped.merge.sam
bwa mem -t 8 -k 15 D75224- ~/programs/HCV/Output1B/unmerged/D75224-_unmerged.fq  > ~/programs/HCV/Output1B/Bash/D75224-/D75224-_ConMapped.un.sam


#2 Hard clipped the aligned reads
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75224-/D75224-_ConMapped.merge.sam > ~/programs/HCV/Output1B/sam/D75224-_me_clipped.sam
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75224-/D75224-_ConMapped.un.sam > ~/programs/HCV/Output1B/sam/D75224-_un_clipped.sam

rm ~/programs/HCV/Output1B/Bash/D75224-/D75224-_ConMapped.merge.sam ~/programs/HCV/Output1B/Bash/D75224-/D75224-_ConMapped.un.sam

