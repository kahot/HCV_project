#!/bin/bash


#1  Align the filtered reads using bwa to the consensus sequence. 

bwa index -p D75383R -a is ~/programs/HCV/Output1B/Consensus/D75383R_consensus.fasta

bwa mem -t 8 -k 15 D75383R ~/programs/HCV/Output1B/merged/D75383R_merged.fq  > ~/programs/HCV/Output1B/Bash/D75383R/D75383R_ConMapped.merge.sam
bwa mem -t 8 -k 15 D75383R ~/programs/HCV/Output1B/unmerged/D75383R_unmerged.fq  > ~/programs/HCV/Output1B/Bash/D75383R/D75383R_ConMapped.un.sam


#2 Hard clipped the aligned reads
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75383R/D75383R_ConMapped.merge.sam > ~/programs/HCV/Output1B/sam/D75383R_me_clipped.sam
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75383R/D75383R_ConMapped.un.sam > ~/programs/HCV/Output1B/sam/D75383R_un_clipped.sam

rm ~/programs/HCV/Output1B/Bash/D75383R/D75383R_ConMapped.merge.sam ~/programs/HCV/Output1B/Bash/D75383R/D75383R_ConMapped.un.sam

