#!/bin/bash


#1  Align the filtered reads using bwa to the consensus sequence. 

bwa index -p D75079- -a is ~/programs/HCV/Output1B/Consensus/D75079-_consensus.fasta

bwa mem -t 8 -k 15 D75079- ~/programs/HCV/Output1B/merged/D75079-_merged.fq  > ~/programs/HCV/Output1B/Bash/D75079-/D75079-_ConMapped.merge.sam
bwa mem -t 8 -k 15 D75079- ~/programs/HCV/Output1B/unmerged/D75079-_unmerged.fq  > ~/programs/HCV/Output1B/Bash/D75079-/D75079-_ConMapped.un.sam


#2 Hard clipped the aligned reads
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75079-/D75079-_ConMapped.merge.sam > ~/programs/HCV/Output1B/sam/D75079-_me_clipped.sam
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75079-/D75079-_ConMapped.un.sam > ~/programs/HCV/Output1B/sam/D75079-_un_clipped.sam

rm ~/programs/HCV/Output1B/Bash/D75079-/D75079-_ConMapped.merge.sam ~/programs/HCV/Output1B/Bash/D75079-/D75079-_ConMapped.un.sam

