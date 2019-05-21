#!/bin/bash


#1  Align the filtered reads using bwa to the consensus sequence. 

bwa index -p D75064- -a is ~/programs/HCV/Output1B/Consensus/D75064-_consensus.fasta

bwa mem -t 8 -k 15 D75064- ~/programs/HCV/Output1B/merged/D75064-_merged.fq  > ~/programs/HCV/Output1B/Bash/D75064-/D75064-_ConMapped.merge.sam
bwa mem -t 8 -k 15 D75064- ~/programs/HCV/Output1B/unmerged/D75064-_unmerged.fq  > ~/programs/HCV/Output1B/Bash/D75064-/D75064-_ConMapped.un.sam


#2 Hard clipped the aligned reads
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75064-/D75064-_ConMapped.merge.sam > ~/programs/HCV/Output1B/sam/D75064-_me_clipped.sam
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output1B/Bash/D75064-/D75064-_ConMapped.un.sam > ~/programs/HCV/Output1B/sam/D75064-_un_clipped.sam

rm ~/programs/HCV/Output1B/Bash/D75064-/D75064-_ConMapped.merge.sam ~/programs/HCV/Output1B/Bash/D75064-/D75064-_ConMapped.un.sam

