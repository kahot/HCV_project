#!/bin/bash

head -n 3 ~/programs/HCV/Output1B/sam/D75366-_me_clipped.sam > h3
cat h3 ~/programs/HCV/Output1B/sam2/D75366-_me-filtered.sam > ~/programs/HCV/Output1B/sam2/D75366-_mefiltered2.sam

#convert to bam file
samtools view -S -b ~/programs/HCV/Output1B/sam2/D75366-_mefiltered2.sam > ~/programs/HCV/Output1B/bam2/D75366-_me_filtered.bam

#sort bam file
samtools sort ~/programs/HCV/Output1B/bam2/D75366-_me_filtered.bam -o ~/programs/HCV/Output1B/bam2/D75366-_me.sort.bam

#index the bam file
samtools index ~/programs/HCV/Output1B/bam2/D75366-_me.sort.bam ~/programs/HCV/Output1B/bam2/D75366-_me.sort.bam.bai

rm ~/programs/HCV/Output1B/sam2/D75366-_mefiltered2.sam ~/programs/HCV/Output1B/bam2/D75366-_me_filtered.bam


head -n 3 ~/programs/HCV/Output1B/sam/D75366-_un_clipped.sam > h3
cat h3 ~/programs/HCV/Output1B/sam2/D75366-_un-filtered.sam > ~/programs/HCV/Output1B/sam2/D75366-_unfiltered2.sam

#convert to bam file
samtools view -S -b ~/programs/HCV/Output1B/sam2/D75366-_unfiltered2.sam > ~/programs/HCV/Output1B/bam2/D75366-_un_filtered.bam

#sort bam file
samtools sort ~/programs/HCV/Output1B/bam2/D75366-_un_filtered.bam -o ~/programs/HCV/Output1B/bam2/D75366-_un.sort.bam

#index the bam file
samtools index ~/programs/HCV/Output1B/bam2/D75366-_un.sort.bam ~/programs/HCV/Output1B/bam2/D75366-_un.sort.bam.bai

rm ~/programs/HCV/Output1B/sam2/D75366-_unfiltered2.sam ~/programs/HCV/Output1B/bam2/D75366-_un_filtered.bam
