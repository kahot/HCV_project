

bwa mem -t 8 -k 15 D75000 ~/programs/HCV/Output/unmerged/D75000_unmerged.fq  > ~/programs/HCV/Output/unmerged/D75000_ConMapped.un.sam

#2 Hard clipped the aligned readsjava -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output/merged/D75000_ConMapped.merge.sam > ~/programs/HCV/Output/sam/D75000_me_clipped.sam
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output/unmerged/D75000_ConMapped.un.sam > ~/programs/HCV/Output/sam/D75000_un_clipped.sam

rm  ~/programs/HCV/Output/unmerged/D75000_ConMapped.un.sam
