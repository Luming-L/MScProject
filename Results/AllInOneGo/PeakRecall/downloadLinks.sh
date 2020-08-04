#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd -V                  
#$ -l h_rt=01:30:00
#$ -l h_vmem=8G
 
#$ -t 14-14

for link in $(cat /home/s1949868/MScProject/Results/AllInOneGo/PeakRecall/bigWigLinks.txt | head -n $SGE_TASK_ID | tail -n1); do 
	wget $link
done
