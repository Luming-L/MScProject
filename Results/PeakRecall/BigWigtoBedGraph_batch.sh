#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd -V                  
#$ -l h_rt=00:10:00
#$ -l h_vmem=8G

#$ -t 1-796

for file in $(ls /exports/eddie/scratch/s1949868/BigWig/* | head -n $SGE_TASK_ID | tail -n1); do
	echo $file
	fileName=`basename -s ".bw" $file`
	echo $fileName

	/home/s1949868/Tools/bigWigToBedGraph $file /exports/eddie/scratch/s1949868/BedGraph/${fileName}.bg
done


