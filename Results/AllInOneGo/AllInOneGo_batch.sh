#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd -V                  
#$ -l h_rt=00:30:00
#$ -pe sharedmem 4
#$ -l h_vmem=4G

#$ -t 1-1

# Configure modules
. /etc/profile.d/modules.sh

for file in $(ls /exports/eddie/scratch/s1949868/BigWig/* | head -n $SGE_TASK_ID | tail -n1); do
	echo $file
	fileName=`basename -s ".bw" $file`
	echo $fileName

	# make a directory for a sample and enter it
	mkdir ./$fileName
	cd ./$fileName

	# run snakemake workflow
	snakemake -s /exports/eddie/scratch/s1949868/AllInOneGo/AllInOneGo.snakemake ${fileName}.bg -j1 -p

done
