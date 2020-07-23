#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd -V                  
#$ -l h_rt=00:30:00
#$ -pe sharedmem 4
#$ -l h_vmem=4G

#$ -t 1-2

# Configure modules
. /etc/profile.d/modules.sh
module load python/3.4.3
module load igmm/apps/BEDTools/2.27.1
module load igmm/apps/meme/4.11.1
module load igmm/apps/MACS2/2.1.1

for file in $(ls /exports/eddie/scratch/s1949868/BigWig/* | head -n $SGE_TASK_ID | tail -n1); do
	echo $file
	fileName=`basename -s ".bw" $file`
	echo $fileName

	# make a directory for a sample and enter it
	mkdir ./$fileName
	cd ./$fileName

	# run snakemake workflow
	snakemake -s /home/s1949868/Results/AllInOneGo/AllInOneGo.snakemake ${fileName}.PRDM9_bound_peaks.bed -j1 -p

done