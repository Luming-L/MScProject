#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd -V                  
#$ -l h_rt=00:30:00
#$ -pe sharedmem 4
#$ -l h_vmem=4G

#$ -t 1-24

# Configure modules
. /etc/profile.d/modules.sh
module load igmm/apps/R/3.6.3

for file in $(ls /exports/eddie/scratch/s1949868/CountMatrices/*2* | head -n $SGE_TASK_ID | tail -n1); do
	echo $file
	# args[1] specify the path of counts file
	# args[2] specify the threshold of PRDM9 expression to group
	Rscript /home/s1949868/Results/GroupByPRDM9Expression/CompareCounts/CompareCounts.R $file 11
done