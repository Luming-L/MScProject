#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd -V                  
#$ -l h_rt=00:20:00
#$ -l h_vmem=8G

#$ -t 1-23

# Configure modules
. /etc/profile.d/modules.sh
module load igmm/apps/R/3.6.3

for file in $(ls /exports/eddie/scratch/s1949868/CountMatrices/*_log2norm.txt | head -n $SGE_TASK_ID | tail -n1); do
	Rscript ~/CompareCounts.R $file
done