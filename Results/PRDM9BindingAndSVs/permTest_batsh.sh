#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd -V                  
#$ -l h_rt=00:30:00
#$ -l h_vmem=8G

#$ -t 1-299

# Configure modules
. /etc/profile.d/modules.sh
module load igmm/apps/R/3.6.3

for sampleID in $(cut -f 1 /exports/eddie/scratch/s1949868/SNPsAndSmallINDELs/MutNumber.noZero.sorted.txt | head -n $SGE_TASK_ID | tail -n1)
do
	echo $sampleID
	Rscript ~/permTest.R $sampleID

done