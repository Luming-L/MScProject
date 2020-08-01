#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd -V                  
#$ -l h_rt=00:30:00
#$ -l h_vmem=8G

#$ -t 1-83

# echo -e "sampleID\tnumSVBs\tnumOverlap\tpval" > permTest_SVB.txt

# Configure modules
. /etc/profile.d/modules.sh
module load igmm/apps/R/3.6.3

for SVBFile in $(ls /exports/eddie/scratch/s1949868/CopyNumber/permTest_SVB/BRCA_CHOL_PCPG/*.masked_cnv_breakpoints.txt | head -n $SGE_TASK_ID | tail -n1); do 
	echo $SVBFile; 
	sampleID=${SVBFile#*BRCA_CHOL_PCPG/}; 
	sampleID=${sampleID%.masked*};
	echo $sampleID

	Rscript /home/s1949868/MScProject/Results/PRDM9BindingAndSVs/permTest.R $sampleID

done