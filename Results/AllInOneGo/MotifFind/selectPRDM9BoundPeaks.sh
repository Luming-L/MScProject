#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd -V                  
#$ -l h_rt=00:05:00 # 4G 1s
#$ -l h_vmem=4G

#$ -t 1-23

# Configure modules
. /etc/profile.d/modules.sh
module load igmm/apps/BEDTools/2.27.1

# select ATAC-seq peaks that contain any PRDM9 motif for each biological sample
for file in $(ls /exports/eddie/scratch/s1949868/MotifFind_fimo_23/allFimoGFF_CaseID_23/*fimo.gff | head -n $SGE_TASK_ID | tail -n1); do
	echo $file

	# obtain the sample name to name the output file
	fileName=`echo ${file#*allFimoGFF_CaseID_23/}`
	fileName=`echo ${fileName%_peakCalls_fimo*}`
	echo $fileName

	# select ATAC-seq peaks that contain any PRDM9 motif by `bedtools intersect`
	bedtools intersect -a /exports/eddie/scratch/s1949868/TCGA-ATAC_Cancer_Type-specific_PeakCalls/sorted/"${fileName}_peakCalls.txt.sorted" -b $file -F 1.0 -u > "${fileName}_PRDM9_bound_peaks.bed"

done