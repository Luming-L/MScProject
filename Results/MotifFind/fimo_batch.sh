#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd -V                  
#$ -l h_rt=00:05:00 # 8G 29s
#$ -l h_vmem=8G # 2G is not enough for `bedtools getfasta`
 
#$ -t 1-1

# Configure modules
. /etc/profile.d/modules.sh
module load igmm/apps/meme/4.11.1
module load igmm/apps/BEDTools/2.27.1

# find PRDM9 motif occurrence
for file in $(ls /exports/eddie/scratch/s1949868/MotifFind_fimo_pan/*PeakSet.bed | head -n $SGE_TASK_ID | tail -n1); do
	echo $file
	fileName=`basename -s ".bed" $file`
	echo $fileName

	# # extract sequences from a FASTA file for each of the intervals in BED file by `bedtools getfasta`
	bedtools getfasta -fi /home/s1949868/Tools/hg38.fa -bed $file -fo "${fileName}.fasta"

	# find PRDM9 Motif Occurrences by `fimo`
	fimo  --verbosity 4 --parse-genomic-coord --max-stored-scores 10000000 --oc "${fileName}_fimo_out" /home/s1949868/PRDM9.pwm.meme "${fileName}.fasta"

done