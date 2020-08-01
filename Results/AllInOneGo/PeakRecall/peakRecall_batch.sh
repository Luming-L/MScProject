#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd -V                  
#$ -l h_rt=00:10:00 # 6G 3min
#$ -l h_vmem=8G # 4G is too small
 
#$ -t 1-796
 
# Configure modules
. /etc/profile.d/modules.sh
module load igmm/apps/MACS2/2.1.1
module load python/3.4.3

# recall peaks from bedgraph
for file in $(ls /exports/eddie/scratch/s1949868/BedGraph/*insertions.bg | head -n $SGE_TASK_ID | tail -n1)
do
echo $file

echo "recall start: $(date)"
python3 /home/s1949868/peakRecall.py $file
echo "recall done: $(date)"

done 
