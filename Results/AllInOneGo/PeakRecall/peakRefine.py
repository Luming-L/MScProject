#!/usr/bin/python

import os, sys
import subprocess as sp

# number of a biological sample
num = int(sys.argv[1])
print(num)

# some paths
Sample_PeakCalls_path = "/exports/eddie/scratch/s1949868/RefineRecalledPeaks"

# use the same ID to identify technical replicates from the same biological sample
l1 = ["_".join(i.split("_")[0:8]) for i in os.listdir(Sample_PeakCalls_path) if i.endswith(".bed.sorted")]
print(len(l1))

# count technical replicates for each ID and store them in a dict
bioSample_dict={}
for item in set(l1):
	bioSample_dict[item]=l1.count(item)
print(list(bioSample_dict.keys())[-1])
print(sorted(list(set(l1)))[-1])


# peakOverlap() function
def peakOverlap(sample, path):
	'''
	Refine recalled peaks from 796 technical replicates by using `bedtools intersect`:
	1.for each technical replicate, output cancer type-specific peaks that can be 100% overlapped by recalled peaks
	2.refined peaks in two replicates from the same sample will be merged. only report peaks obeserved in two replicates
	3.finally get a list of peaks for each of 410 biological samples
	Args:
		sample: ID for replicates from the same biological sample
		path: the path of input files (sorted Cancer type-specific peakCalls and all samples recalled peaks) and output files
	'''
	print(sample)
	if bioSample_dict[sample] == 1: # 1 technical replicate for a ID
		os.system("bedtools intersect -a " + path + "/" + sample.split("_")[0].replace('x', '') + "*txt.sorted -b " + path + "/" + sample + "*bed.sorted -f 0.5 -u > " + path + "/" + sample+"_peakCalls.bed")
		print("bedtools intersect -a " + path + "/" + sample.split("_")[0].replace('x', '') + "*txt.sorted -b " + path + "/" + sample + "*bed.sorted -f 0.5 -u > " + path + "/" + sample+"_peakCalls.bed")

	elif bioSample_dict[sample] == 2: # 2 technical replicates for a ID
		os.system("bedtools intersect -a " + path + "/" + sample.split("_")[0].replace('x', '') + "*txt.sorted -b " + path + "/" + sample + "*bed.sorted -f 0.5 -c -wa" + " | awk '{FS=OFS=" + r'"\t"' + ";if($5>0){print $1,$2,$3,$4}}'" + " > " + path + "/" + sample+"_peakCalls.bed")
		print("bedtools intersect -a " + path + "/" + sample.split("_")[0].replace('x', '') + "*txt.sorted -b " + path + "/" + sample + "*bed.sorted -f 0.5 -c -wa" + " | awk '{FS=OFS=" + r'"\t"' + ";if($5>0){print $1,$2,$3,$4}}'" + " > " + path + "/" + sample+"_peakCalls.bed")


# ======================================= Main =======================================
# call peakOverlap() function
print(sorted(list(set(l1)))[num-1]) # set and dict in python don't have order!
peakOverlap(sample = sorted(list(set(l1)))[num-1], path=Sample_PeakCalls_path)
