#!/usr/bin/python3
import numpy as np
import pandas as pd
import os, sys
import subprocess as sp

# store the path of the bedGraph file
file_path  = sys.argv[1]
print(file_path)

# get the name of the bedGraph file without leading directory components and a trailing SUFFIX.
cmd='basename -s ".bg" '+file_path
result = sp.check_output(cmd, shell=True)
fileName = str(result).split("'",1)[1].split("\\",1)[0]
print(fileName)

# read bedgraph file containing ATAC-seq signal tracks
df = pd.read_csv(file_path, sep="\t",header=None,names=["chr","start","end","score"])

# obtain the length of each chromosome, and use a dictionary to store them
Chrs={}
for chrName in df['chr'].drop_duplicates().values.tolist():
    Chrs[chrName]=df[df['chr']==chrName].iloc[-1:,]['end'].values[0]

# For each track in the bedGraph, the score is the signal in each region (usually 100-bp bin)
# sum up the length of each chromosome as genome length
# calculate (sum_of_signals_in_all_bins/genome_zise)*bin_size, as bias from the whole genome background
# It will be the lambda in poisson distribution
Lambda = df['score'].sum()/sum(Chrs.values())*100
print(df['score'].sum())
print(sum(Chrs.values()))

# the first column is chromosome name
s1 = pd.Series(list(Chrs.keys()))
# the second column is the first position in each chromosome
s2 = pd.Series(np.zeros(len(Chrs.keys()),int))
# the third column is the last position in each chromosome
s3 = pd.Series(list(Chrs.values()))
# the fourth column is Lambda
s4 = pd.Series(np.full(len(Chrs.keys()),Lambda))

# make a dataframe to store tracks
df2=pd.DataFrame({'chr':s1,'start':s2,'end':s3,'Lambda':s4})
df2 = df2[['chr','start','end','Lambda']]

df2.to_csv(fileName+".lambda.bg",sep="\t",header=False,index=False)

# Compare ATAC-seq signal and local lambda to get the scores in pvalue
os.system("macs2 bdgcmp -t "+file_path+" -c " +fileName+".lambda.bg"+" -m ppois -o "+fileName+".pvalue.bg")
print("macs2 bdgcmp -t "+file_path+" -c " +fileName+".lambda.bg"+" -m ppois -o "+fileName+".pvalue.bg")
# Call peaks on score track using a cutoff p-value=0.001
os.system("macs2 bdgpeakcall -i "+fileName+".pvalue.bg"+" -c 2 -l 150 -g 75 -o "+fileName+".peaks.bed")
print("macs2 bdgpeakcall -i "+fileName+".pvalue.bg"+" -c 2 -l 150 -g 75 -o "+fileName+".peaks.bed")