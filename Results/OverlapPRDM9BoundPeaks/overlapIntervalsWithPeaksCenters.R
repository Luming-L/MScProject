### Input
# get the filenames of intervalsFile and peakCentersFile
# args[1]: intervalsFile
# args[2]: peakCentersFile
args=commandArgs(T)
intervalsFile=args[1]
peakCentersFile=args[2]


# get genome size
genomeSize = 3031070000

# chrs that will be kept
Chr <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")


## read and process files
# read peakCentersFile
A_peakCenters_file <- read.delim(file = peakCentersFile, sep = "\t", header = FALSE)
colnames(A_peakCenters_file) <- c("chrom","start","end")
# filter chromosomes, only keep chr1-22 and chrX
A_peakCenters_file <- A_peakCenters_file[A_peakCenters_file$chrom %in% Chr,]

# write peakCenters left after filtering
write.table(A_peakCenters_file,
	file=paste0(peakCentersFile,".temp"),
	sep = "\t",
	append=FALSE,row.names = FALSE,col.names = FALSE,
	quote =FALSE
	)

# read intervalsFile
B_intervals_file <- read.delim(file = intervalsFile, sep = "\t", header = FALSE)
colnames(B_intervals_file) <- c("chrom","start","end")
# filter chromosomes, only keep chr1-22 and chrX
B_intervals_file <- B_intervals_file[B_intervals_file$chrom %in% Chr,]
# write intervals left after filtering
write.table(B_intervals_file,
	file=paste0(intervalsFile,".temp"),
	sep = "\t",
	append=FALSE,row.names = FALSE,col.names = FALSE,
	quote =FALSE
	)

# read filtered peakCenters
A_peakCenters <- read.delim(file = paste0(peakCentersFile,".temp"), sep = "\t", header = FALSE)
colnames(A_peakCenters) <- c("chrom","start","end")
# read filtered intervals
B_intervals <- read.delim(file = paste0(intervalsFile,".temp"), sep = "\t", header = FALSE)
colnames(B_intervals) <- c("chrom","start","end")





### compute the corrected overlap fraction

# this method is only suitable when c is smaller than f.
# t: total number of genomic intervals; wi: width of each interval
# n: total number of single-base peak centers
# g: genome size
# c: chance overlaps between n centers and t intervals
# f: observed overlaps between n peakCenters and t intervals
# o: corrected overlaps

## calculate chance overlaps
t <- nrow(B_intervals)
n <- nrow(A_peakCenters)
g <- genomeSize
c <- sum(1 - ((g - (B_intervals$end - B_intervals$start))/g)^n)

## get observed overlaps
f <- as.integer(system(paste0("bedtools intersect -a ", paste0(intervalsFile,".temp"), " -b ", paste0(peakCentersFile,".temp"), " -u | wc -l"),intern = TRUE))

## compute the corrected overlap fraction
# (1−f/t)=(1−o/t)(1−c/t)
# o/t=1−(1−f/t)/(1−c/t)
if (c<f) {
  correctedOverlapFraction <- 1 - (1-f/t)/(1-c/t)
} else {
  correctedOverlapFraction <- NA
}





### Output
result <- data.frame(
	"numPeakCenters" = n,
	"numIntervals" = t,
	"numObservedOverlaps" = f,
	"numChanceOverlaps" = c,
	"correctedOverlapFraction" = correctedOverlapFraction)

write.table(result,
	file = paste0("overlap_",gsub('.bed','',intervalsFile),"_VS_",gsub('.bed','',peakCentersFile),".txt"),
	sep = "\t",
	append=FALSE,row.names = FALSE,col.names = TRUE,
	quote =FALSE
	)

# remove temporary files
system("rm *.temp")



