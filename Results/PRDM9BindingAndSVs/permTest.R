library(regioneR,lib.loc="/exports/eddie/scratch/s1949868/R/library")
library(BSgenome.Hsapiens.UCSC.hg38,lib.loc="/exports/eddie/scratch/s1949868/R/library")

# get ID from arg
ID=commandArgs(T)

# set B
SVBs <- toGRanges(paste0("/exports/eddie/scratch/s1949868/CopyNumber/SVB/",ID,".masked_cnv_breakpoints.txt"))
# set A
peaks <- toGRanges(paste0("/exports/eddie/scratch/s1949868/SelectPRDM9BoundPeaks_404/",ID,"_PRDM9_bound_peaks.bed"))

# permTest
set.seed(123)
pt <- permTest(A=peaks, B=SVBs,
	genome=BSgenome.Hsapiens.UCSC.hg38,
	randomize.function=circularRandomizeRegions, evaluate.function=numOverlaps,
	ntimes=1000, count.once=TRUE, alternative="greater", mc.set.seed=FALSE,
	force.parallel=TRUE, mc.cores=4,
	verbose=TRUE
	)

# write the results
write.table(data.frame(ID,length(SVBs),pt$numOverlaps[4]$observed,pt$numOverlaps[1]$pval),
	file="/exports/eddie/scratch/s1949868/CopyNumber/permTest_SVB/permTest_SVB.txt",
	sep = "\t",
	append=TRUE,row.names = FALSE,col.names = FALSE,
	quote =FALSE
	)

# save .Rdata
save.image(file = paste0(ID,".RData"))