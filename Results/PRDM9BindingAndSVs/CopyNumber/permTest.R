library(regioneR,lib.loc="/exports/eddie/scratch/s1949868/R/library")
library(BSgenome.Hsapiens.UCSC.hg38,lib.loc="/exports/eddie/scratch/s1949868/R/library")

# set B
mutect2_snv <- readRDS(file = "/exports/eddie/scratch/s1949868/SNPsAndSmallINDELs/mutect2_snv.rds")
mutations <- toGRanges(mutect2_snv[mutect2_snv$Sample_ID==ID,c(3,4,5)])
# set A
peaks <- toGRanges(paste0("/exports/eddie/scratch/s1949868/PRDM9BoundPeaks_410_Case_ID/",ID,"_PRDM9_bound_peaks.bed"))

# permTest
set.seed(123)
pt <- permTest(A=peaks, B=mutations,
	genome=BSgenome.Hsapiens.UCSC.hg38,
	randomize.function=circularRandomizeRegions, evaluate.function=numOverlaps,
	ntimes=1000, count.once=TRUE, alternative="greater", mc.set.seed=FALSE,
	force.parallel=TRUE, mc.cores=4,
	verbose=TRUE
	)

# write the results
write.table(data.frame(ID,cancerType,numMut,numOverlap,pt$numOverlaps[1]$pval),
	file="/exports/eddie/scratch/s1949868/SNPsAndSmallINDELs/MutNumber.noZero.pval.txt",
	sep = "\t",
	append=TRUE,row.names = FALSE,col.names = FALSE,
	quote =FALSE
	)