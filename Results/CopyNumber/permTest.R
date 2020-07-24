library(regioneR,lib.loc="/exports/eddie/scratch/s1949868/R/library")
library(BSgenome.Hsapiens.UCSC.hg38,lib.loc="/exports/eddie/scratch/s1949868/R/library")

# get ID from arg
ID=commandArgs(T)
# get cancerType,numMut,numOverlap belonging to this ID
info=read.delim("/exports/eddie/scratch/s1949868/SNPsAndSmallINDELs/MutNumber.noZero.sorted.txt",header = FALSE)
cancerType=as.character(info[info$V1==ID,2])
numMut=as.character(info[info$V1==ID,3])
numOverlap=as.character(info[info$V1==ID,4])

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