library(dplyr)
library(plyr)

## INPUT
# get the path of counts file
CountsFile=commandArgs(T)

# get cancerType; read normalized counts and PRDM9 bound peaks
if (strsplit(CountsFile,split='m.')[[1]][2] == "txt") {
	# for cancer type-specific data
	cType <- strsplit(strsplit(CountsFile,split='CountMatrices/')[[1]][2],split='_')[[1]][1]
	# read normalized counts
	norm_ct <- read.delim(file = CountsFile,sep = "\t",header=TRUE)
	# read PRDM9 bound peaks
	PRDM9BoundPeaks_df <- read.delim(file = paste0("/exports/eddie/scratch/s1949868/SelectPRDM9BoundPeaks_23/",cType,"_PRDM9_bound_peaks.bed"),sep = "\t",header = FALSE)
} else {
	# for pan cancer data
	cType <- "pan"
	# read normalized counts
	norm_ct <- readRDS(file = CountsFile)
	norm_ct <- norm_ct[1,c(-6,-7)]
	# read PRDM9 bound peaks
	PRDM9BoundPeaks_df <- read.delim(file = "/exports/eddie/scratch/s1949868/SelectPRDM9BoundPeaks_pan/TCGA-ATAC_PanCancer_PRDM9_bound_peaks.bed",sep = "\t",header = FALSE)
}

# read ID conversion table
samples.ids <- read.delim(file = "/exports/eddie/scratch/s1949868/TCGA_identifier_mapping",sep = "\t",header=TRUE)

# read PRDM9 expression
PRDM9 <- read.delim("/exports/eddie/scratch/s1949868/geneExpression/PRDM9Expression.txt", sep = "\t")





## PROCESS
# manage samples.ids
samples.ids$bam_prefix <- as.character(samples.ids$bam_prefix)
samples.ids$Case_ID <- as.character(samples.ids$Case_ID)
# add ID
samples.ids$cancerType <- substr(samples.ids$bam_prefix,1,4)
samples.ids$sample <- substr(samples.ids$Case_ID,1,16)
samples.ids$cancerType_sample <- paste0(samples.ids$cancerType,"_",samples.ids$sample)
# rename norm_ct
norm_ct$name <- as.character(norm_ct$name)
rownames(norm_ct) <- norm_ct$name
# get PRDM9BoundPeaks
PRDM9BoundPeaks <- as.character(PRDM9BoundPeaks_df$V4)

# select counts within PRDM9 bound peaks
norm_ct_PRDM9BoundPeaks <- norm_ct[PRDM9BoundPeaks,]
# get index of sampleID in norm_ct_PRDM9BoundPeaks
idx <- match(gsub("_","-",colnames(norm_ct_PRDM9BoundPeaks)[-c(1:5)]),samples.ids$bam_prefix)
# rename colnames of norm_ct_PRDM9BoundPeaks
colnames(norm_ct_PRDM9BoundPeaks)[-c(1:5)] <- samples.ids[idx,"cancerType_sample"]

# mark replicates
duplicated.idx <- duplicated(colnames(norm_ct_PRDM9BoundPeaks)[-c(1:5)])
# mark rep1
colnames(norm_ct_PRDM9BoundPeaks)[-c(1:5)][!duplicated.idx] <- 
  paste0(colnames(norm_ct_PRDM9BoundPeaks)[-c(1:5)][!duplicated.idx],"_rep1")
# mark rep2
colnames(norm_ct_PRDM9BoundPeaks)[-c(1:5)][duplicated.idx] <- 
  paste0(colnames(norm_ct_PRDM9BoundPeaks)[-c(1:5)][duplicated.idx],"_rep2")

# calculate the mean of replicates
# This function will calculate the Means of the peaks for a given group
groupMeans <- function(df, groups, na.rm = TRUE){
  gm <- lapply(unique(groups), function(x){
    rowMeans(df[,grepl(x,colnames(df)),drop = F],na.rm=TRUE)
  }) %>% Reduce("cbind",.) # calculate rowMean and combine
  
  colnames(gm) <- unique(groups) # specify colnames to sampleID
  
  return(gm)
}
# calculate the mean of the replicates of each patient.  
matMerged <- groupMeans(df = norm_ct_PRDM9BoundPeaks, 
                        groups =  unique(samples.ids[idx,"cancerType_sample"]))
norm_ct_PRDM9BoundPeaks.merge <- matMerged

# Group
# group by PRDM9 expression
withExpression.idx <- intersect(x = rownames(PRDM9[PRDM9$PRDM9Expression!=0,,drop = FALSE]), 
                                y = colnames(norm_ct_PRDM9BoundPeaks.merge))
withoutExpression.idx <- intersect(x = rownames(PRDM9[PRDM9$PRDM9Expression==0,,drop = FALSE]), 
                                y = colnames(norm_ct_PRDM9BoundPeaks.merge))

# t-test
# with_minus_without = log2NormCountsWith-log2NormCountsWithout = log2(With/Without)
result <- plyr::adply(norm_ct_PRDM9BoundPeaks.merge,.margins = 1,.fun = function(peak){
  results <- t.test(peak[withExpression.idx],peak[withoutExpression.idx],conf.level = TRUE)
  
  return(tibble::tibble("raw_p_value"= results$p.value,
                        "with_minus_without" = results$estimate[1] - results$estimate[2]
                        ))
  
}, .progress = "time", .id = "peak")
# use Benjamini & Hochberg (1995) method ("BH" or its alias "fdr") to adjust pvalue
result$p.adj <- stats::p.adjust(result$raw_p_value,method = "BH")
rownames(result) <- result$peak
output <- cbind(norm_ct_PRDM9BoundPeaks[,c(1:3)],result)





## OUTPUT
write.table(output,
	file=paste0(cType,"_CompareCounts_WithAndWithoutPRDM9.txt"),
	sep = "\t",
	append=TRUE,row.names = FALSE,col.names = TRUE,
	quote =FALSE
	)


