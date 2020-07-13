# wilcoxon test for with and without group


library(plyr)


# read expression, numPRDM9MotifInPeaks and numPRDM9BoundPeaks
PRDM9.expression <- read.delim("/exports/eddie/scratch/s1949868/geneExpression/PRDM9Expression.txt", sep = "\t", header = TRUE, row.names = 1)
numPRDM9MotifInPeaks <- read.delim("/exports/eddie/scratch/s1949868/MotifFind_fimo_404/allFimoGFF_CaseID/numPRDM9MotifInPeaks.txt", sep = "\t", header = TRUE, row.names = 1)
numPRDM9BoundPeaks <- read.delim("/exports/eddie/scratch/s1949868/SelectPRDM9BoundPeaks_404/numPRDM9BoundPeaks.txt", sep = "\t", header = TRUE, row.names = 1)


# combine expression, numPRDM9MotifInPeaks and numPRDM9BoundPeaks
PRDM9 <- cbind(PRDM9.expression, numPRDM9MotifInPeaks[rownames(PRDM9.expression),], numPRDM9BoundPeaks[rownames(PRDM9.expression),])
colnames(PRDM9) <- c("PRDM9Expression", "numPRDM9MotifInPeaks", "numPRDM9BoundPeaks")
# add cancerType
PRDM9$cancerType <- substr(rownames(PRDM9),1,4)


# wilcoxon test for numPRDM9MotifInPeaks and numPRDM9BoundPeaks
results <- plyr::ldply(.data = unique(PRDM9$cancerType),.fun = function(cType){

  with.size <- nrow(PRDM9[PRDM9$cancerType == cType & PRDM9$PRDM9Expression != 0,])
  without.size <- nrow(PRDM9[PRDM9$cancerType == cType & PRDM9$PRDM9Expression == 0,])

  if (with.size > 0 & without.size > 0) {

    result1 <- wilcox.test(
      x = PRDM9[PRDM9$cancerType == cType & PRDM9$PRDM9Expression != 0, "numPRDM9MotifInPeaks"],
      y = PRDM9[PRDM9$cancerType == cType & PRDM9$PRDM9Expression == 0, "numPRDM9MotifInPeaks"],
      paired = FALSE, alternative = "greater")
    result2 <- wilcox.test(
      x = PRDM9[PRDM9$cancerType == cType & PRDM9$PRDM9Expression != 0, "numPRDM9BoundPeaks"],
      y = PRDM9[PRDM9$cancerType == cType & PRDM9$PRDM9Expression == 0, "numPRDM9BoundPeaks"],
      paired = FALSE, alternative = "greater")

    return(tibble::tibble(
      "cancerType" = cType,
      "withPRDM9Size" = with.size,
      "withoutPRDM9Size" = without.size,
      "pvalueNumPRDM9MotifInPeaks" = result1$p.value,
      "pvalueNumPRDM9BoundPeaks" = result2$p.value))
 }

}
            )


# write results
write.table(results,
  file="ComparePRDM9MotifAndPeaks_WithAndWithoutPRDM9.txt",
  sep = "\t",
  append=TRUE,row.names = FALSE,col.names = TRUE,
  quote =FALSE
  )
