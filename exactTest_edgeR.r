#!/usr/bin/env Rscript
#Usage: Rscript exactTest_edgeR.r countsFile.csv startColumn endColumn
#Usage Ex: Rscript exactTest_edgeR.r daphnia_countsFile.csv 1 6
#R script to perform statistical analysis of gene count tables using edgeR exact test

#Install edgeR, this should only need to be done once
#bioLite("edgeR")

#Load the edgeR library
library("edgeR")

#Import gene count data
countsTable <- read.csv(file=args[1], row.names="gene")[ ,args[2]:args[3]]
#Add grouping factor
group <- factor(c(rep("ctrl",3),rep("treat",3)))
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)

#There is no purpose in analysing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
list <- list[keep, , keep.lib.sizes=FALSE]
#Calculate normalized factors
list <- calcNormFactors(list)
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)
write.table(normList, file="normalizedCounts.csv", sep=",", row.names=TRUE)

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)

#Perform an exact test for treat vs ctrl
tested <- exactTest(list, pair=c("ctrl", "treat"))
topTags(tested)
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table
#Output resulting table
write.table(resultsTbl, file="exactTestResults.csv", sep=",", row.names=TRUE)
