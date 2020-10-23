#!/usr/bin/env Rscript
#Usage: Rscript exactTest_edgeR.r countsFile.csv startColumn endColumn
#Usage Ex: Rscript exactTest_edgeR.r daphnia_genesCounted.csv 1 6
#R script to perform statistical analysis of gene count tables using edgeR exact test

#Install edgeR, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")

#Load the edgeR library
library("edgeR")

#Set index range for selected genotype
start=args[2]
end=args[3]
outputName=paste(start,end,"DEResults.csv",sep="_")

#Import gene count data
countsTable <- read.csv(file=args[1], row.names="gene")[ ,start:end]
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

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)

#Perform an exact test for treat vs ctrl
tested <- exactTest(list, pair=c("ctrl", "treat"))
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table
#Write test results to csv file
write.table(resultsTbl, file=outputName, sep=",", row.names=TRUE)

