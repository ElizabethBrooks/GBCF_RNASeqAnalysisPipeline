#!/usr/bin/env Rscript
#R script to perform statistical analysis of gene count tables using edgeR exact test
#Usage: Rscript edgeR_projects.r mergedCountsFile
#Usage Ex: Rscript edgeR_projects.r 220705_Yoon_Adipocyte_Pool2_RNAseq_merged_counts_formatted.txt
#Usage Ex: Rscript edgeR_projects.r 220707_Yoon_Jurkat_Pool1_RNAseq_merged_counts_formatted.txt

#Install edgeR, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")

#Load the edgeR library
library("edgeR")

#Set working directory
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq")
setwd("/Users/bamflappy/GBCF/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq/subset")
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq")

#Import gene count data
#inputTable <- read.csv(file=args[1], header = TRUE, sep = "\t", row.names="gene")
#inputTable <- read.table(file="220705_Yoon_Adipocyte_Pool2_RNAseq_merged_counts_formatted.txt", header = TRUE, sep = "\t", row.names="gene")
inputTable <- read.table(file="220705_Yoon_Adipocyte_Pool2_RNAseq_merged_counts_formatted.txt", header = TRUE, sep = "\t", row.names="gene")[,1:12]
#inputTable <- read.table(file="220707_Yoon_Jurkat_Pool1_RNAseq_merged_counts_formatted.txt", header = TRUE, sep = "\t", row.names="gene")

#Trim the data table
countsTable <- head(inputTable, - 5)

#Set number of samples
#numSamples <- 13
numSamples <- 12

#Add grouping factor
#group <- factor(c(rep("100mV",4), rep("180mV",4), rep("CTL",4), "Undetermined"))
group <- factor(c(rep("100mV",4), rep("180mV",4), rep("CTL",4)))

#Create DGE list object
list <- DGEList(counts=countsTable,group=group)


#Prep Stage

#Plot the library sizes before normalization
jpeg("exactTest_plotBarsBefore.jpg")
barplot(list$samples$lib.size*1e-6, names=1:numSamples, ylab="Library size (millions)")
dev.off()
#Draw a MDS plot to show the relative similarities of the samples
# and to view batch and treatment effects before normalization
jpeg("exactTest_plotMDSBefore.jpg")
plotMDS(list, col=rep(1:3, each=4))
dev.off()
#Draw a heatmap of individual RNA-seq samples using moderated 
# log-counts-per-million before normalization
#jpeg("exactTest_plotHeatMapBefore.jpg")
#logcpm <- cpm(list, log=TRUE)
#heatmap(logcpm)
#dev.off()

#There is no purpose in analysing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
table(keep)
list <- list[keep, , keep.lib.sizes=FALSE]
#Calculate normalized factors
list <- calcNormFactors(list)
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)
write.table(normList, file="exactTest_normalizedCounts.csv", sep=",", row.names=TRUE)
#View normalization factors
list$samples
dim(list)

#Plot the library sizes after normalization
jpeg("exactTest_plotBarsAfter.jpg")
barplot(list$samples$lib.size*1e-6, names=1:numSamples, ylab="Library size (millions)")
dev.off()
#Draw a MDS plot to show the relative similarities of the samples
# and to view batch and treatment effects after normalization
jpeg("exactTest_plotMDSAfter.jpg")
plotMDS(list, col=rep(1:3, each=4))
dev.off()
#Draw a heatmap of individual RNA-seq samples using moderated
# log-counts-per-million after normalization
jpeg("exactTest_plotHeatMapAfter.jpg")
logcpm <- cpm(list, log=TRUE)
heatmap(logcpm)
dev.off()

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)
list$common.dispersion
#View dispersion estimates and biological coefficient of variation
jpeg("exactTest_plotBCV.jpg")
plotBCV(list)
dev.off()


#DEA Stage

#Perform an exact test for treat vs ctrl
tested <- exactTest(list, pair=c("ctrl", "treat"))
write.table(tested, file="exactTest.csv", sep=",", row.names=TRUE)
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table
write.table(resultsTbl, file="exactTest_topTags.csv", sep=",", row.names=TRUE)
#Create filtered results table of DE genes
resultsTbl.keep <- resultsTbl$FDR <= fdrCut
resultsTblFiltered <- resultsTbl[resultsTbl.keep,]
write.table(resultsTblFiltered, file="exactTest_topTags_filtered.csv", sep=",", row.names=TRUE)

#Look at the counts-per-million in individual samples for the top genes
o <- order(tested$table$PValue)
cpm(list)[o[1:10],]
#View the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested))

#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
jpeg("exactTest_plotMD.jpg")
plotMD(tested)
abline(h=c(-1, 1), col="blue")
dev.off()

#Make a mean-difference plot of two libraries of count data with smearing of points
#  with very low counts, especially those that are zero for one of the columns
jpeg("exactTest_plotMA.jpg")
plotSmear(tested)
dev.off()
