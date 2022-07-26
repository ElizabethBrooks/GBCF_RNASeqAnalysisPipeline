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
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq/subset")
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq")
setwd("/Users/bamflappy/GBCF/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq/subset")

#Set number of samples
numSamples <- 13

#Import gene count data
#inputTable <- read.csv(file=args[1], header = TRUE, sep = "\t", row.names="gene")
#inputTable <- read.table(file="220705_Yoon_Adipocyte_Pool2_RNAseq_merged_counts_formatted.txt", header = TRUE, sep = "\t", row.names="gene")
#inputTable <- read.table(file="220707_Yoon_Jurkat_Pool1_RNAseq_merged_counts_formatted.txt", header = TRUE, sep = "\t", row.names="gene")
inputTable <- read.table(file="220707_Yoon_Jurkat_Pool1_RNAseq_merged_counts_formatted.txt", header = TRUE, sep = "\t", row.names="gene")

#Trim the data table
countsTable <- head(inputTable, - 5)

#Add grouping factor
group <- factor(c(rep("100mV",4), rep("180mV",4), rep("CTL",4), "Undetermined"))

#Create DGE list object
list <- DGEList(counts=countsTable,group=group)


#Prep Stage

#Plot the library sizes before normalization
jpeg("plotBars_librarySizes.jpg")
barplot(list$samples$lib.size*1e-6, names=1:numSamples, ylab="Library size (millions)")
dev.off()
#Draw a MDS plot to show the relative similarities of the samples
# and to view batch and treatment effects before normalization
jpeg("plotMDS_beforeNormalize.jpg")
plotMDS(list, col=rep(1:3, each=4))
dev.off()

#There is no purpose in analysing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
table(keep)
list <- list[keep, , keep.lib.sizes=FALSE]
#Calculate normalized factors
list <- calcNormFactors(list)
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)
write.table(normList, file="normalized_counts.csv", sep=",", row.names=TRUE)
#View normalization factors
list$samples
dim(list)

#Draw a MDS plot to show the relative similarities of the samples
# and to view batch and treatment effects after normalization
jpeg("plotMDS_afterNormalize.jpg")
plotMDS(list, col=rep(1:3, each=4))
dev.off()
#Draw a heatmap of individual RNA-seq samples using moderated
# log-counts-per-million after normalization
jpeg("plotHeatMap_afterNormalize.jpg")
logcpm <- cpm(list, log=TRUE)
heatmap(logcpm)
dev.off()

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)
list$common.dispersion
#View dispersion estimates and biological coefficient of variation
jpeg("plotBCV.jpg")
plotBCV(list)
dev.off()
