#!/usr/bin/env Rscript
#R script to perform statistical analysis of gene count tables using edgeR exact test
#Usage: Rscript edgeR_projects.r mergedCountsFile
#Usage Ex: Rscript edgeR_projects.r 220705_Yoon_Adipocyte_Pool2_RNAseq_merged_counts_formatted.txt
#Usage Ex: Rscript edgeR_projects.r 220707_Yoon_Jurkat_Pool1_RNAseq_merged_counts_formatted.txt

#Install edgeR, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")

#Load libraries
library(edgeR)
library(ggplot2)

#Set working directory
setwd("/Users/bamflappy/GBCF/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq/gene_counts")
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq/gene_counts")

#Import gene count data
#inputTable <- read.csv(file=args[1], header = TRUE, sep = "\t", row.names="gene")
#Full set
inputTable <- read.table(file="220705_Yoon_Adipocyte_Pool2_RNAseq_merged_counts_formatted.txt", header = TRUE, sep = "\t", row.names="gene")
#inputTable <- read.table(file="220707_Yoon_Jurkat_Pool1_RNAseq_merged_counts_formatted.txt", header = TRUE, sep = "\t", row.names="gene")

#Set working directory
#Full set
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq/gene_counts")
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq/gene_counts")
#No undetermined
setwd("/Users/bamflappy/GBCF/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq/subset_noUndetermined")
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq/subset_noUndetermined")
#DE
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq/differential_expression")
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq/differential_expression")

#Subset the input counts
#Full set
#subsetTable <- inputTable
#No undetermined
subsetTable <- inputTable[,1:12]
#Adipocyte
#subsetTable <- inputTable[ , -which(names(inputTable) %in% c("S0", "S1", "S6", "S10"))]
#Jurkat
#subsetTable <- inputTable[ , -which(names(inputTable) %in% c("S0", "S4", "S7", "S9"))]

#Trim the data table
countsTable <- head(subsetTable, - 5)

#Set number of samples
#Full set
#numSamples <- 13
#No undetermined
numSamples <- 12
#DE subset
#numSamples <- 9

#Add grouping factor
#Full set
#group <- factor(c(rep("100mV",4), rep("180mV",4), rep("CTL",4), "Undetermined"))
#No undetermined
group <- factor(c(rep("100mV",4), rep("180mV",4), rep("CTL",4)))
#DE subset
#group <- factor(c(rep("100mV",3), rep("180mV",3), rep("CTL",3)))

#Create DGE list object
list <- DGEList(counts=countsTable,group=group)


#Prep Stage

#Retrieve library sizes
libraries <- data.frame(
  samples = names(countsTable),
  sizes = list$samples$lib.size*1e-6)
#Plot the library sizes before normalization
#barplot(list$samples$lib.size*1e-6, names=1:numSamples, ylab="Library size (millions)")
jpeg("plotBars_librarySizes.jpg")
ggplot(data = libraries, aes(x = samples, y = sizes)) + 
  geom_bar(stat="identity") +
  labs(x = "Sample", y="Library size (millions)")
dev.off()

#There is no purpose in analyzing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
list <- list[keep, , keep.lib.sizes=FALSE]
#Calculate normalized factors
list <- calcNormFactors(list)
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)
write.table(normList, file="normalized_counts.csv", sep=",", row.names=TRUE)

#Draw a MDS plot to show the relative similarities of the samples
# and to view batch and treatment effects after normalization
jpeg("plotMDS_afterNormalize.jpg")
plotMDS(list, col=rep(1:3, each=3))
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
#View dispersion estimates and biological coefficient of variation
jpeg("plotBCV.jpg")
plotBCV(list)
dev.off()
