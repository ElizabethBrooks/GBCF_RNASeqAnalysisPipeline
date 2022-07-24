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
setwd("/Users/bamflappy/GBCF/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq")
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq/subset")
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq")
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq/subset")

#Import gene count data
#inputTable <- read.csv(file=args[1], header = TRUE, sep = "\t", row.names="gene")
inputTable <- read.table(file="220705_Yoon_Adipocyte_Pool2_RNAseq_merged_counts_formatted.txt", header = TRUE, sep = "\t", row.names="gene")
#inputTable <- read.table(file="220705_Yoon_Adipocyte_Pool2_RNAseq_merged_counts_formatted.txt", header = TRUE, sep = "\t", row.names="gene")[,1:12]
#inputTable <- read.table(file="220707_Yoon_Jurkat_Pool1_RNAseq_merged_counts_formatted.txt", header = TRUE, sep = "\t", row.names="gene")
#inputTable <- read.table(file="220707_Yoon_Jurkat_Pool1_RNAseq_merged_counts_formatted.txt", header = TRUE, sep = "\t", row.names="gene")[,1:12]

#Trim the data table
countsTable <- head(inputTable, - 5)

#Set number of samples
numSamples <- 13
#numSamples <- 12

#Add grouping factor
group <- factor(c(rep("100mV",4), rep("180mV",4), rep("CTL",4), "Undetermined"))
#group <- factor(c(rep("100mV",4), rep("180mV",4), rep("CTL",4)))

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

#Setup color vector for plotting
zis_subset <- c("#0000FF", "#000000", "#FF0000")

#Perform an exact test for 100mV vs CTL
tested <- exactTest(list, pair=c("CTL", "100mV"))
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table
write.table(resultsTbl, file="100mV_CTL.csv", sep=",", row.names=TRUE)
#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
jpeg("100mV_CTL_plotMD.jpg")
plotMD(tested)
abline(h=c(-1, 1), col="blue")
dev.off()
#Identify significantly DE genes
resultsTbl$topDE <- "NA"
resultsTbl$topDE[resultsTbl$logFC > 1 & resultsTbl$FDR < 0.05] <- "UP"
resultsTbl$topDE[resultsTbl$logFC < -1 & resultsTbl$FDR < 0.05] <- "DOWN"
#Create volcano plot
jpeg("100mV_CTL_plotVolcano.jpg")
ggplot(data=resultsTbl, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = zis_subset)
dev.off()

#Perform an exact test for 180mV vs CTL
tested <- exactTest(list, pair=c("CTL", "180mV"))
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table
write.table(resultsTbl, file="180mV_CTL.csv", sep=",", row.names=TRUE)
#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
jpeg("180mV_CTL_plotMD.jpg")
plotMD(tested)
abline(h=c(-1, 1), col="blue")
dev.off()
#Identify significantly DE genes
resultsTbl$topDE <- "NA"
resultsTbl$topDE[resultsTbl$logFC > 1 & resultsTbl$FDR < 0.05] <- "UP"
resultsTbl$topDE[resultsTbl$logFC < -1 & resultsTbl$FDR < 0.05] <- "DOWN"
#Create volcano plot
jpeg("180mV_CTL_plotVolcano.jpg")
ggplot(data=resultsTbl, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = zis_subset)
dev.off()

#Perform an exact test for 180mV vs 100mV
tested <- exactTest(list, pair=c("100mV", "180mV"))
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table
write.table(resultsTbl, file="180mV_100mV.csv", sep=",", row.names=TRUE)
#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
jpeg("180mV_100mV_plotMD.jpg")
plotMD(tested)
abline(h=c(-1, 1), col="blue")
dev.off()
#Identify significantly DE genes
resultsTbl$topDE <- "NA"
resultsTbl$topDE[resultsTbl$logFC > 1 & resultsTbl$FDR < 0.05] <- "UP"
resultsTbl$topDE[resultsTbl$logFC < -1 & resultsTbl$FDR < 0.05] <- "DOWN"
#Create volcano plot
jpeg("180mV_100mV_plotVolcano.jpg")
ggplot(data=resultsTbl, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = zis_subset)
dev.off()
