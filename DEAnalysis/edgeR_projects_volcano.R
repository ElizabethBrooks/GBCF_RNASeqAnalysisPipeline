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
#No undetermined
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq/subset_noUndetermined")
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq/subset_noUndetermined")
#DE
setwd("/Users/bamflappy/GBCF/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq/differential_expression")
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq/differential_expression")

#Subset the input counts
#Full set
#subsetTable <- inputTable
#No undetermined
#subsetTable <- inputTable[,1:12]
#DE Adipocyte
subsetTable <- inputTable[ , -which(names(inputTable) %in% c("S00", "S01", "S06", "S10"))]
#DE Jurkat
#subsetTable <- inputTable[ , -which(names(inputTable) %in% c("S00", "S04", "S07", "S09"))]

#Trim the data table
countsTable <- head(subsetTable, - 5)

#Set number of samples
#Full set
#numSamples <- 13
#No undetermined
#numSamples <- 12
#DE subset
numSamples <- 9

#Add grouping factor
#Full set
#group <- factor(c(rep("100mV",4), rep("180mV",4), rep("CTL",4), "Undetermined"))
#No undetermined
#group <- factor(c(rep("100mV",4), rep("180mV",4), rep("CTL",4)))
#DE subset
group <- factor(c(rep("100mV",3), rep("180mV",3), rep("CTL",3)))

#Create DGE list object
list <- DGEList(counts=countsTable,group=group)


#Prep Stage

#Retrieve library sizes
libraries <- data.frame(
  samples = names(countsTable),
  sizes = list$samples$lib.size*1e-6)

#There is no purpose in analyzing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
list <- list[keep, , keep.lib.sizes=FALSE]
#Calculate normalized factors
list <- calcNormFactors(list)
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)

#Retrieve log-counts-per-million after normalization
logcpm <- cpm(list, log=TRUE)

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)


#DEA Stage

#Setup color vector for plotting
color_subset <- c("#0000FF", "#000000", "#FF0000")

#Perform an exact test for 100mV vs CTL
tested <- exactTest(list, pair=c("CTL", "100mV"))
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table
#Identify significantly DE genes
resultsTbl$sigDE <- "NA"
resultsTbl$sigDE[resultsTbl$logFC > 1 & resultsTbl$FDR < 0.05] <- "UP"
resultsTbl$sigDE[resultsTbl$logFC < -1 & resultsTbl$FDR < 0.05] <- "DOWN"
#Create volcano plot
jpeg("100mV_CTL_plotVolcano.jpg")
ggplot(data=resultsTbl, aes(x=logFC, y=-log10(FDR), color = sigDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = color_subset)
dev.off()

#Perform an exact test for 180mV vs CTL
tested <- exactTest(list, pair=c("CTL", "180mV"))
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table
#Identify significantly DE genes
resultsTbl$sigDE <- "NA"
resultsTbl$sigDE[resultsTbl$logFC > 1 & resultsTbl$FDR < 0.05] <- "UP"
resultsTbl$sigDE[resultsTbl$logFC < -1 & resultsTbl$FDR < 0.05] <- "DOWN"
#Create volcano plot
jpeg("180mV_CTL_plotVolcano.jpg")
ggplot(data=resultsTbl, aes(x=logFC, y=-log10(FDR), color = sigDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = color_subset)
dev.off()

#Perform an exact test for 180mV vs 100mV
tested <- exactTest(list, pair=c("100mV", "180mV"))
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table
#Identify significantly DE genes
resultsTbl$sigDE <- "NA"
resultsTbl$sigDE[resultsTbl$logFC > 1 & resultsTbl$FDR < 0.05] <- "UP"
resultsTbl$sigDE[resultsTbl$logFC < -1 & resultsTbl$FDR < 0.05] <- "DOWN"
#Create volcano plot
jpeg("180mV_100mV_plotVolcano.jpg")
ggplot(data=resultsTbl, aes(x=logFC, y=-log10(FDR), color = sigDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = color_subset)
dev.off()
