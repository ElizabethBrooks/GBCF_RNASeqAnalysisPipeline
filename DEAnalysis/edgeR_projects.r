#!/usr/bin/env Rscript
#R script to perform statistical analysis of gene count tables using edgeR exact test
#Usage: Rscript edgeR_projects.r mergedCountsFile
#Usage Ex: Rscript edgeR_projects.r 220705_Yoon_Adipocyte_Pool2_RNAseq_merged_counts_formatted.txt
#Usage Ex: Rscript edgeR_projects.r 220707_Yoon_Jurkat_Pool1_RNAseq_merged_counts_formatted.txt

#Install edgeR, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#Install packages, this should only need to be done once
#BiocManager::install("edgeR")
#BiocManager::install("biomaRt")
#install.packages("msigdbr")
#install.packages("ggplot2")
#install.packages("tibble")
#install.packages("dplyr")
#install.packages("pheatmap")
#install.packages("ggplotify")

#Load libraries
library(edgeR)
library(biomaRt)
library(msigdbr)
library(ggplot2)
library(tibble)
library(dplyr)
library(pheatmap)
library(ggplotify)

#Connect to a an ensembl website mart by specifying a BioMart and dataset parameters
#Mus_musculus.GRCm39
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
#Homo_sapiens.GRCh38
#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

#View available attributes
#View(listAttributes(ensembl))

#Retrieve gene attributes
genes_att <- getBM(attributes=c('ensembl_gene_id',
                                'external_gene_name',
                                #'uniprot_gn_id',
                                'description'), 
                   mart = ensembl)

#Retrieve all hallmark gene sets
h_gene_sets = msigdbr(species = "mouse", category = "H")
#h_gene_sets = msigdbr(species = "human", category = "H")

#Retrieve the apoptosis gene set
h_apoptosis <- h_gene_sets %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS")

#Retrieve the TNF alpha signaling gene set
h_tnfa <- h_gene_sets %>%
  dplyr::filter(gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

#Retrieve the inflammatory response gene set
h_inflammatory <- h_gene_sets %>%
  dplyr::filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE")

#Retrieve the interferon alpha response gene set
h_interferon <- h_gene_sets %>%
  dplyr::filter(gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE")

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
#DE Adipocyte
setwd("/Users/bamflappy/GBCF/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq/differential_expression")
#DE Jurkat
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
#Full set
#plotMDS(list)
#No undetermined
#plotMDS(list, col=rep(1:3, each=4))
#DE
plotMDS(list, col=rep(1:3, each=3))
dev.off()
#Draw a heatmap of individual RNA-seq samples using moderated
# log-counts-per-million after normalization
logcpm <- cpm(list, log=TRUE)
jpeg("plotHeatMap_cpm_afterNormalize.jpg")
heatmap(logcpm)
dev.off()

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)
#View dispersion estimates and biological coefficient of variation
jpeg("plotBCV.jpg")
plotBCV(list)
dev.off()


#DEA Stage

#Setup color vector for plotting
color_subset <- c("#0000FF", "#000000", "#FF0000")

#Prepare a gene counts table
sample_counts <- data.frame(logcpm) %>% rownames_to_column(var="gene")

#Create data frame with the experimental design layout
exp_factor <- data.frame(group)
rownames(exp_factor) <- colnames(logcpm)

#Perform an exact test for 100mV vs CTL
tested <- exactTest(list, pair=c("CTL", "100mV"))
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table %>% rownames_to_column(var="gene")
write.csv(resultsTbl, file="100mV_CTL_topTags.csv")
#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
jpeg("100mV_CTL_plotMD.jpg")
plotMD(tested)
abline(h=c(-1, 1), col="grey")
dev.off()
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
#Remove the sigDE column
resultsTbl_clean <- resultsTbl[ , -which(names(resultsTbl) %in% c("sigDE"))]
#Return all rows from resultsTbl_clean and genes_att
resultsTbl_att <- resultsTbl_clean %>% left_join(genes_att, by = c("gene" = "ensembl_gene_id"))
#Sort the results by FDR and write to file
resultsTbl_ann <- resultsTbl_att[order(resultsTbl_att$FDR),]
write.csv(resultsTbl_ann, file="annotated/100mV_CTL_annotated.csv")
#Filter for significantly DE genes
resultsTbl_filt <- resultsTbl_ann[resultsTbl_ann$FDR < 0.05,]
#Return all rows from resultsTbl_filt and sample_counts
resultsTbl_counts <- resultsTbl_filt %>% inner_join(sample_counts, by = "gene")
#Return all rows from resultsTbl with a match in h_apoptosis and format for plotting
resultsTbl_apoptosis <- resultsTbl_counts %>% semi_join(h_apoptosis, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_apoptosis_filt <- resultsTbl_apoptosis[ , -which(names(resultsTbl_apoptosis) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_apoptosis_filt) <- resultsTbl_apoptosis$external_gene_name
resultsTbl_apoptosis_filt <- head(resultsTbl_apoptosis_filt, n = 20)
#Return all rows from resultsTbl with a match in h_inflammatory and format for plotting
resultsTbl_inflammatory <- resultsTbl_counts %>% semi_join(h_inflammatory, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_inflammatory_filt <- resultsTbl_inflammatory[ , -which(names(resultsTbl_inflammatory) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_inflammatory_filt) <- resultsTbl_inflammatory$external_gene_name
resultsTbl_inflammatory_filt <- head(resultsTbl_inflammatory_filt, n = 20)
#Return all rows from resultsTbl with a match in h_tnfa and format for plotting
resultsTbl_tnfa <- resultsTbl_counts %>% semi_join(h_tnfa, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_tnfa_filt <- resultsTbl_tnfa[ , -which(names(resultsTbl_tnfa) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_tnfa_filt) <- resultsTbl_tnfa$external_gene_name
resultsTbl_tnfa_filt <- head(resultsTbl_tnfa_filt, n = 20)
#Return all rows from resultsTbl with a match in h_interferon and format for plotting
resultsTbl_interferon <- resultsTbl_counts %>% semi_join(h_interferon, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_interferon_filt <- resultsTbl_interferon[ , -which(names(resultsTbl_interferon) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_interferon_filt) <- resultsTbl_interferon$external_gene_name
resultsTbl_interferon_filt <- head(resultsTbl_interferon_filt, n = 20)
#Create heatmap for resultsTbl_apoptosis_filt
as.ggplot(pheatmap(resultsTbl_apoptosis_filt, scale="row", annotation_col = exp_factor,
                   main="100mV vs CTL Apoptosis Response"))
ggsave("annotated/100mV_CTL_heatmap_cpm_apoptosis.png", bg = "white")
#Create heatmap for resultsTbl_inflammatory_filt
as.ggplot(pheatmap(resultsTbl_inflammatory_filt, scale="row", annotation_col = exp_factor,
                   main="100mV vs CTL Inflammatory Response"))
ggsave("annotated/100mV_CTL_heatmap_cpm_inflammatory.png", bg = "white")
#Create heatmap for resultsTbl_tnfa_filt
as.ggplot(pheatmap(resultsTbl_tnfa_filt, scale="row", annotation_col = exp_factor,
                   main="100mV vs CTL TNF Alpha Response"))
ggsave("annotated/100mV_CTL_heatmap_cpm_tnfa.png", bg = "white")
#Create heatmap for resultsTbl_interferon_filt
as.ggplot(pheatmap(resultsTbl_interferon_filt, scale="row", annotation_col = exp_factor,
                   main="100mV vs CTL Interferon Alpha Response"))
ggsave("annotated/100mV_CTL_heatmap_cpm_interferon.png", bg = "white")

#Perform an exact test for 180mV vs CTL
tested <- exactTest(list, pair=c("CTL", "180mV"))
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table %>% rownames_to_column(var="gene")
write.csv(resultsTbl, file="180mV_CTL_topTags.csv")
#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
jpeg("180mV_CTL_plotMD.jpg")
plotMD(tested)
abline(h=c(-1, 1), col="grey")
dev.off()
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
#Remove the sigDE column
resultsTbl_clean <- resultsTbl[ , -which(names(resultsTbl) %in% c("sigDE"))]
#Return all rows from resultsTbl_clean and genes_att
resultsTbl_att <- resultsTbl_clean %>% left_join(genes_att, by = c("gene" = "ensembl_gene_id"))
#Sort the results by FDR and write to file
resultsTbl_ann <- resultsTbl_att[order(resultsTbl_att$FDR),]
write.csv(resultsTbl_ann, file="annotated/180mV_CTL_annotated.csv")
#Filter for significantly DE genes
resultsTbl_filt <- resultsTbl_ann[resultsTbl_ann$FDR < 0.05,]
#Return all rows from resultsTbl_filt and sample_counts
resultsTbl_counts <- resultsTbl_filt %>% inner_join(sample_counts, by = "gene")
#Return all rows from resultsTbl with a match in h_apoptosis and format for plotting
resultsTbl_apoptosis <- resultsTbl_counts %>% semi_join(h_apoptosis, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_apoptosis_filt <- resultsTbl_apoptosis[ , -which(names(resultsTbl_apoptosis) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_apoptosis_filt) <- resultsTbl_apoptosis$external_gene_name
resultsTbl_apoptosis_filt <- head(resultsTbl_apoptosis_filt, n = 20)
#Return all rows from resultsTbl with a match in h_inflammatory and format for plotting
resultsTbl_inflammatory <- resultsTbl_counts %>% semi_join(h_inflammatory, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_inflammatory_filt <- resultsTbl_inflammatory[ , -which(names(resultsTbl_inflammatory) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_inflammatory_filt) <- resultsTbl_inflammatory$external_gene_name
resultsTbl_inflammatory_filt <- head(resultsTbl_inflammatory_filt, n = 20)
#Return all rows from resultsTbl with a match in h_tnfa and format for plotting
resultsTbl_tnfa <- resultsTbl_counts %>% semi_join(h_tnfa, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_tnfa_filt <- resultsTbl_tnfa[ , -which(names(resultsTbl_tnfa) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_tnfa_filt) <- resultsTbl_tnfa$external_gene_name
resultsTbl_tnfa_filt <- head(resultsTbl_tnfa_filt, n = 20)
#Return all rows from resultsTbl with a match in h_interferon and format for plotting
resultsTbl_interferon <- resultsTbl_counts %>% semi_join(h_interferon, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_interferon_filt <- resultsTbl_interferon[ , -which(names(resultsTbl_interferon) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_interferon_filt) <- resultsTbl_interferon$external_gene_name
resultsTbl_interferon_filt <- head(resultsTbl_interferon_filt, n = 20)
#Create heatmap for resultsTbl_apoptosis_filt
as.ggplot(pheatmap(resultsTbl_apoptosis_filt, scale="row", annotation_col = exp_factor,
                   main="180mV vs CTL Apoptosis Response"))
ggsave("annotated/180mV_CTL_heatmap_cpm_apoptosis.png", bg = "white")
#Create heatmap for resultsTbl_inflammatory_filt
as.ggplot(pheatmap(resultsTbl_inflammatory_filt, scale="row", annotation_col = exp_factor,
                   main="180mV vs CTL Inflammatory Response"))
ggsave("annotated/180mV_CTL_heatmap_cpm_inflammatory.png", bg = "white")
#Create heatmap for resultsTbl_tnfa_filt
as.ggplot(pheatmap(resultsTbl_tnfa_filt, scale="row", annotation_col = exp_factor,
                   main="180mV vs CTL TNF Alpha Response"))
ggsave("annotated/180mV_CTL_heatmap_cpm_tnfa.png", bg = "white")
#Create heatmap for resultsTbl_interferon_filt
as.ggplot(pheatmap(resultsTbl_interferon_filt, scale="row", annotation_col = exp_factor,
                   main="180mV vs CTL Interferon Alpha Response"))
ggsave("annotated/180mV_CTL_heatmap_cpm_interferon.png", bg = "white")

#Perform an exact test for 180mV vs 100mV
tested <- exactTest(list, pair=c("100mV", "180mV"))
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table %>% rownames_to_column(var="gene")
write.csv(resultsTbl, file="180mV_100mV_topTags.csv")
#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
jpeg("180mV_100mV_plotMD.jpg")
plotMD(tested)
abline(h=c(-1, 1), col="grey")
dev.off()
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
#Remove the sigDE column
resultsTbl_clean <- resultsTbl[ , -which(names(resultsTbl) %in% c("sigDE"))]
#Return all rows from resultsTbl_clean and genes_att
resultsTbl_att <- resultsTbl_clean %>% left_join(genes_att, by = c("gene" = "ensembl_gene_id"))
#Sort the results by FDR and write to file
resultsTbl_ann <- resultsTbl_att[order(resultsTbl_att$FDR),]
write.csv(resultsTbl_ann, file="annotated/180mV_100mV_annotated.csv")
