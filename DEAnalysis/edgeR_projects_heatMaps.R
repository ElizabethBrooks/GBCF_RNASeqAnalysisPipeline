#!/usr/bin/env Rscript
#R script to perform statistical analysis of gene count tables using edgeR exact test
#Usage: Rscript edgeR_projects.r mergedCountsFile
#Usage Ex: Rscript edgeR_projects.r 220705_Yoon_Adipocyte_Pool2_RNAseq_merged_counts_formatted.txt
#Usage Ex: Rscript edgeR_projects.r 220707_Yoon_Jurkat_Pool1_RNAseq_merged_counts_formatted.txt

#Install BiocManager, this should only need to be done once
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

#Retrieve gene attributes
genes_att <- getBM(attributes=c('ensembl_gene_id',
                                'external_gene_name',
                                #'uniprot_gn_id',
                                'description'), 
                    mart = ensembl)

#Retrieve all hallmark gene sets
h_gene_sets = msigdbr(species = "mouse", category = "H")

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
#DE
setwd("/Users/bamflappy/GBCF/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq/differential_expression/annotated")
#setwd("/Users/bamflappy/GBCF/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq/differential_expression/annotated")

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
#Calculate normalized counts
normList <- cpm(list, normalized.lib.sizes=TRUE)

#Retrieve log-counts-per-million after normalization
logcpm <- cpm(list, log=TRUE)

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)


#DEA Stage

#Setup color blind appropriate vector for plotting
color_subset <- c("#0000FF", "#000000", "#FF0000")

#Prepare a gene counts table
sample_counts <- data.frame(list$counts) %>% rownames_to_column(var="gene")

#Create data frame with the experimental design layout
exp_factor <- data.frame(sample = c("100mV", "100mV", "100mV", "180mV", "180mV", "180mV", "CTL", "CTL", "CTL"))
rownames(exp_factor) <- colnames(list$counts)

#Perform an exact test for 100mV vs CTL
tested <- exactTest(list, pair=c("CTL", "100mV"))
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table %>% rownames_to_column(var="gene")
#Filter for significantly DE genes
resultsTbl_filt <- resultsTbl[resultsTbl$FDR < 0.05,]
#Return all rows from genes_att with a match in resultsTbl_filt
resultsTbl_att <- genes_att %>% semi_join(resultsTbl_filt, by = c("ensembl_gene_id" = "gene"))
#Return all rows from resultsTbl_att and sample_counts
resultsTbl_counts <- resultsTbl_att %>% inner_join(sample_counts, by = c("ensembl_gene_id" = "gene"))
#Return all rows from resultsTbl with a match in h_apoptosis and format for plotting
resultsTbl_apoptosis <- resultsTbl_counts %>% semi_join(h_apoptosis, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_apoptosis_filt <- resultsTbl_apoptosis[ , -which(names(resultsTbl_apoptosis) %in% c("ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_apoptosis_filt) <- resultsTbl_apoptosis$external_gene_name
#Return all rows from resultsTbl with a match in h_inflammatory and format for plotting
resultsTbl_inflammatory <- resultsTbl_counts %>% semi_join(h_inflammatory, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_inflammatory_filt <- resultsTbl_inflammatory[ , -which(names(resultsTbl_inflammatory) %in% c("ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_inflammatory_filt) <- resultsTbl_inflammatory$external_gene_name
#Return all rows from resultsTbl with a match in h_tnfa and format for plotting
resultsTbl_tnfa <- resultsTbl_counts %>% semi_join(h_tnfa, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_tnfa_filt <- resultsTbl_tnfa[ , -which(names(resultsTbl_tnfa) %in% c("ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_tnfa_filt) <- resultsTbl_tnfa$external_gene_name
#Return all rows from resultsTbl with a match in h_interferon and format for plotting
resultsTbl_interferon <- resultsTbl_counts %>% semi_join(h_interferon, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_interferon_filt <- resultsTbl_interferon[ , -which(names(resultsTbl_interferon) %in% c("ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_interferon_filt) <- resultsTbl_interferon$external_gene_name
#Create heatmap for resultsTbl_apoptosis_filt
as.ggplot(pheatmap(resultsTbl_apoptosis_filt, scale="row", annotation_col = exp_factor,
         main="100mV vs CTL Apoptosis Response"))
ggsave("100mV_CTL_heatmap_apoptosis.png", bg = "white")
#Create heatmap for resultsTbl_inflammatory_filt
as.ggplot(pheatmap(resultsTbl_inflammatory_filt, scale="row", annotation_col = exp_factor,
         main="100mV vs CTL Inflammatory Response"))
ggsave("100mV_CTL_heatmap_inflammatory.png", bg = "white")
#Create heatmap for resultsTbl_tnfa_filt
as.ggplot(pheatmap(resultsTbl_tnfa_filt, scale="row", annotation_col = exp_factor,
         main="100mV vs CTL TNF Alpha Response"))
ggsave("100mV_CTL_heatmap_tnfa.png", bg = "white")
#Create heatmap for resultsTbl_interferon_filt
as.ggplot(pheatmap(resultsTbl_interferon_filt, scale="row", annotation_col = exp_factor,
         main="100mV vs CTL Interferon Alpha Response"))
ggsave("100mV_CTL_heatmap_interferon.png", bg = "white")

#Perform an exact test for 180mV vs CTL
tested <- exactTest(list, pair=c("CTL", "180mV"))
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table %>% rownames_to_column(var="gene")
#Filter for significantly DE genes
resultsTbl_filt <- resultsTbl[resultsTbl$FDR < 0.05,]
#Return all rows from genes_att with a match in resultsTbl_filt
resultsTbl_att <- genes_att %>% semi_join(resultsTbl_filt, by = c("ensembl_gene_id" = "gene"))
#Return all rows from resultsTbl_att and sample_counts
resultsTbl_counts <- resultsTbl_att %>% inner_join(sample_counts, by = c("ensembl_gene_id" = "gene"))
#Return all rows from resultsTbl with a match in h_apoptosis and format for plotting
resultsTbl_apoptosis <- resultsTbl_counts %>% semi_join(h_apoptosis, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_apoptosis_filt <- resultsTbl_apoptosis[ , -which(names(resultsTbl_apoptosis) %in% c("ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_apoptosis_filt) <- resultsTbl_apoptosis$external_gene_name
#Return all rows from resultsTbl with a match in h_inflammatory and format for plotting
resultsTbl_inflammatory <- resultsTbl_counts %>% semi_join(h_inflammatory, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_inflammatory_filt <- resultsTbl_inflammatory[ , -which(names(resultsTbl_inflammatory) %in% c("ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_inflammatory_filt) <- resultsTbl_inflammatory$external_gene_name
#Return all rows from resultsTbl with a match in h_tnfa and format for plotting
resultsTbl_tnfa <- resultsTbl_counts %>% semi_join(h_tnfa, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_tnfa_filt <- resultsTbl_tnfa[ , -which(names(resultsTbl_tnfa) %in% c("ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_tnfa_filt) <- resultsTbl_tnfa$external_gene_name
#Return all rows from resultsTbl with a match in h_interferon and format for plotting
resultsTbl_interferon <- resultsTbl_counts %>% semi_join(h_interferon, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_interferon_filt <- resultsTbl_interferon[ , -which(names(resultsTbl_interferon) %in% c("ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_interferon_filt) <- resultsTbl_interferon$external_gene_name
#Create heatmap for resultsTbl_apoptosis_filt
as.ggplot(pheatmap(resultsTbl_apoptosis_filt, scale="row", annotation_col = exp_factor,
                   main="180mV vs CTL Apoptosis Response"))
ggsave("180mV_CTL_heatmap_apoptosis.png", bg = "white")
#Create heatmap for resultsTbl_inflammatory_filt
as.ggplot(pheatmap(resultsTbl_inflammatory_filt, scale="row", annotation_col = exp_factor,
                   main="180mV vs CTL Inflammatory Response"))
ggsave("180mV_CTL_heatmap_inflammatory.png", bg = "white")
#Create heatmap for resultsTbl_tnfa_filt
as.ggplot(pheatmap(resultsTbl_tnfa_filt, scale="row", annotation_col = exp_factor,
                   main="180mV vs CTL TNF Alpha Response"))
ggsave("180mV_CTL_heatmap_tnfa.png", bg = "white")
#Create heatmap for resultsTbl_interferon_filt
as.ggplot(pheatmap(resultsTbl_interferon_filt, scale="row", annotation_col = exp_factor,
                   main="180mV vs CTL Interferon Alpha Response"))
ggsave("180mV_CTL_heatmap_interferon.png", bg = "white")

#Perform an exact test for 180mV vs 100mV
#tested <- exactTest(list, pair=c("100mV", "180mV"))
#Create results table of DE genes
#resultsTbl <- topTags(tested, n=nrow(tested$table))$table %>% rownames_to_column(var="gene")
#No significantly DE genes
