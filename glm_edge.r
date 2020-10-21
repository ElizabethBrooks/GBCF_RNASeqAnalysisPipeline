#!/usr/bin/env Rscript
#Usage: Rscript glmQLF_edgeR.r countsFile startColumn endColumn factorGroupingFile
#Usage Ex: Rscript glmQLF_edgeR.r cleaned.csv 1 24 expDesign_Olympics_GRP1.csv
#R script to perform statistical analysis of gene count tables using edgeR GLM

#Install edgeR and statmod, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("statmod")

#Load the edgeR library
library("edgeR")
library("statmod")

#Import gene count data
countsTable <- read.csv(file=args[1], row.names="gene")[ ,args[2]:args[3]]
#Import grouping factor
targets <- read.csv(file=args[4], row.names="sample")

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$genotype,sep="."))
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
colnames(list) <- targets$sample

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
summary(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases
# between libraries
list <- calcNormFactors(list)
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)

#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

#Next, the NB dispersion is estimated
list <- estimateDisp(list, design, robust=TRUE)

#Now, estimate and plot the QL dispersions
fit <- glmQLFit(list, design, robust=TRUE)

#Test whether the average across all UV groups is equal to the average across
#all VIS groups, to examine the overall effect of treatment
con.UVvsVIS <- makeContrasts(UVvsVIS = (UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
  - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4,
  levels=design)
#Look at genes with significant expression across all UV groups
treat.anov.UVVIS <- glmTreat(fit, contrast=con.UVvsVIS, lfc=log2(1.2))
#Generate table of DE genes
tagsTblANOVA.filtered <- topTags(treat.anov.UVVIS, n=nrow(treat.anov.UVVIS$table))$table
#Keep significantly DE genes
tagsTblANOVA.filtered.keep <- tagsTblANOVA.filtered$FDR <= 0.05

#Test whether the average across all tolerant groups is equal to the average across
#all not tolerant groups, to examine the overall effect of tolerance
con.TvsN <- makeContrasts(TvsN = (UV.Y05 + VIS.Y05 + UV.E05 + VIS.E05)/4
  - (UV.Y023 + VIS.Y023 + UV.R2 + VIS.R2)/4,
  levels=design)
#Look at genes with significant expression across all UV groups
treat.anov.TN <- glmTreat(fit, contrast=con.TvsN, lfc=log2(1.2))
#Generate table of DE genes
tagsTblANOVATN.filtered <- topTags(treat.anov.TN, n=nrow(treat.anov.TN$table))$table
#Keep significantly DE genes
tagsTblANOVATN.filtered.keep <- tagsTblANOVATN.filtered$FDR <= 0.05

#Test whether there is an interaction effect
con.Inter <- makeContrasts(Inter = ((UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
  - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4)
  - ((UV.Y05 + VIS.Y05 + UV.E05 + VIS.E05)/4
  - (UV.Y023 + VIS.Y023 + UV.R2 + VIS.R2)/4),
  levels=design)
#Look at genes with significant expression
treat.anov.Inter <- glmTreat(fit, contrast=con.Inter, lfc=log2(1.2))
#Generate table of DE genes
tagsTblANOVAInter.filtered <- topTags(treat.anov.Inter, n=nrow(treat.anov.Inter$table))$table
#Keep significantly DE genes
tagsTblANOVAInter.filtered.keep <- tagsTblANOVAInter.filtered$FDR <= 0.05