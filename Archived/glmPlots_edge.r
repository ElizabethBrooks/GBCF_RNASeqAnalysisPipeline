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
#head(countsTable)
#Import grouping factor
targets <- read.csv(file=args[4], row.names="sample")

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$genotype,sep="."))
#cbind(targets,Group=group)
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
colnames(list) <- targets$sample

#Plot the library sizes before normalization
jpeg("glmQLF_plotBarsBefore.jpg")
barplot(list$samples$lib.size*1e-6, names=1:24, ylab="Library size (millions)")
dev.off()

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
summary(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases
# between libraries
list <- calcNormFactors(list)
#list$samples
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)
write.table(normList, file="glmQLF_normalizedCounts.csv", sep=",", row.names=TRUE)

#Verify TMM normalization using a MD plot
#Write plot to file
jpeg("glmQLF_plotMDBefore.jpg")
plotMD(cpm(list, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)
dev.off()

#Use a MDS plot to visualizes the differences
# between the expression profiles of different samples
points <- c(0,1,2,3,15,16,17,18)
colors <- rep(c("blue", "darkgreen", "red", "black"), 2)
#Write plot with legend to file
jpeg("glmQLF_plotMDS.jpg")
plotMDS(list, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)
dev.off()
#Write plot without legend to file
jpeg("glmQLF_plotMDS_noLegend.jpg")
plotMDS(list, col=colors[group], pch=points[group])
dev.off()

#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
#design

#Next, the NB dispersion is estimated
list <- estimateDisp(list, design, robust=TRUE)
#list$common.dispersion
#Visualize the dispersion estimates with a BCV plot
#Write plot to file
jpeg("glmQLF_plotBCV.jpg")
plotBCV(list)
dev.off()

#Now, estimate and plot the QL dispersions
fit <- glmQLFit(list, design, robust=TRUE)
#head(fit$coefficients)
#Write plot to file
jpeg("glmQLF_plotQLDisp.jpg")
plotQLDisp(fit)
dev.off()

#Test whether the average across all UV groups is equal to the average across
#all VIS groups, to examine the overall effect of treatment
con.UVvsVIS <- makeContrasts(UVvsVIS = (UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
  - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4,
  levels=design)
#Look at genes with significant expression across all UV groups
treat.anov.UVVIS <- glmTreat(fit, contrast=con.UVvsVIS, lfc=log2(1.2))
summary(decideTests(treat.anov.UVVIS))
#Write plot to file
jpeg("glmQLF_2WayANOVA_UVvsVIS_plotMD_filtered.jpg")
plotMD(treat.anov.UVVIS)
abline(h=c(-1, 1), col="blue")
dev.off()
#Write tags table of DE genes to file
tagsTblANOVA.filtered <- topTags(treat.anov.UVVIS, n=nrow(treat.anov.UVVIS$table))$table
tagsTblANOVA.filtered.keep <- tagsTblANOVA.filtered$FDR <= 0.05
tagsTblANOVA.filtered.out <- tagsTblANOVA.filtered[tagsTblANOVA.filtered.keep,]
write.table(tagsTblANOVA.filtered.out, file="glmQLF_2WayANOVA_UVvsVIS_topTags_filtered.csv", sep=",", row.names=TRUE)

#Test whether the average across all tolerant groups is equal to the average across
#all not tolerant groups, to examine the overall effect of tolerance
con.TvsN <- makeContrasts(TvsN = (UV.Y05 + VIS.Y05 + UV.E05 + VIS.E05)/4
  - (UV.Y023 + VIS.Y023 + UV.R2 + VIS.R2)/4,
  levels=design)
#Look at genes with significant expression across all UV groups
treat.anov.TN <- glmTreat(fit, contrast=con.TvsN, lfc=log2(1.2))
summary(decideTests(treat.anov.TN))
#Write plot to file
jpeg("glmQLF_2WayANOVA_TvsN_plotMD_filtered.jpg")
plotMD(treat.anov.TN)
abline(h=c(-1, 1), col="blue")
dev.off()
#Write tags table of DE genes to file
tagsTblANOVATN.filtered <- topTags(treat.anov.TN, n=nrow(treat.anov.TN$table))$table
tagsTblANOVATN.filtered.keep <- tagsTblANOVATN.filtered$FDR <= 0.05
tagsTblANOVATN.filtered.out <- tagsTblANOVATN.filtered[tagsTblANOVATN.filtered.keep,]
write.table(tagsTblANOVATN.filtered.out, file="glmQLF_2WayANOVA_TvsN_topTags_filtered.csv", sep=",", row.names=TRUE)

#Test whether there is an interaction effect
con.Inter <- makeContrasts(Inter = ((UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
  - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4)
  - ((UV.Y05 + VIS.Y05 + UV.E05 + VIS.E05)/4
  - (UV.Y023 + VIS.Y023 + UV.R2 + VIS.R2)/4),
  levels=design)
#Look at genes with significant expression
treat.anov.Inter <- glmTreat(fit, contrast=con.Inter, lfc=log2(1.2))
summary(decideTests(treat.anov.Inter))
#Write plot to file
jpeg("glmQLF_2WayANOVA_interaction_plotMD_filtered.jpg")
plotMD(treat.anov.Inter)
abline(h=c(-1, 1), col="blue")
dev.off()
#Generate table of DE genes
tagsTblANOVAInter.filtered <- topTags(treat.anov.Inter, n=nrow(treat.anov.Inter$table))$table
tagsTblANOVAInter.filtered.keep <- tagsTblANOVAInter.filtered$FDR <= 0.05
tagsTblANOVAInter.filtered.out <- tagsTblANOVAInter.filtered[tagsTblANOVAInter.filtered.keep,]
write.table(tagsTblANOVAInter.filtered.out, file="glmQLF_2WayANOVA_interaction_topTags_filtered.csv", sep=",", row.names=TRUE)
