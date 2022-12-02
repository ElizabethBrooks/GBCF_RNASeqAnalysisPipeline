#!/usr/bin/env Rscript

# R script to perform statistical analysis of gene count tables using edgeR exact tests
# usage: Rscript edgeR_prompts.r mergedCountsFile
# usage Ex: Rscript edgeR_prompts.r /YOUR/PATH/Adipocyte/counted/merged_counts_formatted.txt
# usage Ex: Rscript edgeR_prompts.r /YOUR/PATH/Jurkat/counted/merged_counts_formatted.txt

# set working directory


# install packages, this should only need to be done once


# load libraries


# connect to a an ensembl website mart by specifying a BioMart and dataset parameters for Mus_musculus.GRCm39


# retrieve gene attributes


# retrieve all hallmark gene sets


# retrieve the apoptosis gene set


# retrieve the TNF alpha signaling gene set


# retrieve the inflammatory response gene set


# retrieve the interferon alpha response gene set


# import gene count data


# subset the input counts


# trim the data table


# set number of samples


# add grouping factor


# create DGE list object


#
## Prep Stage
#

# retrieve library sizes


# plot the library sizes before normalization


# there is no purpose in analyzing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels


# keep only the genes that meet our criteria


# calculate normalized factors


# calculate cpm of normalized counts


# write normalized cpm to file


# draw a MDS plot to show the relative similarities of the samples
# and to view batch and treatment effects after normalization


# draw a heatmap of individual RNA-seq samples using moderated
# log-counts-per-million after normalization


# produce a matrix of pseudo-counts using estimates of the common and tagwise dispersions


# view dispersion estimates and biological coefficient of variation



#
## DE Analysis Stage
#

# setup color blind safe color vectors for plotting


# prepare a gene counts table


# create data frame with the experimental design layout


# attach rownames to the data frame



##
### 100mV vs CTL
##

# perform an exact test for 100mV vs CTL


# create results table of DE genes


# write results table of DE genes to a CSV file


# plot log-fold change against log-counts per million, with DE genes highlighted
# the blue lines indicate 2-fold changes


# identify significantly DE genes


# create volcano plot


# remove the sigDE column


# return all rows from resultsTbl_clean and genes_att


# sort the results by FDR


# write the sorted DE genes to a CSV file


# filter for significantly DE genes


# return all rows from resultsTbl_filt and sample_counts


# return all rows from resultsTbl with a match in h_apoptosis and format for plotting


# return all rows from resultsTbl with a match in h_inflammatory and format for plotting


# return all rows from resultsTbl with a match in h_tnfa and format for plotting


# return all rows from resultsTbl with a match in h_interferon and format for plotting


# create heatmap for resultsTbl_apoptosis_filt


# create heatmap for resultsTbl_inflammatory_filt


# create heatmap for resultsTbl_tnfa_filt


# create heatmap for resultsTbl_interferon_filt



##
### 180mV vs CTL
##

# perform an exact test for 180mV vs CTL


# create results table of DE genes


# write results table of DE genes to a CSV file


# plot log-fold change against log-counts per million, with DE genes highlighted
# the blue lines indicate 2-fold changes


# identify significantly DE genes


# create volcano plot


# remove the sigDE column


# return all rows from resultsTbl_clean and genes_att


# sort the results by FDR


# write the sorted DE genes to a CSV file


# filter for significantly DE genes


# return all rows from resultsTbl_filt and sample_counts


# return all rows from resultsTbl with a match in h_apoptosis and format for plotting


# return all rows from resultsTbl with a match in h_inflammatory and format for plotting


# return all rows from resultsTbl with a match in h_tnfa and format for plotting


# return all rows from resultsTbl with a match in h_interferon and format for plotting


# create heatmap for resultsTbl_apoptosis_filt


# create heatmap for resultsTbl_inflammatory_filt


# create heatmap for resultsTbl_tnfa_filt


# create heatmap for resultsTbl_interferon_filt



##
### 180mV vs 100mV
##

# perform an exact test for 180mV vs 100mV


# create results table of DE genes


# write results table of DE genes to a CSV file


# plot log-fold change against log-counts per million, with DE genes highlighted
# the blue lines indicate 2-fold changes


# identify significantly DE genes


# create volcano plot


# remove the sigDE column


# return all rows from resultsTbl_clean and genes_att


# sort the results by FDR


# write the sorted DE genes to a CSV file

