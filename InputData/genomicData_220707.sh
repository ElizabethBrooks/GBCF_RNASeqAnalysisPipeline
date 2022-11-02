#!/bin/bash

# BASH script to retrieve genomic data for project 220707
# usage: bash genomicData_220707.sh outputsPath
# usage ex: bash genomicData_220707.sh /scratch365/ebrooks5/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq

# retrieve outputs absolute path
outputsPath="$1"

# move to the outputs directory
cd $outputsPath

# retrieve the genome assembly
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# retrieve the gene annotations
wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz
