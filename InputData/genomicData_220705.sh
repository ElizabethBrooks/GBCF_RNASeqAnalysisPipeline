#!/bin/bash

# BASH script to retrieve genomic data for project 220705
# usage: bash genomicData_220705.sh outputsPath
# usage ex: bash genomicData_220705.sh /scratch365/ebrooks5/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq

# retrieve outputs absolute path
outputsPath="$1"

# move to the outputs directory
cd $outputsPath

# retrieve the genome assembly
wget https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# retrieve the gene annotations
wget https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.gtf.gz
