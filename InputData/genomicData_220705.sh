#!/bin/bash

# BASH script to retrieve genomic data for project 220705
# usage: bash genomicData_220705.sh

# retrieve the genome assembly
wget https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# retrieve the gene annotations
wget https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.104.gtf.gz
