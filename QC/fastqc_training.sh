#!/bin/bash
## job script header

# script to perform fastqc quality control of paired end reads
# usage: qsub fastqc_training.sh readPath outputsPath
# usage Ex: qsub fastqc_training.sh /afs/crc.nd.edu/group/genomics/R2D2/220705_Yoon_Adipocyte_Pool2_RNAseq /scratch365/ebrooks5/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq
# usage Ex: qsub fastqc_training.sh /afs/crc.nd.edu/group/genomics/R2D2/220707_Yoon_Jurkat_Pool1_RNAseq /scratch365/ebrooks5/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq

# required software for OSCER servers
##

# retrieve paired reads absolute path for alignment
readPath="$1"

# retrieve analysis outputs absolute path
outputsPath="$2"

# make the new directory for project analysis files
mkdir $outputsPath

# name of a new directory for outputs of this analysis stage
qcOut=$outputsPath"/qc"

# make the outputs directory
mkdir $qcOut

# move to the outputs directory
cd $qcOut

# loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in "$readPath"/*_R1_001.fastq.gz; do

	# trim extension from current file name
	curSample=$(echo $f1 | sed 's/_R._001\.fastq\.gz//')

	# set paired file name
	f2=$curSample"_R2_001.fastq.gz"

	# print start status message
	echo "Processing $curSample"

	# perform QC on both paired end reads for the current sample
	fastqc $f1 -o $qcOut --extract
	fastqc $f2 -o $qcOut --extract
done

# print final status message
echo "Analysis complete!"
