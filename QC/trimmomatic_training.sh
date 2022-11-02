#!/bin/bash
## job script header that requests 8 threads

# script to perform trimmomatic trimming of paired end reads
# usage: qsub trimmomatic_training.sh readPath outputsPath
# usage Ex: qsub trimmomatic_training.sh /afs/crc.nd.edu/group/genomics/R2D2/220705_Yoon_Adipocyte_Pool2_RNAseq /scratch365/ebrooks5/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq
# usage Ex: qsub trimmomatic_training.sh /afs/crc.nd.edu/group/genomics/R2D2/220707_Yoon_Jurkat_Pool1_RNAseq /scratch365/ebrooks5/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq

#Required software for OSCER
##

# retrieve paired reads absolute path for alignment
readPath="$1"

# retrieve analysis outputs absolute path
outputsPath="$2"

# retrieve adapter absolute path for alignment
##adapterPath=""

# make a new directory for project analysis files
mkdir $outputsPath

# name of a new directory for outputs of this analysis stage
trimOut=$outputsPath"/trimmed"

# make the new outputs directory
mkdir $trimOut

# move to the outputs directory
cd $trimOut

# loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in "$readPath"/*_R1_001.fastq.gz; do

	# trim extension from current file name
	curSample=$(echo $f1 | sed 's/_R._001\.fastq\.gz//')

	# set paired file name
	f2=$curSample"_R2_001.fastq.gz"

	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_R._001\.fastq\.gz//')

	# print status message
	echo "Processing $sampleTag"

	# determine phred score for trimming
	if grep -iF "Illumina 1.5" $outputsPath"/qc/"$sampleTag"_R1_001_fastqc/fastqc_data.txt"; then
		score=64
	elif grep -iF "Illumina 1.9" $outputsPath"/qc/"$sampleTag"_R1_001_fastqc/fastqc_data.txt"; then
		score=33
	else
		echo "ERROR: Illumina encoding not found... exiting"
		exit 1
	fi

	# perform adapter trimming on paired reads using 8 threads
	trimmomatic PE -threads 8 -phred"$score" $f1 $f2 $sampleTag"_pForward.fq.gz" $sampleTag"_uForward.fq.gz" $sampleTag"_pReverse.fq.gz" $sampleTag"_uReverse.fq.gz" ILLUMINACLIP:"$adapterPath" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

	# clean up
	rm -r $noPath"_R1_001_fastqc.zip"
	rm -r $noPath"_R1_001_fastqc/"
	rm -r $noPath"_R2_001_fastqc.zip"
	rm -r $noPath"_R2_001_fastqc/"
done

#Print status message
echo "Analysis complete!"
