#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N fastqc_yoon_jobOutput
#$ -pe smp 8
#Script to perform fastqc quality control of paired end reads
#Usage: qsub fastqc_yoon.sh
#Usage Ex: qsub fastqc_yoon.sh

#Required modules for ND CRC servers
module load bio
#Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../InputData/inputPaths_yoon.txt | tr -d " " | sed "s/pairedReads://g")
#Retrieve adapter absolute path for alignment
adapterPath=$(grep "adapter:" ../InputData/inputPaths_yoon.txt | tr -d " " | sed "s/adapter://g")
#Retrieve trimming outputs absolute path
outputsPath=$(grep "qc:" ../InputData/outputPaths_yoon.txt | tr -d " " | sed "s/qc://g")
#Move to outputs directory
cd $(dirname $outputsPath)
#Make a new directory for trimming
qcOut=$outputsPath
mkdir $qcOut
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $qcOut directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to the new directory
cd $qcOut
#Report software version
fastqc -version > "summary.txt"
#Loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in "$readPath"/*_R1_001.fastq.gz; do
	#Trim extension from current file name
	curSample=$(echo $f1 | sed 's/._001\.fastq\.gz//')
	#Set paired file name
	#f2=$curSample"2_001.fastq.gz"
	#Print status message
	echo "Processing $curSample"
	#Perform QC on both paired end reads for the current sample
	#fastqc $f1 -o $qcOut --extract
	#fastqc $f2 -o $qcOut --extract
	#Quickly check the QC results of the first pair
	#if grep -iF "WARN" $curSample"1_001_fastqc/summary.txt"; then
	#	grep -iF "WARN" $curSample"1_001_fastqc/summary.txt" > $curSample"1_001_fastqc_report.txt"
	#fi
	#if grep -iF "FAIL" $curSample"1_001_fastqc/summary.txt"; then
	#	grep -iF "FAIL" $curSample"1_001_fastqc/summary.txt" > $curSample"1_001_fastqc_report.txt"
	#fi
	#Quickly check the QC results of the second pair
	#if grep -iF "WARN" $curSample"2_001_fastqc/summary.txt"; then
	#	grep -iF "WARN" $curSample"2_001_fastqc/summary.txt" > $curSample"2_001_fastqc_report.txt"
	#fi
	#if grep -iF "FAIL" $curSample"2_001_fastqc/summary.txt"; then
	#	grep -iF "FAIL" $curSample"2_001_fastqc/summary.txt" > $curSample"2_001_fastqc_report.txt"
	#fi
	echo "Processed!"
done
#Print status message
echo "Analysis complete!"
