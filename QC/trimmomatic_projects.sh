#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N trimmomatic_projects_jobOutput
#$ -pe smp 8
#Script to perform trimmomatic trimming of paired end reads
#Usage: qsub trimmomatic_projects.sh inputsFile
#Usage Ex: qsub trimmomatic_projects.sh inputPaths_yoon_adipocyte_July2022.txt
#Usage Ex: qsub trimmomatic_projects.sh inputPaths_yoon_junkrat_July2022.txt

#Required modules for ND CRC servers
module load bio

#Retrieve input argument of a inputs file
inputsFile=$1

#Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")
#Retrieve adapter absolute path for alignment
adapterPath=$(grep "adapter:" ../"InputData/"$inputsFile | tr -d " " | sed "s/adapter://g")
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

#Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir
mkdir $outputsPath

#Make a new directory for analysis
trimOut=$outputsPath"/trimmed"
mkdir $trimOut
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $trimOut directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to the new directory
cd $trimOut

#Name output file of inputs
inputOutFile=$outputsPath"/pipeline_summary.txt"
versionFile=$outputsPath"/version_summary.txt"

#Add software version to outputs
echo "Trimmomatic:" >> $versionFile
trimmomatic -version >> $versionFile

#Loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in "$readPath"/*_R1_001.fastq.gz; do
	#Trim extension from current file name
	curSample=$(echo $f1 | sed 's/_R._001\.fastq\.gz//')
	#Set paired file name
	f2=$curSample"_R2_001.fastq.gz"
	#Trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_R._001\.fastq\.gz//')
	#Print status message
	echo "Processing $sampleTag"
	#Determine phred score for trimming
	if grep -iF "Illumina 1.5" $outputsPath"/qc/"$sampleTag"_R1_001_fastqc/fastqc_data.txt"; then
		score=64
	elif grep -iF "Illumina 1.9" $outputsPath"/qc/"$sampleTag"_R1_001_fastqc/fastqc_data.txt"; then
		score=33
	else
		echo "ERROR: Illumina encoding not found... exiting"
		#echo "ERROR: Illumina encoding not found for $curSample" >> $inputOutFile
		exit 1
	fi
	#Perform adapter trimming on paired reads
	#using 8 threads
	#removed HEADCROP:13
	trimmomatic PE -threads 8 -phred"$score" $f1 $f2 $sampleTag"_pForward.fq.gz" $sampleTag"_uForward.fq.gz" $sampleTag"_pReverse.fq.gz" $sampleTag"_uReverse.fq.gz" ILLUMINACLIP:"$adapterPath" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	#Add run inputs to output summary file
	echo trimmomatic PE -threads 8 -phred"$score" $f1 $f2 $sampleTag"_pForward.fq.gz" $sampleTag"_uForward.fq.gz" $sampleTag"_pReverse.fq.gz" $sampleTag"_uReverse.fq.gz" ILLUMINACLIP:"$adapterPath" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >> $inputOutFile
	#Clean up
	rm -r $noPath"_R1_001_fastqc.zip"
	rm -r $noPath"_R1_001_fastqc/"
	rm -r $noPath"_R2_001_fastqc.zip"
	rm -r $noPath"_R2_001_fastqc/"
	#Print status message
	echo "Processed!"
done

#Print status message
echo "Analysis complete!"
