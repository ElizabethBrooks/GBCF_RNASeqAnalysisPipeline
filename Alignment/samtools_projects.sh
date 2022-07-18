#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N samtools_projects_jobOutput
#$ -pe smp 8
#Script to perform samtools sorting of trimmed and aligned paired end reads
#Usage: qsub samtools_projects.sh inputsFile
#Usage Ex: qsub samtools_projects.sh inputPaths_yoon_adipocyte_July2022.txt name
#Usage Ex: qsub samtools_projects.sh inputPaths_yoon_junkrat_July2022.txt name

#Required modules for ND CRC servers
module load bio

#Retrieve input argument of a inputs file
inputsFile=$1

#Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

#Retrieve sorting method flags from input
if [[ "$2" == "name" || "$2" == "Name" || "$2" == "n" || "$2" == "N" ]]; then
	#Name sorted flag with num threads flag
	flags="-@ 8 -n"
	methodTag="Name"
elif [[ "$2" == "coordinate" || "$2" == "Coordinate" || "$2" == "c" || "$2" == "C" ]]; then
	#Coordinate sorted with num threads flag
	flags="-@ 8"
	methodTag="Coordinate"
else
	#Report error with input flag
	echo "ERROR: a flag for sorting method (name or coordinate) is expected... exiting"
	exit 1
fi

#Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir

#Retrieve aligned reads input absolute path
inputsPath=$outputsPath"/aligned"

#Make an outputs directory for analysis
anOut=$outputsPath"/sorted"
mkdir $anOut
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $anOut directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to the outputs directory
cd $anOut

#Name output file of inputs
inputOutFile="summary.txt"
#Add software versions to outputs
samtools --version > $inputOutFile

#Loop through all reads and sort sam/bam files for input to samtools
for f1 in "$inputsPath"/*/; do
	#Name of aligned file
	curAlignedSample="$f1"accepted_hits.bam
	#Trim file path from current folder name
	curSampleNoPath=$(basename "$f1")
	#Create directory for current sample outputs
	mkdir "$anOut"/"$curSampleNoPath"
	#Output current sample name to summary file
	echo "$curSampleNoPath" >> $inputOutFile
	#Run samtools to prepare mapped reads for sorting
	#using 8 threads
	echo "Sample $curSampleNoPath is being name sorted..."
	samtools sort -@ 8 -n -o "$anOut"/"$curSampleNoPath"/sortedName.bam -T /tmp/"$curSampleNoPath".sortedName.bam "$curAlignedSample"
	echo "Sample $curSampleNoPath has been name sorted!"
	#Add run inputs to output summary file
	echo samtools sort -@ 8 -n -o "$anOut"/"$curSampleNoPath"/sortedName.bam -T /tmp/"$curSampleNoPath".sortedName.bam "$curAlignedSample" >> "$inputOutFile"
	#Determine which sorting method is to be performed
	if [[ "$methodTag" == "Coordinate" ]]; then
		#Run fixmate -m to update paired-end flags for singletons
		echo "Sample $curSampleNoPath singleton flags are being updated..."
		samtools fixmate -m "$anOut"/"$curSampleNoPath"/sortedName.bam "$anOut"/"$curSampleNoPath"/sortedFixed.bam
		echo "Sample $curSampleNoPath singleton flags have been updated!"
		#Clean up
		rm "$anOut"/"$curSampleNoPath"/sortedName.bam
		#Run samtools to prepare mapped reads for sorting by coordinate
		#using 8 threads
		echo "Sample $curSampleNoPath is being sorted..."
		samtools sort "$flags" -o "$anOut"/"$curSampleNoPath"/accepted_hits.bam -T /tmp/"$curSampleNoPath".sorted.bam "$anOut"/"$curSampleNoPath"/sortedFixed.bam
		echo "Sample $curSampleNoPath has been sorted!"
		rm "$anOut"/"$curSampleNoPath"/sortedFixed.bam
		#Remove duplicate reads
		samtools markdup -r "$anOut"/"$curSampleNoPath"/accepted_hits.bam "$anOut"/"$curSampleNoPath"/noDups.bam
		#Add run inputs to output summary file
		echo samtools fixmate -m "$anOut"/"$curSampleNoPath"/sortedName.bam "$anOut"/"$curSampleNoPath"/sortedFixed.bam >> "$inputOutFile"
		echo samtools sort "$flags" -o "$anOut"/"$curSampleNoPath"/accepted_hits.bam -T /tmp/"$curSampleNoPath".sorted.bam "$anOut"/"$curSampleNoPath"/sortedFixed.bam >> "$inputOutFile"
		echo samtools markdup -r "$anOut"/"$curSampleNoPath"/accepted_hits.bam "$anOut"/"$curSampleNoPath"/noDups.bam >> "$inputOutFile"
	else
		#Run fixmate -m to update paired-end flags for singletons
		echo "Sample $curSampleNoPath singleton flags are being updated..."
		samtools fixmate -m "$anOut"/"$curSampleNoPath"/sortedName.bam "$anOut"/"$curSampleNoPath"/accepted_hits.bam
		echo "Sample $curSampleNoPath singleton flags have been updated!"
		rm "$anOut"/"$curSampleNoPath"/sortedName.bam
		#Remove duplicate reads
		samtools markdup -r "$anOut"/"$curSampleNoPath"/accepted_hits.bam "$anOut"/"$curSampleNoPath"/noDups.bam
		#Add run inputs to output summary file
		echo samtools fixmate -m "$anOut"/"$curSampleNoPath"/sortedName.bam "$anOut"/"$curSampleNoPath"/accepted_hits.bam >> "$inputOutFile"
		echo samtools markdup -r "$anOut"/"$curSampleNoPath"/accepted_hits.bam "$anOut"/"$curSampleNoPath"/noDups.bam >> "$inputOutFile"
	fi
done

