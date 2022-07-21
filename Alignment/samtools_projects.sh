#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N samtools_projects_jobOutput
#$ -pe smp 8
#Script to perform samtools sorting of trimmed and aligned paired end reads
#Usage: qsub samtools_projects.sh inputsFile
#Usage Ex: qsub samtools_projects.sh inputPaths_yoon_adipocyte_July2022.txt
#Usage Ex: qsub samtools_projects.sh inputPaths_yoon_junkrat_July2022.txt

#Required modules for ND CRC servers
module load bio

#Retrieve input argument of a inputs file
inputsFile=$1

#Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

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
inputOutFile=$outputsPath"/pipeline_summary.txt"
versionFile=$outputsPath"/version_summary.txt"
#Add software versions to outputs
samtools --version >> $versionFile

#Loop through all reads and sort sam/bam files for input to samtools
for f1 in "$inputsPath"/*/; do
	#Name of aligned file
	curAlignedSample=$f1"accepted_hits.bam"
	#Trim file path from current folder name
	curSampleNoPath=$(basename $f1)
	#Create directory for current sample outputs
	mkdir $curSampleNoPath
	#Output current sample name to summary file
	echo $curSampleNoPath >> $inputOutFile
	#Run samtools to prepare mapped reads for sorting
	#using 8 threads
	echo "Sample $curSampleNoPath is being name sorted..."
	samtools sort -@ 8 -n -o $curSampleNoPath"/sortedName.bam" -T "/tmp/"$curSampleNoPath".sortedName.bam" $curAlignedSample
	echo "Sample $curSampleNoPath has been name sorted!"
	#Add run inputs to output summary file
	echo "samtools sort -@ 8 -n -o "$curSampleNoPath"/sortedName.bam -T /tmp/"$curSampleNoPath".sortedName.bam "$curAlignedSample >> $inputOutFile
	#Run fixmate -m to update paired-end flags for singletons
	samtools fixmate -m $curSampleNoPath"/sortedName.bam" $curSampleNoPath"/accepted_hits.bam"
	rm "$curSampleNoPath"/sortedName.bam
	#Remove duplicate reads
	samtools markdup -r $curSampleNoPath"/accepted_hits.bam" $curSampleNoPath"/noDups.bam"
	#Add run inputs to output summary file
	echo "samtools fixmate -m "$curSampleNoPath"/sortedName.bam "$curSampleNoPath"/accepted_hits.bam" >> $inputOutFile
	echo "samtools markdup -r "$curSampleNoPath"/accepted_hits.bam "$curSampleNoPath"/noDups.bam" >> $inputOutFile
	echo "Sample $curSampleNoPath has been name sorted!"
done

#Clean up
rm -r $inputsPath

#Print status message
echo "Analysis complete!"
