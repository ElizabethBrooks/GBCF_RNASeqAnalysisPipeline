#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N htseq_projects_jobOutput
#Script to run htseq-count on trimmed, aligned, and name sorted paired end reads
#Usage Ex: qsub htseq_projects.sh inputPaths_yoon_adipocyte_July2022.txt
#Usage Ex: qsub htseq_projects.sh inputPaths_yoon_junkrat_July2022.txt

#Required modules for ND CRC servers
module load bio

#Retrieve input argument of a inputs file
inputsFile=$1

#Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")
#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeFeatures:" ../"InputData/"$inputsFile | tr -d " " | sed "s/genomeFeatures://g")

#Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir

#Retrieve aligned reads input absolute path
inputsPath=$outputsPath"/sorted"

#Make an outputs directory for analysis
anOut=$outputsPath"/counted"
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
htseq-count --help | tail -1 >> $versionFile

#Loop through all sorted forward and reverse paired reads and store the file locations in an array
for f1 in "$inputsPath"/*/; do
	#Name of aligned file
	curAlignedSample=$f1"accepted_hits.bam"
	#Trim file path from current folder name
	curSampleNoPath=$(basename $f1)
	#Create directory for current sample outputs
	mkdir $curSampleNoPath
	#Output current sample name to summary file
	echo $curSampleNoPath >> $inputOutFile
	#Output status message
	echo "Sample $curSampleNoPath is being counted..."
	#Count reads using htseq-count
	htseq-count -f bam -s no -t gene $curAlignedSample $genomeFile > $curSampleNoPath"/counts.txt"
	#Add run inputs to output summary file
	echo "$curSampleNoPath" >> "$inputOutFile"
	echo "htseq-count -f bam -s no -t gene "$curAlignedSample" "$genomeFile" > "$curSampleNoPath"/counts.txt" >> $inputOutFile
	#Output status message
	echo "Sample $curSampleNoPath has been counted!"
done

#Clean up
rm -r $inputsPath

#Print status message
echo "Analysis complete!"

