#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N htseq_projects_jobOutput
#Script to run htseq-count on trimmed, aligned, then name sorted paired end reads


#Required modules for ND CRC servers
module load bio

#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")

	#Retrieve sorted reads input absolute path
	inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
	inputsDir="$inputsPath"/"$1"
	outputsPath="$inputsPath"


#Loop through all sorted forward and reverse paired reads and store the file locations in an array
for f1 in "$inputsDir"/*/*.bam; do
	#Name of sorted and aligned file
	curAlignedSample="$f1"
	#Trim file paths from current sample folder name
	curSampleNoPath=$(echo $f1 | sed 's/accepted\_hits\.bam//g')
	curSampleNoPath=$(basename $curSampleNoPath)
	#Create directory for current sample outputs
	mkdir "$outputFolder"/"$curSampleNoPath"
	#Count reads using htseq-count
	echo "Sample $curSampleNoPath is being counted..."
	#Determine which flags to use based on sorting method
	if [[ "$1"  == sortedName* ]]; then
		#Use name sorted flag (default)
		#https://github.com/simon-anders/htseq/issues/37
		#--secondary-alignments ignore --supplementary-alignments ignore
		#Flag to output features in sam format
		#-o "$outputFolder"/"$curSampleNoPath"/counted.sam
		htseq-count -f bam -a 60 -s no -m union -t gene -i ID "$curAlignedSample" "$genomeFile" > "$outputFolder"/"$curSampleNoPath"/counts.txt
		#Add run inputs to output summary file
		echo "$curSampleNoPath" >> "$inputOutFile"
		echo "htseq-count -f bam -a 60 -s no -m union -t gene -i ID" "$curAlignedSample" "$genomeFile" ">" "$outputFolder"/"$curSampleNoPath"/counts.txt >> "$inputOutFile"
	elif [[ "$1"  == sortedCoordinate* ]]; then
		#Use coordinate sorted flag
		#https://github.com/simon-anders/htseq/issues/37
		#--secondary-alignments ignore --supplementary-alignments ignore
		#Flag to output features in sam format
		#-o "$outputFolder"/"$curSampleNoPath"/counted.sam
		htseq-count -f bam -a 60 -r pos -s no -m union -t gene -i ID "$curAlignedSample" "$genomeFile" > "$outputFolder"/"$curSampleNoPath"/counts.txt
		#Add run inputs to output summary file
		echo "$curSampleNoPath" >> "$inputOutFile"
		echo "htseq-count -f bam -a 60 -r pos -s no -m union -t gene -i ID" "$curAlignedSample" "$genomeFile" ">" "$outputFolder""/""$curSampleNoPath""/counts.txt" >> "$inputOutFile"
	else
		echo "ERROR: The bam file "$f1" was not found... exiting"
		exit 1
	fi
	echo "Sample $curSampleNoPath has been counted!"
done
#Copy previous summaries
cp "$inputsDir"/*.txt "$outputFolder"
