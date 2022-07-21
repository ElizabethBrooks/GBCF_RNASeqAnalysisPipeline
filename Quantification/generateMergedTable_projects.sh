#!/bin/bash
#Script to generate guide file and merge gene counts using the merge_tables.py script
#Usage: bash generateMergedTable_projects.sh inputsFile
#Usage ex: bash generateMergedTable_projects.sh inputPaths_yoon_adipocyte_July2022.txt
#Usage ex: bash generateMergedTable_projects.sh inputPaths_yoon_junkrat_July2022.txt

#Load necessary modules for ND CRC servers
module load bio

#Retrieve input argument of a inputs file
inputsFile=$1

#Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")
#Retrieve genome features absolute path for alignment
guideFile=$(grep "guideFile:" ../"InputData/"$inputsFile | tr -d " " | sed "s/guideFile://g")

#Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir

#Retrieve aligned reads input absolute path
inputsPath=$outputsPath"/counted"

#Name output file of inputs
inputOutFile="pipeline_summary.txt"

#Name tmp guide file for merging
tmpGuide=$inputsPath"/tmp_guideFile.txt"

#Loop through all counted paired reads and create a guide file
for f1 in "$inputsPath"/*/; do
	currTag=$(echo $f1 | sed 's/.$//')
	echo "'$f1'counts.txt $currTag" >> $tmpGuide
done

#Move to location of merge_tagles.py script
cd ../util

#Merge gene counts based on generated guide file
python merge_tables.py $tmpGuide
echo "python merge_tables.py $tmpGuide" >> $inputOutFile

#Move the output merged counts file
mv merged_counts.txt $inputsPath

#Clean up
rm $tmpGuide
