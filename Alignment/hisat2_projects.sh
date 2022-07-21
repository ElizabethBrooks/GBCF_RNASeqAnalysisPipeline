#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N hisat2_projects_jobOutput
#$ -pe smp 8
#Script to perform hisat2 alignment of trimmed paired end reads
#Usage: qsub hisat2_projects.sh inputsFile
#Usage Ex: qsub hisat2_projects.sh inputPaths_yoon_adipocyte_July2022.txt
#Usage Ex: qsub hisat2_projects.sh inputPaths_yoon_junkrat_July2022.txt

#Required modules for ND CRC servers
module load bio

#Retrieve input argument of a inputs file
inputsFile=$1

#Retrieve genome reference absolute path for alignment
ref=$(grep "genomeReference:" ../"InputData/"$inputsFile | tr -d " " | sed "s/genomeReference://g")
#Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

#Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir

#Make an outputs directory for analysis
anOut=$outputsPath"/aligned"
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
hisat2 --version >> $versionFile
samtools --version >> $versionFile

#Set trimmed reads absolute path
inputsPath=$outputsPath"/trimmed"

#Create build output directory for Hisat reference
buildOut="build"
mkdir $buildOut
#Trim path and file extension from build file
refNoPath=$(basename $ref)
refNoEx=$(echo $refNoPath | sed 's/\.fasta//' | sed 's/\.fa//')
#Copy genome build fasta file to hisat2 build folder
cp "$ref" "$buildOut"/"$refNoEx"
#Trim file extension
refNoPath=$(echo $refNoEx | sed 's/\.fa//g')
#Begin hisat2 build
echo "Beginning hisat2 build... "
hisat2-build -p 8 -f "$buildOut"/"$refNoEx" "$buildOut"/"$refNoPath"
echo hisat2-build -p 8 -f "$buildOut"/"$refNoEx" "$buildOut"/"$refNoPath" >> $inputOutFile
echo "hisat2 build complete!"

#Loop through all forward and reverse paired reads and run Hisat2 on each pair
# using 8 threads and samtools to convert output sam files to bam
for f1 in "$inputsPath"/*pForward.fq.gz; do
	#Trim extension from current file name
	curSample=$(echo $f1 | sed 's/.pForward\.fq\.gz//')
	#Trim file path from current file name
	curSampleNoPath=$(basename $f1)
	curSampleNoPath=$(echo $curSampleNoPath | sed 's/.pForward\.fq\.gz//')
	#Create directory for current sample outputs
	mkdir "$curSampleNoPath"
	#Print status message
	echo "Processing $curSampleNoPath"
	#Run hisat2 with default settings
	hisat2 -p 8 -q -x $buildOut"/"$refNoEx -1 $f1 -2 $curSample"_pReverse.fq.gz" -S $curSampleNoPath"/accepted_hits.sam" --summary-file $curSampleNoPath"/alignedSummary.txt"
	#Add sample and hisat2 run inputs to output summary file
	echo $curSampleNoPath >> $inputOutFile
	echo "hisat2 -p 8 -q -x "$buildOut"/"$refNoEx" -1 "$f1" -2 "$curSample"_pReverse.fq.gz -S "$curSampleNoPath"/accepted_hits.sam --summary-file "$curSampleNoPath"/alignedSummary.txt" >> "$inputOutFile"
	#Convert output sam files to bam format for downstream analysis
	samtools view -@ 8 -bS $curSampleNoPath"/accepted_hits.sam" > $curSampleNoPath"/accepted_hits.bam"
	#Remove the now converted .sam file
	rm $curSampleNoPath"/accepted_hits.sam"
	#Add samtools run inputs to output summary file
	echo "samtools view -@ 8 -bS "$curSampleNoPath"/accepted_hits.sam > "$curSampleNoPath"/accepted_hits.bam" >> $inputOutFile
	#Print status message
	echo "Processed!"
done

#Clean up
rm -r "build"
rm -r $inputsPath

#Print status message
echo "Analysis complete!"
