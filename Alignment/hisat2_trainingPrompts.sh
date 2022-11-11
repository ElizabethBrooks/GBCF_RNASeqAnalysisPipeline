#!/bin/bash
## job script header that requests 8 threads

#Script to perform hisat2 alignment of trimmed paired end reads
#Usage: sbatch hisat2_trainingPrompts.sh inputsFile outputsFile referenceGenome
#Usage Ex: sbatch hisat2_trainingPrompts.sh /YOUR/PATH/adipocyte/aligned /YOUR/PATH/adipocyte /YOUR/PATH/adipocyte/....fa
#Usage Ex: sbatch hisat2_trainingPrompts.sh /YOUR/PATH/jurkat/aligned /YOUR/PATH/jurkat /YOUR/PATH/jurkat/....fa

# required software for OSCER


# retrieve paired reads absolute path for alignment


# retrieve analysis outputs absolute path


# retrieve adapter absolute path for alignment


# make a new directory for project analysis files


# name of a new directory for outputs of this analysis stage


# make the new outputs directory


# move to the outputs directory


# name of a new directory for outputs of genome build step


# create build output directory for genome reference


# trim path from build file


# trim file extension from build file


# copy genome build fasta file to hisat2 build folder


# trim file extension


# begin hisat2 build


# loop through all forward and reverse paired reads and run Hisat2 on each pair
# using 8 threads and samtools to convert output sam files to bam


	#trim extension from current file name


	# trim file path from current file name


	#  trim extension from current file name


	# create directory for current sample outputs


	# print status message


	# run hisat2 with default settings


	# convert output sam files to bam format for downstream analysis


	# remove the now converted .sam file


	# print status message


# end loop


# clean up and remove build directory


# clean up and remove extra files


# print status message


