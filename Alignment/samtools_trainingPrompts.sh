#!/bin/bash
## job script header that requests 8 threads

# script to perform samtools sorting of trimmed and aligned paired end reads
# usage: sbatch samtools_trainingPrompts.sh inputsFile outputsFile
# usage Ex: sbatch samtools_trainingPrompts.sh /YOUR/PATH/adipocyte/aligned /YOUR/PATH/adipocyte
# usage Ex: sbatch samtools_trainingPrompts.sh /YOUR/PATH/jurkat/aligned /YOUR/PATH/jurkat

# required software for OSCER


# retrieve paired reads absolute path for alignment


# retrieve analysis outputs absolute path


# make a new directory for project analysis files


# name of a new directory for outputs of this analysis stage


# make the new outputs directory


# move to the outputs directory


# loop through all reads and sort sam/bam files for input to samtools


	# name of aligned file


	# trim file path from current folder name


	# create directory for current sample outputs


	# print status message


	# run samtools to prepare mapped reads for sorting using 8 threads


	# run samtools fixmate to update paired-end flags for singletons
	

	# remove the current sortedName.bam file


	# remove duplicate reads using samtools markdup


	# print status message


# end loop


# clean up and remove extra files


# print final status message


