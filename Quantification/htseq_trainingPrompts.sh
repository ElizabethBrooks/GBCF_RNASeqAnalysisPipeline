#!/bin/bash
## job script header that requests 1 thread

# script to run htseq-count on trimmed, aligned, and name sorted paired end reads
# usage: sbatch htseq_trainingPrompts.sh inputsFile outputsFile genomeFeatures
# usage Ex: sbatch htseq_trainingPrompts.sh /YOUR/PATH/Adipocyte/sorted /YOUR/PATH/Adipocyte /YOUR/PATH/Adipocyte/....gtf
# usage Ex: sbatch htseq_trainingPrompts.sh /YOUR/PATH/Jurkat/sorted /YOUR/PATH/Jurkat /YOUR/PATH/Jurkat/....gtf

# required software for OSCER


# retrieve paired reads absolute path for alignment


# retrieve analysis outputs absolute path


# retrieve reference genome absolute path


# make a new directory for project analysis files


# name of a new directory for outputs of this analysis stage


# make the new outputs directory


# move to the outputs directory


# loop through all sorted forward and reverse paired reads and store the file locations in an array


	# name of aligned file


	# trim file path from current folder name


	# create directory for current sample outputs


	# output current sample name to summary file


	# output status message


	# count reads using htseq-count


	# output status message


# end loop


#clean up


# print status message


