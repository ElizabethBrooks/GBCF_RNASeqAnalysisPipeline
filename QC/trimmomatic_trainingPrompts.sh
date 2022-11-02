#!/bin/bash
## job script header that requests 8 threads

# script to perform trimmomatic trimming of paired end reads
# usage: qsub trimmomatic_training.sh readPath outputsPath
# usage Ex: qsub trimmomatic_training.sh /afs/crc.nd.edu/group/genomics/R2D2/220705_Yoon_Adipocyte_Pool2_RNAseq /scratch365/ebrooks5/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq
# usage Ex: qsub trimmomatic_training.sh /afs/crc.nd.edu/group/genomics/R2D2/220707_Yoon_Jurkat_Pool1_RNAseq /scratch365/ebrooks5/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq

# required software for OSCER


# retrieve paired reads absolute path for alignment


# retrieve analysis outputs absolute path


# retrieve adapter absolute path for alignment


# make a new directory for project analysis files


# name of a new directory for outputs of this analysis stage


# make the new outputs directory


# move to the outputs directory


# loop through all forward and reverse reads and run trimmomatic on each pair


	# trim extension from current file name


	# set paired file name


	# trim to sample tag


	# print status message


	# determine phred score for trimming


	# perform adapter trimming on paired reads using 8 threads


# end loop


# print final status message

