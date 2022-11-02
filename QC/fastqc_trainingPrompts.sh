#!/bin/bash
## job script header

# script to perform fastqc quality control of paired end reads
# usage: qsub fastqc_training.sh readPath outputsPath
# usage Ex: qsub fastqc_training.sh /afs/crc.nd.edu/group/genomics/R2D2/220705_Yoon_Adipocyte_Pool2_RNAseq /scratch365/ebrooks5/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq
# usage Ex: qsub fastqc_training.sh /afs/crc.nd.edu/group/genomics/R2D2/220707_Yoon_Jurkat_Pool1_RNAseq /scratch365/ebrooks5/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq

# required software for OSCER servers


# retrieve paired reads absolute path for alignment


# retrieve analysis outputs absolute path


# make the new directory for project analysis files


# name of a new directory for outputs of this analysis stage


# make the outputs directory


# move to the outputs directory


# loop through all forward and reverse reads and run trimmomatic on each pair


	# trim extension from current file name


	# set paired file name


	# print status message


	# perform QC on the first paired end read for the current sample


	# perform QC on the second paired end read for the current sample


# end loop


# print final status message

