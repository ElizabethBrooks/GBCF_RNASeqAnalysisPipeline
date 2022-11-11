#!/bin/bash
## job script header that requests 1 thread

# script to generate guide file and merge gene counts using the merge_tables.py script
# usage: bash generateMergedTable_trainingPrompts.sh inputsFile outputsFile
# usage ex: bash generateMergedTable_trainingPrompts.sh /YOUR/PATH/Adipocyte/counted /YOUR/PATH/Adipocyte
# usage ex: bash generateMergedTable_trainingPrompts.sh /YOUR/PATH/Jurkat/counted /YOUR/PATH/Jurkat

# required software for OSCER


# retrieve paired reads absolute path for alignment


# retrieve analysis outputs absolute path


# retrieve reference genome absolute path


# make a new directory for project analysis files


# name of a new directory for outputs of this analysis stage


# make the new outputs directory


# move to the outputs directory


# name tmp guide file for merging


# loop through all counted paired reads


	# trim extension from current file


	# create a guide file


# end loop


# move to location of merge_tagles.py script in the util director


# merge gene counts using the generated guide file and merge_tables.py


# trim / from the inputs path


# remove file paths from sample tags


# clean up sample names for downstream analysis with R


# add the header to the final re-formatted merged counts file


# add the count data to the final re-formatted merged counts file


# clean up and remove the tmp guide file


# clean up and remove the un-formatted merged counts file


