#!/bin/bash
#Script to run merge gene IDs from DE analysis results
#Usage: bash mergeDEResults.sh mergedGeneIDsFile
#Usage Ex: bash mergeDEResults.sh mergedGeneIDs

#Loop over DE analysis results
header="START"
for f in *_DEResults.csv; do
	#Determine genotype based on column numbers
	if [[ $f == "1_6"* ]]; then #Y05
		genotype="Y05"
	elif [[ $f == "7_12"* ]]; then #Y023
		genotype="Y023"
	elif [[ $f == "13_18"* ]]; then #E05
		genotype="E05"
	elif [[ $f == "19_24"* ]]; then #R2
		genotype="R2"
	elif [[ $f == "25_30"* ]]; then #PA
		genotype="PA"
	elif [[ $f == "31_36"* ]]; then #Sierra
		genotype="Sierra"
	fi
	#Retrieve first column without header
	# and write to a temporary file
	cut -d"," -f1 $f | tail -n+2 > "$genotype"_geneIDs.csv
	#Add genotype to header
	header="$header,$genotype"
done

#Add gene IDs to file with header
echo $header | sed 's/START,//g' > mergedGeneIDs.csv
#Merge gene Ids across genotypes
paste -d "," *_geneIDs.csv >> mergedGeneIDs.csv
#Clean up
rm *_geneIDs.csv