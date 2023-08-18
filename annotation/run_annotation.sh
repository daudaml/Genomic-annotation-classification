#!/bin/bash

# RUN PYTHON DATA DOWNLOAD

python3 data_collection.py


# RUN PRODIGAL ON FNA FILES

if [ -d prodigal_results ]; then

    rm -r prodigal_results

fi



mkdir prodigal_results

mkdir prodigal_results/biocontrol/

mkdir prodigal_results/pathogen/



find . -maxdepth 2 -name "*.fna" | while read file; do

    prodigal -i "$file" -o prodigal_results/${file%.fna}.gbk -a prodigal_results/${file%.fna}.faa

done



# REMOVE ASTERIX FROM OUTPUT FILES

find . -maxdepth 3 -name "*.faa" | while read file; do

    sed 's/*//' "$file" > "${file%.faa}_no_asterix.faa"

done


# RUN INTERPROSCAN ON FAA FILES

if [ -d interproscan_results ]; then

    rm -r interproscan_results

fi



mkdir interproscan_results

mkdir interproscan_results/biocontrol/

mkdir interproscan_results/pathogen/



find . -maxdepth 3 -name "*_no_asterix.faa" | while read file; do

    interproscan.sh -i "$file" -f tsv -b interproscan_results/$(dirname "$file" | sed 's/prodigal_results\///')/$(basename "$file" _no_asterix.faa) -appl Pfam,TIGRFAM

done
