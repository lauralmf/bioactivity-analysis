#!/bin/bash

# Make output folder and change directory to output folder.
mkdir -p busco_results && cd busco_results

for file in ../assembly_files/*.fasta; do busco -i "$file" -m genome --lineage_dataset bacteria_odb12 -c 8; done

# You might need to do 
# for file in ./*.fasta; do busco -i "$file" -m genome --lineage_dataset bacteria_odb12 -c 8; done
# from in the assembly_files folder and then copy the results to a busco_results folder from the assembly_files folder.
