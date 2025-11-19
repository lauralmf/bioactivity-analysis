#!/bin/bash

# Make folder where the antiSMASH results should be saved and change directory to the results folder.
mkdir -p antismash_results && cd antismash_results

# Run antiSMASH on each MAG:
for file in ../assembly_files/*.fasta; do
    base=$(basename "$file" .fasta)  # Extract filename without extension
    antismash "$file" --output-dir "../antismash_results/$base" \
        -t bacteria \
        --genefinding-tool prodigal-m --asf \
        --cc-mibig --cb-knownclusters --pfam2go --rre --tfbs --cpu 16
done

# For each iteration, an antiSMASH output folder is made for a MAG and added to the folder, antismash_results.
