#!/bin/bash

for file in ../assembly_files/*.fasta; do
    base=$(basename "$file" .fasta)  # Extract filename without extension
    antismash "$file" --output-dir "../antismash_results/$base" \
        -t bacteria \
        --genefinding-tool prodigal-m --asf \
        --cc-mibig --cb-knownclusters --pfam2go --rre --tfbs --cpu 16
done