#!/bin/bash

# Create folder for bakta results if not present and change directory to the results folder.
mkdir -p bakta_results && cd bakta_results

# Run script
for file in ../assembly_files/*.fasta; do bakta --db db "$file"; done

# For each iteration, a bakta output folder is made for a single MAG and saved to the bakta_results folder.
