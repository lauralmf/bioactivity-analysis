#!/bin/bash

# Do this first
mkdir bakta_results
cd bakta_results

for file in ../assembly_files/*.fasta; do bakta --db db "$file"; done

# For each iteration, a bakta output folder is made for a single MAG and saved to the bakta_results folder.
