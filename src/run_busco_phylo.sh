#!/bin/bash

# Command for running BUSCO_phylogentics from command line
BUSCO_phylogenetics.py -i ../busco_results -o ../busco_phylo_results -t 8

# If you experience an error running BUSCO_phylogenetics, run MUSCLE and TrimAL manually
for file in ../busco_phylo_results/supermatrix/sequences/*; do muscle -in muscle "$file" -out "${file%.faa}"_alignment.fasta; done
for file in ../busco_phylo_results/supermatrix/*alignment.fasta; do trimal -in "$file" -automated1 -out "${file%.fasta}"_trimmed.fasta; done

# Run AMAS to concatenate trimmed alignments
AMAS.py concat -i ../assembly_files/*.fasta -f fasta -d aa -p ../busco_phylo_results/supermatrix/concatenated_partition.txt -t ../busco_phylo_results/supermatrix/concatenated_alignment.fasta

# Run AMAS to create files for tree
AMAS.py convert -i ../busco_phylo_results/supermatrix/concatenated_alignment.fasta -f fasta -d aa -u phylip # Optional

# Run fasttree to make tree - either from phylip or fasta
fasttree ../busco_phylo_results/supermatrix/concatenated_alignment.fasta > ../busco_phylo_results/supermatrix/busco_tree.tree
