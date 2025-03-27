# Command for running BUSCO_phylogentics from command line
BUSCO_phylogenetics.py -i busco_results -o busco_phylo_results -t 8

# If you experience an error running BUSCO_phylogenetics, run MUSCLE and TrimAL manually
for file in ./*; do muscle -in muscle "$file" -out "${file%.faa}"_alignment.fasta; done
for file in ./*alignment.fasta; do trimal -in "$file" -automated1 -out "${file%.fasta}"_trimmed.fasta; done

# Run AMAS to concatenate trimmed alignments
AMAS.py concat -i ./*.fasta -f fasta -d aa -p concatenated_partition.txt -t concatenated_alignment.fasta

# Run AMAS to create files for tree
AMAS.py convert -i concatenated_alignment.fasta -f fasta -d aa -u phylip # Optional

# Run fasttree to make tree - either from phylip or fasta
fasttree concatenated_alignment.fasta > busco_tree.tree
