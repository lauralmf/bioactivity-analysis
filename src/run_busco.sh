#!/bin/bash

for file in ./*.fasta; do busco -i "$file" -m genome --lineage_dataset bacteria_odb12 -c 8; done