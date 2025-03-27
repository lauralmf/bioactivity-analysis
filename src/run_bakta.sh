#!/bin/bash

for file in ../assembly_files/*.fasta; do bakta --db db "$file"; done
