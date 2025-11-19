#!/bin/bash

export GTDBTK_DATA_PATH=path/to/mash_db
gtdbtk classify_wf --genome_dir ./assembly_files/ -x fasta --out_dir ./gtdbtk_results/ --mash_db path/to/mash_db
