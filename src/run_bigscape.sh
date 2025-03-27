#!/bin/bash

# Collect all BGC region files outputted by antiSMASH in all folders.
find ../antismash_results -iname *region*.gbk -exec sh -c 'cp $1 $2' sh {} ../gbk_regions \; 2>/dev/null

# Follow instructions from https://github.com/medema-group/BiG-SCAPE/wiki/, then:
bigscape cluster -i ../gbk_regions/ -o ../bigscape_results/ -p ../BiG-SCAPE/big_scape/Pfam-A.hmm -m 3.1 --include-singletons -c 8 --mix --extend-strategy simple_match
