# bioactivity-analysis

## Discovery of novel antibiotics

- Analysis of bioactivity patterns of metagenomic isolates when screened for bioactivity against other pathogens.
- Genomic profile analysis of metagenomic isolates.
- Correlation between bioactivity screenings performed in lab and genomic profiles.

This program includes several steps, which can be run individually or in a complete workflow.

### Preprocessing scripts:
- The preprocessing part takes an assembly folder with whole-genome sequenced bacterial genomes as input:
1) Running [antiSMASH]([url](https://docs.antismash.secondarymetabolites.org/install/)) in the command line
2) Processing antiSMASH outputs
3) Running [BiG-SCAPE]([url](https://github.com/medema-group/BiG-SCAPE/wiki/))
4) Running [Prokka]([url](https://github.com/tseemann/prokka)) and [bakta]([url](https://github.com/oschwengers/bakta))* (for whole-genome annotations)
5) Running [GTDB-Tk]([url](https://github.com/Ecogenomics/GTDBTk)) for taxonomic classification
6) Running [BUSCO]([url](https://busco.ezlab.org/busco_userguide.html#installation-with-conda)) for extracting BUSCO genes
7) Constructing a multiple sequence alignment (MSA) from BUSCO genes using [BUSCO_phylogenomics]([url](https://github.com/jamiemcg/BUSCO_phylogenomics))

_*Depending on your purpose and number of genomes, I can recommend running Prokka instead of bakta to save time._

### Exploratory data analysis - part 1:
- The exploratory data analysis requires a csv or tsv file in the following format:
  
  |  Isolate  | Pathogen 1 | Pathogen 2 | ... | Pathogen N |  Phylum  |  Genus  |  Strain  |
  | --------- | ---------- | ---------- | --- | ---------- | -------- | ------- | -------- | 
  | Isolate 1 |     111    |     101    | ... |    000     | Phylum 1 | Genus 1 | Strain 1 |
  |    ...    |     ...    |     ...    | ... |    ...     |    ...   |   ...   |    ...   |
  | Isolate N |     011    |     111    | ... |    111     | Phylum N | Genus N | Strain N |

- Where the cells between the isolates and pathogens are results from triplicate bioactivity screens from co-cultures with pathogens.
- "1" is a positive hit, indicating bioactivity, and "0" is a negative hit, indicating no bioactivity.
- "Phylum", "Genus" and "Strain" are obtained from the GTDB-Tk results.

2) Calculating Z-scores Z-scores from bioactivity matrix

### Exploratory data analysis - part 2:
1) Constructing a phylogenetic tree from assembly 16S sequences

To be continued...
