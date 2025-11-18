# bioactivity-analysis
>-- Under construction --
## Discovery of novel antibiotics

The code in this repository makes up a program that is used to discover the antibiotic potential of metagenomics isolates against pathogenic bacteria.
In brief, the program covers:

- Analysis of bioactivity patterns of metagenomic isolates when screened for bioactivity against other pathogens.
- Genomic profile analysis of metagenomic isolates (metagenomics-assembled genomes - MAGs).
- Correlation between bioactivity screenings performed in lab and genomic profiles.

The program includes several steps, which can be run individually or in a complete workflow (_snakemake file tba_). The description assumes that the reader has some experience working with whole-genome sequencing (WGS) data. 
All scripts are found in the source ("src") folder.

## The raw bioactivity data
The raw bioactivity dataset should be a csv or tsv file that looks something like this if the bioactivity is annotated from triplicate results:

  |  SEQID  | Pathogen 1 | Pathogen 2 | ... | Pathogen N |
  | --------- | ---------- | ---------- | --- | ---------- |
  | SEQID.barcode1 |     111    |     101    | ... |    000     |
  |    ...    |     ...    |     ...    | ... |    ...     |
  | SEQID.barcodeN |     011    |     111    | ... |    111     |

Depending on the purpose, you can do a one-hot or triplicate encoding of the data, where the triplicate encoding counts the number of positive hits, and the one-hot encoding assigns a "1" for at least one positive hit, otherwise a "0". The script for triplicate encoding is included in the ```S1``` script but does not write a file as it currently is, since I found the one-hot encoded data most meaningful (but can easily be added to the script).

## Preprocessing
>[Optional] Filter for low-quality MAGs
Before running the various bioinformatic pipelines that this program includes, I recommend running [CheckM2]([url](https://github.com/chklovski/CheckM2)) on all of the MAGs to only include high-quality MAGs before continuing.

>[If you decided to run CheckM2] To filter for low-quality genomes, the CheckM2 results are added to the raw bioactivity data by left-joining on the "SEQID" column (the MAG identifiers) before filtering. The columns from CheckM2 can now be removed from the dataframe.

The preprocessing steps involves running different bioinformatic pipelines. For automation, I set up a script that does this on multiple genomes at onces. The input is an assembly folder with MAGs as input:
1) Running [antiSMASH]([url](https://docs.antismash.secondarymetabolites.org/install/)) in the command line
2) Processing antiSMASH outputs
3) Running [BiG-SCAPE]([url](https://github.com/medema-group/BiG-SCAPE/wiki/))
4) Running [Prokka]([url](https://github.com/tseemann/prokka)) and [bakta]([url](https://github.com/oschwengers/bakta))* (for whole-genome annotations)
5) Running [GTDB-Tk]([url](https://github.com/Ecogenomics/GTDBTk)) for taxonomic classification
6) Running [BUSCO]([url](https://busco.ezlab.org/busco_userguide.html#installation-with-conda)) for extracting BUSCO genes
7) Constructing a multiple sequence alignment (MSA) from BUSCO genes using [BUSCO_phylogenomics]([url](https://github.com/jamiemcg/BUSCO_phylogenomics))

_*Depending on your purpose and number of genomes, I can recommend running Prokka instead of bakta to save time, but personally, I prefer bakta._

### Data wrangling
>[If you decided to run CheckM2] To filter for low-quality genomes, the CheckM2 results are added to the raw bioactivity data by left-joining on the "SEQID" column (the MAG identifiers) before filtering. The columns from CheckM2 can now be removed from the dataframe.

1) First, the GTDB-Tk annotations from the output tsv files are added to the bioactivity data by left-joining on the "SEQID" column. The GTDB-Tk-assigned taxonomy is nicely separated into Phylum, Genus, and Strain in the ```S1_raw_data_cleanup_add_gtdbtk.R```file, and the dataframe now looks like this:

  |  SEQID  | Pathogen 1 | Pathogen 2 | ... | Pathogen N |  Phylum  |  Genus  |  Strain  |
  | --------- | ---------- | ---------- | --- | ---------- | -------- | ------- | -------- | 
  | SEQID.barcode1 |     111    |     101    | ... |    000     | Phylum 1 | Genus 1 | Strain 1 |
  |    ...    |     ...    |     ...    | ... |    ...     |    ...   |   ...   |    ...   |
  | SEQID.barcodeN |     011    |     111    | ... |    111     | Phylum N | Genus N | Strain N |

>-- Under construction --
