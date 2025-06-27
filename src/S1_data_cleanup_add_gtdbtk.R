library(tidyverse)

# ---- Data cleanup ---- #
# Filter for hiqh-quality sequenced genomes
raw_data <- read_csv("./datasets/raw_data.csv")
raw_data_filtered <- raw_data %>% drop_na()

# Remove irrelevant columns:
raw_data_filtered <- raw_data_filtered %>% 
  select(-matches("Coverage|N50|Contamination|Genome|Complete|Mbp"))

# Replace no-break space characters in ESKAPE column names with space:
clean_names <- function(x) gsub("\u00A0", " ", x) 
colnames(raw_data_filtered) <- clean_names(colnames(raw_data_filtered))

# ---- Add GTDB-Tk annotations ---- #
gtdb <- read_tsv("../gtdb_results/gtdbtk.bac120.summary.tsv") # Read GTDB-Tk outputs
gtdb$user_genome <- gsub("_assembly$", "", gtdb$user_genome) # Remove "_assembly" after SEQID for a more streamlined SEQID

for (i in 1:nrow(raw_data_filtered)){
  for (j in 1:nrow(gtdb)){
    if (raw_data_filtered$SEQID[i] == gtdb$user_genome[j]){
      raw_data_filtered$Classification[i] <- gtdb$classification[j]
      break
    }
  }
  if (is.na(raw_data_filtered$Classification[i])) {
    raw_data_filtered$Classification[i] <- raw_data_filtered$Classification[i]
  }
}

# Add separate Phylum, Genus and Strain columns
raw_data_filtered <- raw_data_filtered %>%
  mutate(
    Phylum = ifelse(
      str_detect(Classification, "p__"),
      str_extract(Classification, "p__[^;\\s]+") %>% str_remove("p__"),
      NA),
    
    Genus = ifelse(
      str_detect(Classification, "g__"),
      str_extract(Classification, "g__[^;\\s]+") %>% str_remove("g__") %>% str_extract("^[^_]+"), # Extract before "_"
      ifelse(str_detect(Classification, " sp\\.$"), 
             str_extract(Classification, "^[^\\s]+"), 
             NA)),
    
    Strain = ifelse(
      str_detect(Classification, "s__"),
      str_extract(Classification, "s__[^;]+") %>% str_remove("s__"),
      paste0(Genus, " sp.")
    )
  )

# For those that do not have a specific species identifier, which in my case were:
strain_phylum_map <- data.frame(
  Strain = c("Arthrobacter sp.", "Streptomyces sp.", "Peribacillus sp.", "Lysinibacillus sp.", "Psychrobacillus sp.", "Luteibacter sp.",
            "Kitasatospora sp.", "Curtobacterium sp.", "Variovorax sp.", "Pseudomonas_E sp."),
  Phylum = c("Actinomycetota", "Actinomycetota", "Bacillota", "Bacillota", "Bacillota", "Pseudomonadota",
             "Actinomycetota", "Actinomycetota", "Pseudomonadota", "Pseudomonadota")
)

raw_data_filtered <- raw_data_filtered %>%
  mutate(Phylum = case_when(
    !is.na(Phylum) ~ Phylum,
    Strain %in% strain_phylum_map$Strain ~ strain_phylum_map$Phylum[match(Strain, strain_phylum_map$Strain)], TRUE ~ NA_character_))

# Use Pseudomonas as genus instead of Pseudomonas_E
raw_data_filtered <- raw_data_filtered %>%
  mutate(Genus = ifelse(Genus == "Pseudomonas_E", "Pseudomonas", Genus))

# Remove GTDB-Tk-formatted classification column
raw_data_filtered <- raw_data_filtered %>% select(-Classification)

raw_data_filtered <- raw_data_filtered %>% 
  mutate(Strain_full = paste0(Strain, " (", SEQID, ")"))

write_csv(raw_data_filtered, "./datasets/clean_data.csv") # Write to new csv file

# ---- One-hot encoding of the data ---- #
eskape_cols <- raw_data_filtered %>% select(matches("DSM|AU")) %>% colnames()
encoded <- raw_data_filtered %>% mutate(across(all_of(eskape_cols), ~ str_count(.x, "1")))

write_csv(encoded, "./datasets/clean_encoded_bioactivity.csv")

# ---- Triplicate encoding of the data ---- #
triplicate <- raw_data_filtered

for (col in colnames(raw_data_filtered)) {
  if (!col %in% c("SEQID", "Medium", "Genus", "Strain", "Strain_full", "Phylum")) {
    l <- numeric(nrow(raw_data_filtered))  
    
    for (i in 1:nrow(raw_data_filtered)) {
      l[i] <- str_count(string = raw_data_filtered[i, col], "1")  # Count number of occurrences of "1"
    }
    triplicate[[col]] <- l  
  }
}

write_csv(triplicate, "./datasets/clean_triplicate_bioactivity.csv")
