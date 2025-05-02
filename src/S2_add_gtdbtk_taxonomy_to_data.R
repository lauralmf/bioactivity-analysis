library(tidyverse)

# Read the GTDB-Tk results:
gtdb <- read_tsv("../gtdb_results/gtdbtk.bac120.summary.tsv")

# Ensure that your GTDB-Tk query names match the ones in your bioactivity data: 
# Here, I removed "_assembly", because I didn't have the "_assembly" after my
# sequencing ID (SEQID) in my bioactivity data:
gtdb$user_genome <- gsub("_assembly$", "", gtdb$user_genome) 

for (i in 1:nrow(data)){
  for (j in 1:nrow(gtdb)){
    if (data$SEQID[i] == gtdb$user_genome[j]){
      data$Classification[i] <- gtdb$classification[j]
      break
    }
  }
  if (is.na(data$Classification[i])) {
    data$Classification[i] <- data$Classification[i]
  }
}

# Add separate Phylum, Genus and Strain columns
data <- data %>%
  mutate(
    Phylum = ifelse(
      str_detect(Classification, "p__"),
      str_extract(Classification, "p__[^;\\s]+") %>% str_remove("p__"),
      NA),
    
    Genus = ifelse(
      str_detect(Classification, "g__"),
      str_extract(Classification, "g__[^;\\s]+") %>% 
        str_remove("g__") %>% 
        str_extract("^[^_]+"), # Extract before "_"
      ifelse(str_detect(Classification, " sp\\.$"), 
             str_extract(Classification, "^[^\\s]+"), 
             NA)),
    
    Strain = ifelse(
      str_detect(Classification, "s__"),
      str_extract(Classification, "s__[^;]+") %>% str_remove("s__"),
      paste0(Genus, " sp.")
    )
  )

# For those that do not have a specific species identifier (in my case):
strain_phylum_map <- data.frame(
  Strain = c("Arthrobacter sp.", "Streptomyces sp.", "Peribacillus sp.", 
             "Lysinibacillus sp.", "Psychrobacillus sp.", "Luteibacter sp.",
             "Kitasatospora sp.", "Curtobacterium sp.", "Variovorax sp.", 
             "Pseudomonas_E sp."),
  Phylum = c("Actinomycetota", "Actinomycetota", "Bacillota", "Bacillota", 
             "Bacillota", "Pseudomonadota", "Actinomycetota", "Actinomycetota", 
             "Pseudomonadota", "Pseudomonadota")
)

data <- data %>%
  mutate(Phylum = case_when(
    !is.na(Phylum) ~ Phylum,
    Strain %in% strain_phylum_map$Strain ~ strain_phylum_map$Phylum[match(
      Strain, strain_phylum_map$Strain)], TRUE ~ NA_character_))

# Use Pseudomonas as genus instead of Pseudomonas_E
data <- data %>%
  mutate(Genus = ifelse(Genus == "Pseudomonas_E", "Pseudomonas", Genus))

# Remove GTDB-Tk-formatted classification column
data <- data %>% select(-Classification)

# To get the strain name with the sequencing identifier for unique names:
data <- data %>% 
  mutate(Strain_full = paste0(Strain, " (", SEQID, ")"))
