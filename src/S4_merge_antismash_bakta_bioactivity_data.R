library(tidyverse)
library(stringr)

# File paths:
data <- read_csv("./datasets/clean_encoded_bioactivity.csv") # Clean and one-hot encoded bioactivity data
bakta <- read_csv("../bakta_results/all_bakta_annos.csv")    # bakta annotations
anti_json <- read_csv("../antismash_results/all_bgc_annotations_json.csv") # BGC annotations

eskape_cols <- data %>% select(matches("DSM|AU")) %>% colnames()

# ---- Calculate a bioactivity index and add this to data ---- #
data <- data %>% 
  pivot_longer(cols = all_of(eskape_cols), names_to = "ESKAPE", values_to = "Bioactivity") %>% 
  group_by(Strain_full) %>% 
  mutate(Bioactivity_index = sum(Bioactivity == 1)) %>% 
  select(Phylum, Genus, Strain, Strain_full, Medium, Bioactivity_index, everything()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "ESKAPE", values_from = "Bioactivity")

# ---- Merge antiSMASH and bakta annotations with data ---- #

# -- Create dataframes -- #
all <- anti_json %>% distinct(SEQID)
hq <- data %>% distinct(SEQID)
lq <- setdiff(all, hq)

bgc <- anti_json %>% 
  filter(!SEQID %in% lq$SEQID) %>% 
  mutate(Present = 1) %>% 
  group_by(SEQID, `Most similar known cluster`) %>% 
  distinct() %>% 
  count(Present) %>% 
  pivot_wider(names_from = `Most similar known cluster`, values_from = n, values_fill = 0) %>% 
  select(-Present)

phage <- bakta %>% # For phage annotations only
  filter(!SEQID %in% lq$SEQID & grepl("phage", Product, ignore.case = T)) %>% 
  select(-Gene) %>% 
  group_by(SEQID, Product) %>% 
  distinct(Product) %>% 
  mutate(Present = 1) %>% 
  pivot_wider(names_from = Product, values_from = Present, values_fill = 0)

all_genes <- bakta %>% 
  filter(!SEQID %in% lq$SEQID) %>% 
  select(-Gene) %>% 
  group_by(SEQID, Product) %>% 
  distinct(Product) %>% 
  mutate(Present = 1) %>% 
  pivot_wider(names_from = Product, values_from = Present, values_fill = 0)

# -- Join with data -- #
data_all_genes <- left_join(data, all_genes, by = "SEQID")

data_bgc <- left_join(data, bgc, by = "SEQID")
data_bgc <- data_bgc %>% mutate(across(where(is.numeric), ~ replace_na(., 0)))

data_phage <- left_join(data, phage, by = "SEQID")
data_phage <- data_phage %>% mutate(across(where(is.numeric), ~ replace_na(., 0)))

data_bgc_phage <- left_join(data_bgc, phage, by = "SEQID")
data_bgc_phage <- data_bgc_phage %>% mutate(across(where(is.numeric), ~ replace_na(., 0)))

df_long <- data_bgc_phage %>% 
  pivot_longer(cols = -c(eskape_cols, SEQID, Phylum, Genus, Strain, Strain_full, Medium), names_to = "Phage_or_BGC", values_to = "Present")

# -- Write to csvs -- #
write_csv(data_all_genes, "./datasets/clean_encoded_antismash_bakta.csv")
write_csv(df_long, "./datasets/clean_encoded_antismash_bakta_long.csv")
write_csv(data_bgc, "./datasets/clean_encoded_antismash.csv")
write_csv(data_phage, "./datasets/clean_encoded_prophage.csv")
write_csv(data_bgc_phage, "./datasets/clean_encoded_antismash_prophage.csv")
