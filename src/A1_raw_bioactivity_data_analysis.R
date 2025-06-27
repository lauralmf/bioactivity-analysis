library(tidyverse)
library(patchwork)
library(paletteer)
library(ggtext)
library(ggrepel)
library(umap)
library(xtable)
library(caret)
library(Rtsne)
library(ape)
library(nlme)
library(ggtree)
library(phylolm)
library(ggnewscale)
library(ade4)
library(ape)
library(nlme)
library(scales)
library(viridis)
library(phylolm)

# Load data and extracting ESKAPE colnames
data <- read_csv("./datasets/clean_encoded_bioactivity.csv")
eskape_cols <- data %>% select(matches("DSM|AU")) %>% colnames()

# Define genus colors
gc <- c("#d55d33","#5c6acc","#87bb37","#9957ca","#52c35e",
  "#d35aba","#4d902e","#d74380","#5dc28e","#d14249",
  "#47c3cc","#da9334","#5f8fcc","#c0b338","#bb8dd6",
  "#8f862d","#944f89","#418146","#e3869c","#328d73",
  "#a64b59","#9fb368","#9a5e2d","#686c2c","#da9e6c")

# ---- Exploratory data analysis ----
# -- Proportion of bioactive isolates within phylum
tmp <- data %>% 
  pivot_longer(cols = matches("DSM|epidermidis"), names_to = "ESKAPE", values_to = "Bioactive") %>% 
  mutate(Bioactive = as.factor(Bioactive)) %>% 
  group_by(ESKAPE) %>% 
  mutate(ESKAPE_killed_by = sum(Bioactive == 1)) %>% 
  arrange(ESKAPE_killed_by) %>% 
  group_by(SEQID) %>% 
  mutate(ESKAPEs_killed = sum(Bioactive == 1)) %>% 
  ungroup() %>% 
  group_by(Phylum) %>% 
  mutate(strains_in_phylum = n()/12) %>% 
  group_by(Phylum, ESKAPE) %>% 
  mutate(ESKAPE_killed_by_phylum = sum(Bioactive == 1),
         "Prop. of ESKAPEs killed by phylum" = ESKAPE_killed_by_phylum/ESKAPE_killed_by,
         "Prop. of bioactive isolates in phylum" = ESKAPE_killed_by_phylum/strains_in_phylum) %>% 
  ungroup()

fill_colors <- c("Actinomycetota" = "#F8766DFF", "Bacillota" = "#00BA38FF", "Pseudomonadota" = "#619CFFFF", "None" = "#D3D3D380")

tmp <- tmp %>% 
  mutate(Phylum_fill = factor(ifelse(Bioactive == 1, Phylum, "None"), levels = names(fill_colors)))

p_phylum <- tmp %>% 
  select(Phylum, ESKAPE, ESKAPE_killed_by, ESKAPE_killed_by_phylum, strains_in_phylum, `Prop. of bioactive isolates in phylum`) %>% distinct() %>% 
  mutate(Phylum = paste0(Phylum, " (", strains_in_phylum, ")")) %>% 
  ggplot(aes(y = reorder(ESKAPE, ESKAPE_killed_by_phylum), x = ESKAPE_killed_by_phylum, fill = Phylum, alpha = `Prop. of bioactive isolates in phylum`))+
  geom_col(position = "stack")+
  guides(
    fill = guide_legend(title = "Phylum", nrow = 1, byrow = TRUE),
    alpha = guide_legend(title = "Prop. of bioactive isolates in phylum", nrow = 1, byrow = TRUE))+
  labs(x = "Number of bioactive isolates out of 193",
       y = "ESKAPE")+
  scale_x_continuous(breaks = seq(0, 140, 10), expand = expansion(mult = c(0, 0)))+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.position = "bottom",
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"))

# -- Proportion of bioactive isolates within genus
strain_counts <- data %>% 
  group_by(Genus) %>% 
  count(Genus) %>% 
  mutate(Strains_in_genus = n) %>% 
  mutate(label = paste0(Genus, " (", Strains_in_genus, ")")) %>% 
  select(-n)

genus_data <- left_join(strain_counts, data, by = "Genus") %>% 
  select(c(SEQID, Phylum, Genus, Strain, Strain_full, label, Medium, everything()))

p_genus_data <- genus_data %>% 
  pivot_longer(cols = all_of(eskape_cols), names_to = "ESKAPE", values_to = "Bioactivity") %>% 
  group_by(SEQID, ESKAPE) %>% 
  mutate(n_bioactive = sum(Bioactivity == 1)) %>% 
  mutate(p_bioactive = n_bioactive/Strains_in_genus)

p_genus <- p_genus_data %>% 
  ggplot(aes(x = ESKAPE, y = p_bioactive, fill = ESKAPE))+
  geom_col()+
  facet_wrap(~label)+
  labs(y = "Prop. of bioactive isolates per genus", x = "ESKAPE")+
  scale_fill_paletteer_d("ggthemes::Classic_Green_Orange_12")+
  theme(
    legend.position = "bottom", 
    panel.background = element_rect(fill = "#F7F7F7F7"),
    axis.title = element_text(face = "bold", size = 14), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    plot.title.position = "plot", 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_blank(),
    strip.text = element_text(size = 13, face = "bold"))

# ---- BUSCO phylogenetics ---- #
# Read BUSCO tree
tree <- read.tree("../filtered_busco_phylo_results/supermatrix/sequences/busco_tree.tree")

tree$tip.label <- gsub("^BUSCO_", "", tree$tip.label)
tree$tip.label <- gsub("_assembly\\.fasta$", "", tree$tip.label)
tiplabels <- rep(NA, length(tree$tip.label))

for (i in 1:length(tiplabels)){
  for (row in 1:nrow(data)) {
    if (tree$tip.label[i] == data[row, "SEQID"]) {
      tiplabels[i] <- data[row, "Strain_full"] %>% pull()
    }
  }
}

tree$tip.label <- tiplabels

# Prepare and rename tip data for joining
tip_data <- data %>%
  mutate(
    Bioactivity_index = factor(Bioactivity_index),
    Phylum = factor(Phylum),
    Medium = factor(Medium),
    Strain = factor(Strain),
    short_label = gsub(" \\(.*\\)", "", Strain_full)
  ) %>%
  distinct(Strain_full, Strain, Bioactivity_index, Medium, short_label, Phylum, Genus)

# Plot tree and join metadata
p_tree <- ggtree(tree, layout = "circular", branch.length = "none")
p <- p_tree %<+% tip_data

viridis_colors <- viridis(n = 12, option = "viridis", direction = 1)
my_colors <- c("0" = "white", setNames(viridis_colors, as.character(1:12)))

busco_tree <- p +
  geom_tiplab(aes(label = paste0(short_label, " on ", Medium), color = Genus), size = 4, offset = 1, show.legend = F) +
  scale_color_manual(values = gc)+
  new_scale_color() +  # Start a new color scale
  geom_tippoint(aes(color = Bioactivity_index), size = 3) +
  scale_color_manual(values = my_colors) +  # First scale for tippoint
  guides(color = guide_legend(nrow = 1)) +
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = 100, r = 0, b = 0, l = 0),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 13),
    plot.margin = margin(t = 5, r = 1, b = 2, l = 1, unit = "cm"),
    legend.spacing.x = unit(1, "cm")
  )

# ---- Phylogenetic logistic regression ---- #
tip_data_extended <- tip_data %>%
  select(-c(Phylum, Genus, Medium)) %>% 
  left_join(data, by = "Strain_full")

p <- p %<+% tip_data_extended

df <- p$data
df <- df %>% drop_na()

new_eskape_cols <- gsub(" ", "_", eskape_cols)
df <- df %>% rename_with(~ new_eskape_cols, .cols = all_of(eskape_cols))
df <- df %>% column_to_rownames("label")

tree <- compute.brlen(tree, method = "Grafen", power = 1)
mat <- vcv(tree, corr=TRUE)

results <- tibble()

for (col in new_eskape_cols) {
  formula_strain <- as.formula(paste0(col, " ~ Strain"))
  formula_genus  <- as.formula(paste0(col, " ~ Genus.x"))
  formula_phylum <- as.formula(paste0(col, " ~ Phylum.x"))
  
  formula_i <- as.formula(paste0("Bioactivity_index", " ~ Strain"))

  fit_p <- phyloglm(formula_phylum, phy = tree, data = df, method = "logistic_MLE", btol = 100, log.alpha.bound = 10)
  fit_g <- phyloglm(formula_genus, phy = tree, data = df, method = "logistic_MLE", btol = 100, log.alpha.bound = 10)
  fit_s <- phyloglm(formula_strain, phy = tree, data = df, method = "logistic_MLE", btol = 100, log.alpha.bound = 10)
  
  summary_p <- summary(fit_p)
  summary_g <- summary(fit_g)
  summary_s <- summary(fit_s)

  tmp_df <- tibble(
    Taxonomy = c(rep("Phylum", nrow(summary_p$coefficients)), 
                 rep("Genus", nrow(summary_g$coefficients)), 
                 rep("Strain", nrow(summary_s$coefficients))),
    
    Name = c(gsub("Phylum\\.x", "", rownames(summary_p$coefficients)), 
              gsub("Genus\\.x", "", rownames(summary_g$coefficients)), 
              gsub("Strain", "", rownames(summary_s$coefficients))),
    Pathogen = col,
    Estimate = c(summary_p$coefficients[, 1], summary_g$coefficients[, 1], summary_s$coefficients[, 1]),
    p.value  = c(summary_p$coefficients[, 4], summary_g$coefficients[, 4], summary_s$coefficients[, 4]),
    Significant = ifelse(p.value < 0.05, 1, 0)
  )
  
  tmp_df <- tmp_df %>% 
    mutate(p.value = formatC(p.value, format = "e", digits = 2))

  results <- bind_rows(results, tmp_df)
}

signif_results <- results %>%
  mutate(Name = ifelse(grepl("\\(Intercept\\)", Name), "Actinomycetota", Name)) %>%
  filter(Significant == 1) %>%
  select(Taxonomy, Name, Pathogen, Estimate, p.value) %>% 
  mutate(Name = case_when(
    (Name == "(Intercept)" & Taxonomy == "Phylum") ~ "Actinomycetota",
    (Name == "(Intercept)" & Taxonomy == "Genus") ~ "Agromyces",
    (Name == "(Intercept)" & Taxonomy == "Strain") ~ "Agromyces sp.", 
    TRUE ~ Name
  ))

signif_results <- signif_results %>% 
  mutate(ESKAPE = Pathogen) %>% 
  select(-Pathogen) %>% 
  mutate(ESKAPE = case_when(
    grepl("mori", ESKAPE) ~ "Enterobacter mori DSM 26271",
    grepl("epidermidis", ESKAPE) ~ "Staphylococcus epidermidis AU 24",
    grepl("coli", ESKAPE) ~ "Escherichia coli DSM 498",
    grepl("putida", ESKAPE) ~ "Pseudomonas putida DSM 6125",
    grepl("aureus", ESKAPE) ~ "Staphylococcus aureus DSM 20231",
    grepl("aeruginosa", ESKAPE) ~ "Pseudomonas aeruginosa DSM 19880",
    grepl("cloacae", ESKAPE) ~ "Enterococcus cloacae DSM 30054",
    grepl("baumannii", ESKAPE) ~ "Acinetobacter baumannii DSM 300007",
    TRUE ~ ESKAPE
  )) %>% 
  select(Taxonomy, Name, ESKAPE, Estimate, p.value)

print(xtable(signif_results, caption = "Significant Phylum-level Effects from Phylogenetic GLM"), 
      include.rownames = FALSE, 
      booktabs = TRUE)

# ---- Dimensionality reduction methods ---- #
# -- Genus
m <- data %>% 
  select(-c(SEQID, Phylum, Genus, Strain, Medium, Strain_full, Bioactivity_index)) %>% 
  as.matrix()

set.seed(8)

# -------- # PCA
pca <- prcomp(m)
l <- pca$rotation %>% as.data.frame() %>% rownames_to_column(var = "rowname")
s <- pca$x %>% as.data.frame()
v <- (pca$sdev)^2
p <- v/sum(v)

pca_df <- cbind(data, s)

arrow_scale <- 2

eskape_l <- l %>% filter(rowname %in% c(
    "Enterococcus mundtii DSM 4838", "Enterococcus faecium DSM 20477",
    "Enterobacter mori DSM 26271", "Enterobacter cloacae DSM 30054",
    "Acinetobacter baylyi DSM 14961", "Acinetobacter baumannii DSM 300007",
    "Staphylococcus epidermidis AU 24", "Staphylococcus aureus DSM 20231",
    "Pseudomonas putida DSM 6125", "Pseudomonas aeruginosa DSM 19880",
    "Klebsiella oxytoca DSM 5175", "Escherichia coli DSM 498")
)

pca_g <- ggplot(pca_df, aes(x = PC1, y = PC2))+
  geom_point(aes(color = Genus), size = 4, alpha = 0.8, show.legend = T)+
  labs(x = paste0("PC1", " (", round(p[1]*100, 1), "%)"), 
       y = paste0("PC2", " (", round(p[2]*100, 1), "%)"))+
  scale_color_manual(values = gc)+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  theme(legend.position = "bottom", 
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "#F7F7F7F7"))

# -------- # UMAP
umap <- umap(m)
umap_df <- cbind(data, as.data.frame(umap$layout))

umap_g <- ggplot(umap_df, aes(x = V1, y = V2))+
  geom_point(aes(color = Genus), size = 4, alpha = 0.8, show.legend = T)+
  labs(x = "UMAP1", y = "UMAP2")+
  scale_color_manual(values = gc)+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  theme(legend.position = "bottom", 
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "#F7F7F7F7"))

dim_g <- pca_g + umap_g + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", legend.justification = "center")

# -- Phylum
eskape_l <- l %>% 
  mutate(Phylum = case_when(
    grepl("Enterococcus", rowname) ~ "Bacillota",
    grepl("Escherichia", rowname) ~ "Bacillota",
    grepl("Staphylococcus", rowname) ~ "Bacillota",
    grepl("Acinetobacter", rowname) ~ "Actinomycetota",
    grepl("Pseudomonas", rowname) ~ "Pseudomonadota",
    grepl("Klebsiella", rowname) ~ "Pseudomonadota",
    grepl("Enterobacter", rowname) ~ "Pseudomonadota"
  ))

eskape_l <- as_tibble(eskape_l)

pca_p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Phylum))+
  geom_point(size = 4, alpha = 0.8, show.legend = T)+
  labs(x = paste0("PC1", " (", round(p[1]*100, 1), "%)"), 
       y = paste0("PC2", " (", round(p[2]*100, 1), "%)"))+
  
  geom_segment(data = eskape_l, aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale), 
               arrow = arrow(length = unit(0.2, "cm")), 
               linewidth = 1, show.legend = F) +
  
  # Add labels to loadings
  geom_text(data = eskape_l, aes(x = PC1 * arrow_scale, y = PC2 * arrow_scale, label = rowname), 
            size = 5, hjust = 0.5, vjust = -0.5, show.legend = F, fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  
  theme(legend.position = "bottom", 
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      panel.background = element_rect(fill = "#F7F7F7F7"))

umap_p <- ggplot(umap_df, aes(x = V1, y = V2))+
  geom_point(aes(color = Phylum), size = 4, alpha = 0.8, show.legend = T)+
  labs(x = "UMAP1", y = "UMAP2")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  theme(legend.position = "bottom", 
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      panel.background = element_rect(fill = "#F7F7F7F7"))

dim_p <- pca_p + umap_p + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", ,
        legend.justification = "center",
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

# Combine results
pca_comb <- dim_p / dim_g +
    theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

# ---- MANVOA ---- #
df_test <- data %>% 
  group_by(SEQID) %>% 
  pivot_longer(cols = matches("DSM|AU"), names_to = "ESKAPE", values_to = "Bioactive") %>% 
  ungroup()

# Run MANOVAs
man.p <- manova(as.matrix(data[eskape_cols]) ~ Phylum, data = data)
man.g <- manova(as.matrix(data[eskape_cols]) ~ Genus, data = data)
man.s <- manova(as.matrix(data[eskape_cols]) ~ Strain, data = data)

get_pvals <- function(man, group){
  df <- tibble("Enterococcus mundtii DSM 4838" = summary(aov(man))[[1]]$"Pr(>F)"[1],
               "Enterobacter mori DSM 26271" = summary(aov(man))[[2]]$"Pr(>F)"[1],
               "Staphylococcus epidermidis AU 24" = summary(aov(man))[[3]]$"Pr(>F)"[1],
               "Escherichia coli DSM 498" = summary(aov(man))[[4]]$"Pr(>F)"[1],
               "Pseudomonas putida DSM 6125" = summary(aov(man))[[5]]$"Pr(>F)"[1],
               "Acinetobacter baylyi DSM 14961" = summary(aov(man))[[6]]$"Pr(>F)"[1],
               "Enterococcus faecium DSM 20477" = summary(aov(man))[[7]]$"Pr(>F)"[1],
               "Staphylococcus aureus DSM 20231" = summary(aov(man))[[8]]$"Pr(>F)"[1],
               "Klebsiella oxytoca DSM 5175" = summary(aov(man))[[9]]$"Pr(>F)"[1],
               "Acinetobacter baumannii DSM 300007" = summary(aov(man))[[10]]$"Pr(>F)"[1],
               "Pseudomonas aeruginosa DSM 19880" = summary(aov(man))[[11]]$"Pr(>F)"[1],
               "Enterobacter cloacae DSM 30054" = summary(aov(man))[[12]]$"Pr(>F)"[1])
  df_long <- df %>%
    pivot_longer(cols = everything(), names_to = "ESKAPE", values_to = "p.value") %>%
    mutate(
      sig_label = case_when(
        p.value < 0.0001 ~ " ***",
        p.value < 0.001  ~ " **",
        p.value < 0.01   ~ " *",
        p.value < 0.05   ~ " .",
        p.value < 0.1    ~ " ",
        TRUE             ~ ""
      ),
      !!group := paste0(format(p.value, scientific = TRUE, digits = 3), sig_label)
    ) %>%
    select(ESKAPE, !!group)
}

pvals_phylum <- get_pvals(man.p, "Phylum")
pvals_genus  <- get_pvals(man.g, "Genus")
pvals_strain <- get_pvals(man.s, "Strain")

pvals_df <- left_join(pvals_phylum, pvals_genus, by = "ESKAPE")
pvals_df <- left_join(pvals_df, pvals_strain, by = "ESKAPE")

latex_table <- xtable(pvals_df)
print(latex_table, include.rownames = FALSE)
