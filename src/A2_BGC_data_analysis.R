library(paletteer)
library(tidyverse)
library(colorspace)
library(umap)
library(patchwork)
library(ggrepel)
library(scales)
library(glmnet)
library(Rtsne)
library(randomForest)
library(colorspace)
library(vegan)
library(broom)
library(combinat)
library(xgboost)
library(Matrix)
library(Boruta)
library(vegan)
library(ggtree)
library(viridis)
library(ape)
library(ggnewscale)
library(ggpubr)
library(forcats)
library(phylolm)
library(xtable)

# Read data, get ESKAPE and BGC columns:
data <- read_csv("./datasets/clean_encoded_antismash.csv")
raw_data <- read_csv("./datasets/clean_encoded_bioactivity.csv")
eskape_cols <- raw_data %>% select(matches("DSM|AU")) %>% colnames()

bgc_cols <- data %>%
  select(-c(all_of(eskape_cols), Medium, SEQID, Phylum, Genus, Strain, Strain_full, Bioactivity_index)) %>%
  colnames()

# Pivot longer dataframe
data_long <- data %>% pivot_longer(cols = all_of(bgc_cols), names_to = "BGC", values_to = "Present")

# Define genus colors
gc <- c("#d55d33","#5c6acc","#87bb37","#9957ca","#52c35e",
  "#d35aba","#4d902e","#d74380","#5dc28e","#d14249",
  "#47c3cc","#da9334","#5f8fcc","#c0b338","#bb8dd6",
  "#8f862d","#944f89","#418146","#e3869c","#328d73",
  "#a64b59","#9fb368","#9a5e2d","#686c2c","#da9e6c")

# ---- Exploratory data analysis ---- #
# -- Bioactivity index vs. number of BGCs - phylum
bgc_eda <- data %>% 
  pivot_longer(cols = all_of(bgc_cols), names_to = "BGC", values_to = "Present") %>% 
  group_by(Strain_full) %>% 
  mutate(n_BGCs = sum(Present == 1)) %>% 
  select(Phylum, Strain, Genus, Bioactivity_index, n_BGCs) %>% 
  distinct()

plot_data <- bgc_eda %>%
  count(n_BGCs, Bioactivity_index, Genus, Phylum) %>%
  arrange(n)  # smallest counts first

# Compute Pearson correlation and p-value per Phylum
cor_data <- plot_data %>%
  group_by(Phylum) %>%
  summarise(
    cor_test = list(cor.test(n_BGCs, Bioactivity_index, method = "pearson")),
    x = min(n_BGCs, na.rm = TRUE)
  ) %>%
  mutate(
    tidied = lapply(cor_test, tidy),
    R = sapply(tidied, function(x) round(x$estimate, 2)),
    p = sapply(tidied, function(x) {
      if (x$p.value < 0.001) "p < 0.001" else paste0("p = ", signif(x$p.value, 2))
    }),
    label = paste0("\u03B2 = ", R, ", ", p)
  )

# Plot with correlation text
bgc_eda <- ggplot(plot_data, aes(x = n_BGCs, y = Bioactivity_index, color = Phylum)) +
  geom_point(size = 3, show.legend = FALSE, alpha = 0.7) +
  geom_smooth(method = "glm", method.args = list(family = "poisson"), se = TRUE, color = "black") +
  geom_text(data = cor_data, aes(x = x, y = 12, label = label),
    inherit.aes = FALSE,
    color = "black",
    hjust = 0
  ) +
  facet_wrap(~Phylum, scales = "free_x") +
  scale_y_continuous(breaks = seq(0, 12, 1)) +
  labs(y = "Bioactivity index", x = "Number of BGCs") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "#F7F7F7F7"),
    strip.text = element_text(size = 14, face = "bold")
  )

# -- Bioactivity index vs. number of BGCs - genus
  bgc_eda <- data %>% 
  pivot_longer(cols = all_of(bgc_cols), names_to = "BGC", values_to = "Present") %>% 
  group_by(Strain_full) %>% 
  mutate(n_BGCs = sum(Present == 1)) %>% 
  select(Phylum, Strain, Genus, Bioactivity_index, n_BGCs) %>% 
  distinct()

plot_data <- bgc_eda %>%
  count(n_BGCs, Bioactivity_index, Genus, Phylum) %>%
  arrange(n)  # smallest counts first

# Compute Pearson correlation and p-value per Phylum
cor_data <- plot_data %>%
  group_by(Phylum) %>%
  summarise(
    cor_test = list(cor.test(n_BGCs, Bioactivity_index, method = "pearson")),
    x = min(n_BGCs, na.rm = TRUE)
  ) %>%
  mutate(
    tidied = lapply(cor_test, tidy),
    R = sapply(tidied, function(x) round(x$estimate, 2)),
    p = sapply(tidied, function(x) {
      if (x$p.value < 0.001) "p < 0.001" else paste0("p = ", signif(x$p.value, 2))
    }),
    label = paste0("\u03B2 = ", R, ", ", p)
  )

# Plot with correlation text
bgc_eda <- ggplot(plot_data, aes(x = n_BGCs, y = Bioactivity_index, color = Genus)) +
  geom_point(size = 3, show.legend = FALSE, alpha = 0.7) +
  geom_smooth(method = "glm", method.args = list(family = "poisson"), se = TRUE, color = "black") +
  geom_text(data = cor_data, aes(x = x, y = 12, label = label),
    inherit.aes = FALSE,
    color = "black",
    hjust = 0
  ) +
  facet_wrap(~Genus, scales = "free_x") +
  scale_y_continuous(breaks = seq(0, 12, 1)) +
  labs(y = "Bioactivity index", x = "Number of BGCs") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "#F7F7F7F7"),
    strip.text = element_text(size = 14, face = "bold")
  )

# ---- Dimensionality reduction methods ---- #
# -- PCA and UMAP
m <- data %>% 
  select(-c(SEQID, Phylum, Genus, Strain, Strain_full, Medium, Bioactivity_index)) %>% 
  select(-eskape_cols) %>% 
  as.matrix()

set.seed(8)

# -------- #
pca <- prcomp(m)
l <- pca$rotation %>% as.data.frame() %>% rownames_to_column(var = "rowname")
s <- pca$x %>% as.data.frame()
v <- (pca$sdev)^2
p <- v/sum(v)

pca_df <- cbind(data, s)

pca_g <- ggplot(pca_df, aes(x = PC1, y = PC2))+
  geom_point(aes(color = Genus), size = 3, alpha = 0.8, show.legend = T)+
  labs(x = paste0("PC1", " (", round(p[1]*100, 1), "%)"), 
       y = paste0("PC2", " (", round(p[2]*100, 1), "%)"))+
  scale_color_manual(values = gc)+
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

# -------- #
umap <- umap(m)
umap_df <- cbind(data, as.data.frame(umap$layout))
umap_df <- umap_df %>% 
  mutate(UMAP1 = V1, UMAP2 = V2) %>% 
  select(-c(V1, V2))

umap_g <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2))+
  geom_point(aes(color = Genus), size = 3, alpha = 0.8, show.legend = F)+
  scale_color_manual(values = gc)+
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

dim_g <- pca_g + umap_g + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "none", legend.justification = "center", plot.title = element_text(size = 20, face = "bold"))

# -- t-SNE
set.seed(8)
df_for_tsne <- data %>%
  select(-Bioactivity_index) %>% 
  select(where(is.numeric), Genus)

df_unique <- df_for_tsne %>%
  distinct(across(where(is.numeric)), .keep_all = TRUE)

m_dup <- df_unique %>%
  select(where(is.numeric)) %>%
  as.matrix()

genus_labels <- df_unique$Genus

tsne_result <- Rtsne(m_dup)

tsne_df <- data.frame(
  tSNE1 = tsne_result$Y[,1],
  tSNE2 = tsne_result$Y[,2],
  Genus = genus_labels,
  "Enterococcus mundtii DSM 4838" = df_unique$`Enterococcus mundtii DSM 4838`,
  "Enterobacter mori DSM 26271" = df_unique$`Enterobacter mori DSM 26271`,
  "Escherichia coli DSM 498" = df_unique$`Escherichia coli DSM 498`,
  "Klebsiella oxytoca DSM 5175" = df_unique$`Klebsiella oxytoca DSM 5175`,
  "Staphylococcus epidermidis AU 24" = df_unique$`Staphylococcus epidermidis AU 24`,
  "Staphylococcus aureus DSM 20231" = df_unique$`Staphylococcus aureus DSM 20231`,
  "Acinetobacter baumannii DSM 300007" = df_unique$`Acinetobacter baumannii DSM 300007`,
  "Acinetobacter baylyi DSM 14961" = df_unique$`Acinetobacter baylyi DSM 14961`,
  "Pseudomonas putida DSM 6125" = df_unique$`Pseudomonas putida DSM 6125`,
  "Pseudomonas aeruginosa DSM 19880" = df_unique$`Pseudomonas aeruginosa DSM 19880`,
  "Enterococcus faecium DSM 20477" = df_unique$`Enterococcus faecium DSM 20477`,
  "Enterobacter cloacae DSM 30054" = df_unique$`Enterobacter cloacae DSM 30054`
)

tsne_g <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2)) +
  geom_point(aes(color = Genus), size = 3, alpha = 0.8, show.legend = T)+
  scale_color_manual(values = gc)+
  theme(legend.position = "bottom", 
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      panel.background = element_rect(fill = "#F7F7F7F7"))

# Merge with PCA and UMAP plots:
dim_g_coupled <- dim_g + theme(legend.position = "none")
all_dim <- dim_g_coupled / tsne_g

# ---- PCA logistic regression ---- #
X <- pca_df %>% select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, all_of(eskape_cols))
new_eskape_cols <- gsub(" ", "_", eskape_cols)

for (k in seq_along(new_eskape_cols)) {
  X <- X %>%
    mutate(!!sym(new_eskape_cols[k]) := .data[[eskape_cols[k]]])
}

X <- X %>% select(-all_of(eskape_cols))
pc_cols <- grep("^PC", names(X), value = TRUE)

models <- list()

results <- tibble()

for (col in new_eskape_cols) {
  formula_str <- paste0("`", col, "` ~ ", paste(pc_cols, collapse = " + "))
  models[[col]] <- glm(as.formula(formula_str), data = X, family = "binomial")

  tmp <- tibble(
    ESKAPE = col,
    PC = pc_cols,
    Coefficients = models[[col]]$coefficients[pc_cols],
    Accuracy = sum(round(models[[col]]$fitted.values) == X[[col]]) / length(X[[col]])
)
  results <- rbind(results, tmp)
}

results <- results %>% drop_na()

results %>% 
  select(ESKAPE, Accuracy) %>% distinct() %>% 
  mutate(ESKAPE = gsub("_", " ", ESKAPE)) %>% 
  arrange(ESKAPE)

# ---- Facet t-SNE with color against bioactivity/no bioactivity against respective ESKAPE
bc <- c("#4A81D4", "#F2BF88")
a <- 0.8
titlesize <- 14
dotsize <- 2.5

abaum <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2))+
  geom_point(aes(color = ifelse(factor(Acinetobacter.baumannii.DSM.300007) == 0, "No", "Yes")), 
             size = dotsize, alpha = a, show.legend = T)+
  labs(title = "Acinetobacter baumannii")+
  guides(color = guide_legend(title = "Bioactivity"))+
  scale_color_manual(values = bc)+
  theme(legend.position = "bottom", 
        panel.grid = element_blank(), 
        axis.title.x = element_blank(),
        plot.title = element_text(size = titlesize),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

abaum

abayl <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2))+
  geom_point(aes(color = ifelse(factor(Acinetobacter.baylyi.DSM.14961) == 0, "No", "Yes")),
                 size = dotsize, alpha = a, show.legend = T)+
  labs(title = "Acinetobacter baylyi")+
  guides(color = guide_legend(title = "Bioactivity"))+
  scale_color_manual(values = bc)+
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = titlesize),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

efaec <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2))+
  geom_point(aes(color = ifelse(factor(Enterococcus.faecium.DSM.20477) == 0, "No", "Yes")),
         size = dotsize, alpha = a, show.legend = T)+
  labs(title = "Enterococcus faecium")+
  guides(color = guide_legend(title = "Bioactivity"))+
  scale_color_manual(values = bc)+
  theme(legend.position = "bottom", 
        panel.grid = element_blank(), 
        axis.title.x = element_blank(),
        plot.title = element_text(size = titlesize),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

emundt <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2))+
  geom_point(aes(color = ifelse(factor(Enterococcus.mundtii.DSM.4838) == 0, "No", "Yes")), size = dotsize, alpha = a, show.legend = T)+
  labs(title = "Enterococcus mundtii")+
  guides(color = guide_legend(title = "Bioactivity"))+
  scale_color_manual(values = bc)+
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = titlesize),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

ecloa <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2))+
  geom_point(aes(color = ifelse(factor(Enterobacter.cloacae.DSM.30054) == 0, "No", "Yes")), size = dotsize, alpha = a, show.legend = T)+
  labs(title = "Enterobacter cloacae")+
  guides(color = guide_legend(title = "Bioactivity"))+
  scale_color_manual(values = bc)+
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = titlesize),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

emori <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2))+
  geom_point(aes(color = ifelse(factor(Enterobacter.mori.DSM.26271) == 0, "No", "Yes")), size = dotsize, alpha = a, show.legend = T)+
  labs(title = "Enterobacter mori")+
  guides(color = guide_legend(title = "Bioactivity"))+
  scale_color_manual(values = bc)+
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = titlesize),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

paeru <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2))+
  geom_point(aes(color = ifelse(factor(Pseudomonas.aeruginosa.DSM.19880) == 0, "No", "Yes")),
           size = dotsize, alpha = a, show.legend = T)+
  labs(title = "Pseudomonas aeruginosa")+
  guides(color = guide_legend(title = "Bioactivity"))+
  scale_color_manual(values = bc)+
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = titlesize),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

pputi <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2))+
  geom_point(aes(color = ifelse(factor(Pseudomonas.putida.DSM.6125) == 0, "No", "Yes")), 
             size = dotsize, alpha = a, show.legend = T)+
  labs(title = "Pseudomonas putida")+
  guides(color = guide_legend(title = "Bioactivity"))+
  scale_color_manual(values = bc)+
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = titlesize),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

saure <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2))+
  geom_point(aes(color = ifelse(factor(Staphylococcus.aureus.DSM.20231) == 0, "No", "Yes")),
             size = dotsize, alpha = a, show.legend = T)+
  labs(title = "Staphylococcus aureus")+
  guides(color = guide_legend(title = "Bioactivity"))+
  scale_color_manual(values = bc)+
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = titlesize),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

sepid <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2))+
  geom_point(aes(color = ifelse(factor(Staphylococcus.epidermidis.AU.24) == 0, "No", "Yes")),
             size = dotsize, alpha = a, show.legend = T)+
  labs(title = "Staphylococcus epidermidis")+
  guides(color = guide_legend(title = "Bioactivity"))+
  scale_color_manual(values = bc)+
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = titlesize),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

koxyt <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2))+
  geom_point(aes(color = ifelse(factor(Klebsiella.oxytoca.DSM.5175) == 0, "No", "Yes")), 
             size = dotsize, alpha = a, show.legend = T)+
  labs(title = "Klebsiella oxytoca")+
  guides(color = guide_legend(title = "Bioactivity"))+
  scale_color_manual(values = bc)+
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        plot.title = element_text(size = titlesize),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

ecoli <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2))+
  geom_point(aes(color = ifelse(factor(Escherichia.coli.DSM.498) == 0, "No", "Yes")), 
             size = dotsize, alpha = a, show.legend = T)+
  labs(title = "Escherichia coli")+
  guides(color = guide_legend(title = "Bioactivity"))+
  scale_color_manual(values = bc)+
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = titlesize),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "#F7F7F7F7"))

aa <- abaum+abayl
et <- efaec+emundt
eb <- ecloa+emori
pp <- paeru+pputi
ss <- saure+sepid
ke <- koxyt+ecoli

tsne_bgc_coloring <- aa / et / eb / pp / ss / ke + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", legend.justification = "center")

# ---- Pairwise combinations of isolates with the same strain name ---- #
# Identify columns
exclude_cols <- c("SEQID", "Phylum", "Genus", "Strain", "Strain_full", "Medium", "Bioactivity_index")
exclude_cols <- c(exclude_cols, eskape_cols)

# Only keep columns that are NOT in the exclude list
gene_cols <- setdiff(names(data), exclude_cols)

# Function to compute gene and bioactivity diffs between two rows
compare_strain_pair <- function(df_pair, gene_cols, eskape_cols) {
  s1 <- df_pair[1, , drop = FALSE]
  s2 <- df_pair[2, , drop = FALSE]
  
  # Gene diffs (ignore NAs in comparisons)
  gene_diff <- gene_cols[which((s1[gene_cols] != s2[gene_cols]) & !is.na(s1[gene_cols]) & !is.na(s2[gene_cols]))]
  
  # Bioactivity diffs per ESKAPE (ignore NAs)
  eskape_diff <- sapply(eskape_cols, function(col) {
    if (!is.na(s1[[col]]) && !is.na(s2[[col]])) {
      return(s1[[col]] != s2[[col]])
    } else {
      return(NA)
    }
  })
  
  # Medium diff (check if either of the mediums is NA)
  medium_diff <- ifelse(!is.na(s1[["Medium"]]) & !is.na(s2[["Medium"]]), s1[["Medium"]] != s2[["Medium"]], NA)
  
  # Return comparison row
  tibble(
    Strain = s1[["Strain"]],
    Strain_full_1 = s1[["Strain_full"]],
    Strain_full_2 = s2[["Strain_full"]],
    Medium_1 = s1[["Medium"]],
    Medium_2 = s2[["Medium"]],
    Medium_Diff = medium_diff,
    GeneDiffs = list(gene_diff),
    N_GeneDiffs = length(gene_diff),
    BioactivityDiffs = list(eskape_diff),
    N_BioactivityDiffs = sum(!is.na(eskape_diff) & eskape_diff)
  )
}

# Main logic — group by Strain, pairwise compare Strain_fulls
pairwise_results <- data %>%
  group_by(Strain) %>%
  group_map(~{
    if (nrow(.x) < 2) return(NULL)  # Skip if only one genome

    pairs <- combn(.x$Strain_full, 2, simplify = FALSE)
    
    # For each pair, extract the relevant rows and compare
    map_dfr(pairs, function(pair_ids) {
      df_pair <- .x %>% filter(Strain_full %in% pair_ids)
      compare_strain_pair(df_pair, gene_cols, eskape_cols)
    })
  }) %>%
  bind_rows()

pairwise_results2 <- pairwise_results %>% 
  mutate(Genus = gsub( " .*$", "", Strain_full_1),
         Strain = gsub(" \\(.*\\)", "", Strain_full_1)) %>% 
  unnest(BioactivityDiffs) %>% 
  mutate(ESKAPE = rep(eskape_cols, 343),
         BioactivityDiffs = ifelse(BioactivityDiffs == TRUE, 1, 0))

distinct1 <- pairwise_results %>% 
  count(Strain_full_1) %>% 
  mutate(Strain_full = Strain_full_1) %>% 
  select(-Strain_full_1)

distinct2 <- pairwise_results %>% 
  count(Strain_full_2) %>% 
  mutate(Strain_full = Strain_full_2) %>% 
  select(-Strain_full_2)

distinct <- bind_rows(distinct1, distinct2)

strain_counts <- distinct %>% distinct(Strain_full) %>% 
  mutate(Strain = gsub(" \\(.*\\)", "", Strain_full)) %>% 
  select(-Strain_full) %>% 
  count(Strain) %>% 
  mutate(Label = paste0(Strain, " (", n, ")")) %>% 
  select(-n)

pairwise_results2 <- left_join(pairwise_results2, strain_counts)

rm(distinct1, distinct2, distinct)

pairwise_results2 <- pairwise_results2 %>%
  mutate(Medium_Diff = case_when(
    Medium_Diff == TRUE ~ "Different",
    Medium_Diff == FALSE ~ "Same",
    TRUE ~ NA_character_
  ))

summary_df <- pairwise_results2 %>%
  group_by(Strain, ESKAPE) %>%
  summarise(
    Prop_BioactivityDiff = mean(BioactivityDiffs, na.rm = TRUE),
    Avg_GeneDiffs = mean(N_GeneDiffs, na.rm = TRUE),
    Prop_MediumDiff = mean(Medium_Diff == "Different", na.rm = TRUE),
    Different_medium = ifelse(Prop_MediumDiff > 0, T, F)
  ) %>%
  ungroup() %>%
  mutate(Strain = factor(Strain, levels = sort(unique(Strain))))  # <-- reordering y-axis here

summary_df <- left_join(summary_df, strain_counts)

###

name1 <- "Proportion of pairs\nwith different\nculture medium"
name2 <- "Proportion of pairs\nwith different\nbioactivity"
name3 <- "Average BGC\ndifference\nbetween pairs"

alpha <- summary_df$Prop_MediumDiff
color <- summary_df$Avg_GeneDiffs
size <- summary_df$Prop_BioactivityDiff

summary_df$Label <- factor(summary_df$Label, levels = rev(sort(unique(summary_df$Label))))

# Plot
comp_plot <- ggplot(summary_df, aes(x = ESKAPE, y = Label, shape = Different_medium)) +
  geom_point(aes(size = size, color = color)) +
  scale_color_gradient(low = "#B7D3FE", high = "#F3732D", name = name3) +
  # scale_alpha(range = c(0.3, 1), name = name1) +
  scale_size(range = c(2, 6), name = name2) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17))+
  guides(shape = guide_legend(override.aes = list(size = 3), title = "Isolates with different\nculture medium"))+
  labs(x = "ESKAPE", y = "Strain") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right",
        panel.background = element_rect(fill = "white"))
