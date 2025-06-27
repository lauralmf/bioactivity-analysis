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

# Read data, get ESKAPE and BGC cols
data <- read_csv("./datasets/clean_encoded_antismash.csv")
raw_data <- read_csv("./datasets/clean_encoded_bioactivity.csv")
eskape_cols <- raw_data %>% select(matches("DSM|AU")) %>% colnames()

bgc_cols <- data %>%
  select(-c(all_of(eskape_cols), Medium, SEQID, Phylum, Genus, Strain, Strain_full, Bioactivity_index)) %>%
  colnames()

# Prepare dataframe for modeling
df_long <- data %>% 
  select(all_of(eskape_cols), all_of(bgc_cols)) %>% 
  pivot_longer(cols = all_of(eskape_cols), names_to = "ESKAPE", values_to = "Bioactivity")

X <- as.matrix(df_long %>% select(-c(ESKAPE, Bioactivity)))
y <- df_long$Bioactivity

# ---- Lasso logistic regression ---- #
# -- Global modeling
cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = 1)

betas <- as.matrix(coef(cv_fit, s = "lambda.min"))
llr_all <- data.frame(BGC = rownames(betas), Coef = betas[,1])
llr_all <- llr_all %>% filter(BGC != "(Intercept)" & Coef > 0.05) %>% arrange(desc(abs(Coef))) %>% head(100)

llr_all$BGC <- gsub("γ", "gamma", llr_all$BGC)
llr_all$BGC <- gsub("α", "alpha", llr_all$BGC)
llr_all$BGC <- gsub("β", "beta", llr_all$BGC)
llr_all$BGC <- gsub("ε", "epsilon", llr_all$BGC)

# Plot results
llr_all_p <- ggplot(llr_all, aes(x = reorder(BGC, -Coef), y = Coef, fill = Coef), color = "black")+
  geom_col(show.legend = F)+
  scale_fill_continuous_diverging(palette = "Green-Brown")+
  labs(x = "BGC", y = "LLR coefficient", title = "LLR global model")+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(hjust = 1, angle = 45),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14)
    # plot.margin = unit(c(-1,0,0,1), 'cm')
  )

# -- Local modeling
llr_sep <- map_dfr(eskape_cols, function(strain) {
  bgcs <- data %>% select(all_of(bgc_cols)) %>% as.matrix()
  eskape <- as.factor(data[[strain]])
  
  llr_model <- cv.glmnet(x = bgcs, y = eskape, family = "binomial", alpha = 1)
  
  print(paste("Processing:", strain, "| Unique values:", paste(unique(eskape), collapse = ", ")))
  
  betas <- as.matrix(coef(llr_model, s = "lambda.min"))
  betas_df <- data.frame(BGC = rownames(betas), Coef = betas[,1], row.names = NULL)
  betas_df <- betas_df %>% 
    filter(BGC != "(Intercept)", Coef != 0) %>%
    arrange(desc(abs(Coef))) %>%
    mutate(ESKAPE = strain)
  
  return(betas_df)
})

selected_bgcs <- llr_sep %>%
  group_by(BGC) %>%
  mutate(n_targets = n_distinct(ESKAPE), targets = paste(unique(ESKAPE), collapse = ", ")) %>% 
  ungroup() %>% 
  filter(Coef > 0.05) %>% 
  pull(BGC)

# Plot results
llr_sep_p <- llr_sep %>%
  group_by(BGC) %>%
  filter(BGC %in% selected_bgcs) %>% 
  mutate(n_targets = n_distinct(ESKAPE), targets = paste(unique(ESKAPE), collapse = ", ")) %>% 
  mutate("LLR coefficient" = Coef, BGC = factor(BGC)) %>% 
  ungroup() %>%
  ggplot(aes(x = reorder(BGC, -`LLR coefficient`), y = reorder(ESKAPE, -`LLR coefficient`), fill = `LLR coefficient`), color = "white")+
    geom_tile(show.legend = T)+
    labs(y = "ESKAPE", x = "BGC")+
    scale_fill_continuous_diverging(palette = "Green-Brown")+
    coord_flip()+
    theme(
      panel.border = element_rect(color = "grey", fill = NA),  # border grey, no fill
      panel.background = element_rect(fill = "white"),         # panel stays white
      plot.background = element_rect(fill = "white", color = NA),  # overall background
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.title.y = element_blank(),
      legend.position = "top",
      legend.title.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.direction = "horizontal"
    )

# ---- Random Forest ---- #
# -- Global modeling
set.seed(8)
rf_model <- randomForest(x = X, y = as.factor(y), importance = TRUE, mtry = 100)

rf_all <- rf_model$importance
rf_all <- as.data.frame(rf_all) %>% rownames_to_column("BGC")

# Plot results
rf_all_p <- rf_all %>% 
  arrange(-MeanDecreaseGini) %>% 
  head(50) %>% 
  select(BGC, MeanDecreaseGini) %>% 
  ggplot(aes(x = reorder(BGC, -MeanDecreaseGini), y = MeanDecreaseGini, fill = MeanDecreaseGini))+
  geom_col(show.legend = F)+
  scale_fill_continuous_diverging(palette = "Green-Brown")+
  xlab("BGC")+
  ggtitle("RF global model")+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(hjust = 1, angle = 45),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
  )

# -- Individual modeling
rf_sep <- map_dfr(eskape_cols, function(target) {
  
  X <- data %>% mutate(across(all_of(bgc_cols), as.factor)) %>% select(all_of(bgc_cols))
  y <- data[[target]]
  
  rf <- randomForest(x = X, y = as.factor(y), importance = TRUE)
  
  imp_df <- as.data.frame(rf$importance)
  imp_df <- imp_df %>% rownames_to_column("BGC") %>% select(BGC, MeanDecreaseAccuracy, MeanDecreaseGini) %>% mutate(ESKAPE = target)
})

# Plot results
rf_sep_p <- rf_sep %>% 
  arrange(-MeanDecreaseGini) %>% 
  select(BGC, ESKAPE, MeanDecreaseGini) %>% 
  head(200) %>% 
  ggplot(aes(x = ESKAPE, y = reorder(BGC, -MeanDecreaseGini), fill = MeanDecreaseGini))+
  geom_tile(show.legend = T)+
    labs(y = "ESKAPE", x = "BGC")+
      scale_fill_continuous_diverging(
    palette = "Green-Brown",
    breaks = c(0, 1, 2),  # choose values meaningful to your data
    labels = c(0, 1, 2))+
    theme(
      panel.border = element_rect(color = "grey", fill = NA),  # border grey, no fill
      panel.background = element_rect(fill = "white"),         # panel stays white
      plot.background = element_rect(fill = "white", color = NA),  # overall background
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.title.y = element_blank(),
      legend.position = "top",
      legend.title.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.direction = "horizontal"
    )

# ---- XGBoost ---- #
# -- Global modeling
set.seed(8)
dtrain <- xgb.DMatrix(data = X, label = y)

params <- list(objective = "binary:logistic", eval_metric = "logloss")
xgb_model <- xgb.train(params = params, data = dtrain, nrounds = 1000, verbose = 1, eta = 0.05)

xgb_all <- xgb.importance(model = xgb_model)
xgb_all <- as.data.frame(xgb_all)

# Plot results
xgb_all_p <- xgb_all %>% 
  arrange(-Gain) %>% 
  head(50) %>% 
  select(Feature, Gain) %>% 
  ggplot(aes(x = reorder(Feature, -Gain), y = Gain, fill = Gain))+
  geom_col(show.legend = F)+
  scale_fill_continuous_diverging(palette = "Green-Brown")+
  ggtitle("XGBoost global model")+
  xlab("BGC")+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(hjust = 1, angle = 45),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
  )

# -- Individual modeling
set.seed(8)
xgb_sep <- map_dfr(eskape_cols, function(target) {
  
  X <- data %>% select(all_of(bgc_cols)) %>% as.matrix()
  
  # Convert factor to numeric 0/1
  y <- as.numeric(as.factor(data[[target]])) - 1
  
  # Skip if not binary
  if (anyNA(y) || length(unique(y)) != 2) return(tibble())
  
  dtrain <- xgb.DMatrix(data = X, label = y)
  params <- list(objective = "binary:logistic", eval_metric = "logloss", verbosity = 0)
  
  xgb_model <- xgb.train(params = params, data = dtrain, verbose = 0, nrounds = 100)
  
  imp <- xgb.importance(model = xgb_model) %>%
    as_tibble() %>%
    mutate(ESKAPE = target)
  
  return(imp)
})

# Plot results
xgb_sep_p <- xgb_sep %>% 
  arrange(-Gain) %>% 
  select(Feature, ESKAPE, Gain) %>% 
  head(200) %>% 
  ggplot(aes(x = ESKAPE, y = reorder(Feature, -Gain), fill = Gain))+
  geom_tile(show.legend = T)+
    labs(z = "ESKAPE", y = "BGC")+
    scale_fill_continuous_diverging(
    palette = "Green-Brown",
    breaks = c(0, 0.1, 0.2),  # choose values meaningful to your data
    labels = c(0, 0.1, 0.2))+
    theme(
      panel.border = element_rect(color = "grey", fill = NA),  # border grey, no fill
      panel.background = element_rect(fill = "white"),         # panel stays white
      plot.background = element_rect(fill = "white", color = NA),  # overall background
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.title.y = element_blank(),
      legend.position = "top",
      legend.title.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.direction = "horizontal"
    )

# ---- Consensus BGCs ---- #
# -- Global models 
top_llr_all <- llr_all %>% arrange(-Coef) %>% head(100)
top_xgb_all <- xgb_all %>% arrange(-Gain) %>% head(100)
top_rf_all  <- rf_all %>% arrange(-MeanDecreaseGini) %>% head(100)

intersection <- Reduce(intersect, list(top_llr_all$BGC, top_xgb_all$Feature, top_rf_all$BGC))
tibble(intersection)

# -- Individual models
eskape_intersection <- tibble(ESKAPE = c(), Consensus = c())

for (i in eskape_cols){
  eskape_llr <- llr_sep %>% filter(ESKAPE == i)
  eskape_xgb <- xgb_sep %>% filter(ESKAPE == i)
  eskape_rf <- rf_sep %>% filter(ESKAPE == i)
  
  eskape_top_llr <- eskape_llr %>% arrange(-Coef) %>% head(50)
  eskape_top_xgb <- eskape_xgb %>% arrange(-Gain) %>% head(50)
  eskape_top_rf <- eskape_rf %>% arrange(-MeanDecreaseGini) %>% head(50)
  
  intersection <- Reduce(intersect, list(eskape_top_llr$BGC, eskape_top_xgb$Feature, eskape_top_rf$BGC))
  intersection_df <- tibble(ESKAPE = rep(i, length(intersection)), Consensus = intersection)
  
  eskape_intersection <- rbind(eskape_intersection, intersection_df)
}

# ---- Phylogenetic logistic regression ---- #
data <- data %>% 
  pivot_longer(cols = all_of(eskape_cols), names_to = "ESKAPE", values_to = "Bioactivity") %>% 
  group_by(Strain_full) %>% 
  mutate(Bioactivity_index = sum(Bioactivity == 1)) %>% 
  select(Phylum, Genus, Strain, Strain_full, Medium, Bioactivity_index, everything()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "ESKAPE", values_from = "Bioactivity")

# Read in your tree
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
  for (bgc in bgc_cols) { 
    formula <- as.formula(paste0(col, " ~ ", "`", bgc, "`"))
    tryCatch({fit_s <- phyloglm(formula, phy = tree, data = df, method = "logistic_MLE", btol = 100)}, error = function(e) {
                message("MLE failed for ", formula, e$message)
                fit <- phyloglm(formula, phy = tree, data = df, method = "logistic_IG10", btol = 100)
              })

    summary_s <- summary(fit_s)
    coef_mat <- summary_s$coefficients
    coef_mat <- coef_mat[rownames(coef_mat) != "(Intercept)", , drop = FALSE]
    
    tmp <- tibble(
      ESKAPE = col,
      BGC = bgc,
      Estimate = coef_mat[, "Estimate"],
      p.value = formatC(coef_mat[, "p.value"], format = "e", digits = 2),
    )
    results <- rbind(results, tmp)
  }
}

results_filtered <- results %>% 
  mutate(Significant = ifelse(p.value < 0.001, 1, 0)) %>% 
  mutate(ESKAPE = gsub("_", " ", ESKAPE)) %>% 
  select(ESKAPE, BGC, Estimate, p.value, Significant) %>% 
  filter(Significant == 1) %>% 
  select(-Significant)
  
latex_table <- xtable(results_filtered)
print(latex_table, include.rownames = FALSE)

# ---- Proportion-based approach to model bioactivity: proportion of isolates with BGCs that exhibit bioactivity against ESKAPE ---- #
bgc_df <- data %>%
  select(-c(SEQID, Phylum, Genus, Strain, Strain_full, Medium, Bioactivity_index), -all_of(eskape_cols))

# Count total number of individual BGCs present
frqs_df <- tibble(
  BGC = colnames(bgc_df),
  frq = colSums(bgc_df, na.rm = TRUE)
)

frqs_df[eskape_cols] <- NA_real_

for (strain in eskape_cols) {
  for (i in seq_len(nrow(frqs_df))) {
    gene <- frqs_df$BGC[i]
    is_present <- data[[gene]] == 1
    strain_vals <- data[[strain]][is_present]
    frqs_df[[strain]][i] <- sum(strain_vals, na.rm = TRUE) / frqs_df$frq[i]
  }
}

# Do a strict filtering of BGCs to exclude BGCs for which < 1 of isolates with a specific BGC exhibit bioactivity:
frqs_df <- frqs_df %>% 
  mutate(across(all_of(eskape_cols), ~ ifelse(.x < 1, 0, 1))) %>% 
  filter(rowSums(across(-BGC, ~ . != 0)) > 0)

# BGCs with proposed activity against combinations:
eskape_short <- c("E. coli", "E. mundtii", "E. mori", "E. cloacae", "E. faecium", 
                  "A. baumannii", "A. baylyi", "S. aureus", "S. epidermidis",
                  "P. putida", "P. aeruginosa", "K. oxytoca")

df <- frqs_df %>% 
  mutate("E. coli" = `Escherichia coli DSM 498`,
         "E. mundtii" = `Enterococcus mundtii DSM 4838`,
         "E. mori" = `Enterobacter mori DSM 26271`,
         "E. cloacae" = `Enterobacter cloacae DSM 30054`,
         "E. faecium" = `Enterococcus faecium DSM 20477`,
         "A. baumannii" = `Acinetobacter baumannii DSM 300007`,
         "A. baylyi" = `Acinetobacter baylyi DSM 14961`,
         "S. aureus" = `Staphylococcus aureus DSM 20231`,
         "S. epidermidis" = `Staphylococcus epidermidis AU 24`,
         "P. putida" = `Pseudomonas putida DSM 6125`,
         "P. aeruginosa" = `Pseudomonas aeruginosa DSM 19880`) %>% 
  select(-eskape_cols, -frq)

bgc_cols <- colnames(df)[1]

# Assuming df contains BGC column + 12 ESKAPE pathogen columns
eskape_short <- setdiff(names(df), "BGC")  # Adjust if you have other metadata columns

# Result list
narrow_hits <- list()

# Loop over combinations of 2 or more pathogens
for (k in 1:length(eskape_short)) {
  combos_k <- combn(eskape_short, k, simplify = FALSE)
  
  for (combo in combos_k) {
    # Get all other ESKAPE columns
    all_other_eskape <- setdiff(eskape_short, combo)
    
    # Filter: only rows where combo == 1 and all other eskape == 0
    filtered_df <- df %>%
      filter(
        rowSums(select(., all_of(combo)) == 1) == length(combo),
        rowSums(select(., all_of(all_other_eskape)) == 0) == length(all_other_eskape)
      )
    
    if (nrow(filtered_df) > 0) {
      narrow_hits[[paste(combo, collapse = " + ")]] <- filtered_df$BGC
    }
  }
}

narrow_summary <- enframe(narrow_hits, name = "Combo", value = "BGCs") %>%
  mutate(Num_BGCs = map_int(BGCs, length))

filtered_summary <- narrow_summary %>%
  mutate(Combo_size = stringr::str_count(Combo, "\\+") + 1) %>%
  arrange(desc(Num_BGCs))

# Plot results
single_eskape_n_bgc <- filtered_summary %>% 
  filter(Combo_size == 1) %>% 
  ggplot(aes(x = reorder(Combo, -Num_BGCs), y = Num_BGCs))+
  geom_col()+
  labs(y = "Number of potential narrow-spectrum BGCs against ESKAPE combinations",
       x = "ESKAPE")+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "#F7F7F7F7"),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11))
