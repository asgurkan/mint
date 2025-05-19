#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
library(patchwork)
library(cowplot)
library(grid)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

fam_file <- args[1]
q_prefix <- args[2]
max_k <- as.integer(args[3])
out_dir <- args[4]

# Read individual IDs from FAM
fam <- read.table(fam_file, header = FALSE)
individual_ids <- fam$V1

# Initialize list for all Q files
admix_list <- list()

# Loop over K values
for (K in 1:max_k) {
  q_file <- file.path(q_prefix, paste0("K", K), paste0("pruned.", K, ".Q"))
  q_data <- read.table(q_file, header = FALSE)
  q_data <- cbind(Individual = individual_ids, q_data)

  tidy_data <- q_data %>%
    pivot_longer(cols = -Individual,
                 names_to = "Cluster",
                 values_to = "Proportion") %>%
    mutate(K = K,
           Cluster = as.integer(gsub("V", "", Cluster)))

  admix_list[[as.character(K)]] <- tidy_data
}

# Combine into one dataframe
admix_data <- bind_rows(admix_list)

# Color palette (up to max_k colors)
color_palette <- colorRampPalette(brewer.pal(8, "Set3"))(max_k)

# Create output directory if not exists
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Plot loop
for (K in 1:max_k) {
  plot_data <- admix_data %>% filter(K == !!K)

  p <- ggplot(plot_data, aes(x = factor(Individual), y = Proportion, fill = factor(Cluster))) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = color_palette[1:K]) +
    labs(title = paste("ADMIXTURE Plot for K =", K), x = "Individual", y = "Proportion") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  ggsave(filename = file.path(out_dir, paste0("admix_plot_K_", K, ".png")),
         plot = p, width = 10, height = 6)
}


