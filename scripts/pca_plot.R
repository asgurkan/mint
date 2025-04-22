#!/usr/bin/env Rscript

# logging 
if (exists("snakemake") && !is.null(snakemake@log) && length(snakemake@log) > 0) {
  log_file <- snakemake@log[[1]]
  
  # file connection in write text mode
  log_con <- file(log_file, open = "wt")
  
  # directing both outputs and messages to the connection
  sink(log_con, type = "output")
  sink(log_con, type = "message")
}

# Load libraries
library(ggplot2)
library(ggrepel)

# input
eigenvec_path <- snakemake@input[["eigenvec_path"]]
eigenval_path <- snakemake@input[["eigenval_path"]]
sample_tsv_path <- snakemake@input[["sample_tsv"]]

# output
output_dir <- snakemake@output[["pca_output_dir"]]

# plot parameters
scatter_size <- 4
scatter_text_size <- 3
overlap_count <- 3

#  PCA eigenvalues to compute % variance explained
cat("Reading eigenvalues from:", eigenval_path, "\n")
eigs <- scan(eigenval_path)
total_var <- sum(eigs)
pct_var <- (eigs / total_var) * 100  

# PCA eigenvectors (Format: FID, IID, PC1, PC2)
cat("Reading PCA eigenvectors from:", eigenvec_path, "\n")
pca <- read.table(eigenvec_path, header = FALSE, stringsAsFactors = FALSE)
num_pcs <- ncol(pca) - 2  # subtract 2 for FID, IID
pc_cols <- paste0("PC", seq_len(num_pcs))
colnames(pca) <- c("FID", "IID", pc_cols)

cat("Reading sample TSV from:", sample_tsv_path, "\n")
metadata <- read.table(sample_tsv_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata$IID <- metadata$sample
metadata$species <- metadata$population

# columns for merging
metadata <- metadata[, c("IID", "species")]

cat("Merging PCA data with metadata...\n")
df <- merge(pca, metadata, by = "IID")

# -- PC1 vs PC2 --
pc1_var <- round(pct_var[1], 2)
pc2_var <- round(pct_var[2], 2)

cat("Generating PCA (PC1 vs PC2) plot...\n")
p1 <- ggplot(df, aes(x = PC1, y = PC2, color = species, label = IID)) +
  geom_point(size = scatter_size, alpha = 0.9) +
  geom_text_repel(size = scatter_text_size, max.overlaps = overlap_count) +
  labs(
    title = "PCA (PC1 vs PC2)",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    color = "Species"
  ) +
  theme_minimal()

cat("Saving PC1 vs PC2 plots to:", output_dir, "\n")
# save white background
ggsave(
  filename = file.path(output_dir, "PC1_PC2/pca_PC1_PC2_white.png"),
  plot = p1,
  width = 10,
  height = 8,
  dpi = 600,
  bg = "white",
  create.dir = TRUE
)

# save transparent background
ggsave(
  filename = file.path(output_dir, "PC1_PC2/pca_PC1_PC2_transparent.png"),
  plot = p1,
  width = 10,
  height = 8,
  dpi = 600,
  bg = "transparent"
)

# -- PC2 vs PC3 (only if >= 3 PCs) -- 
if (num_pcs >= 3) {
  pc3_var <- round(pct_var[3], 2)
  cat("Generating PCA (PC2 vs PC3) plot...\n")
  p2 <- ggplot(df, aes(x = PC2, y = PC3, color = species, label = IID)) +
    geom_point(size = scatter_size, alpha = 0.8) +
    geom_text_repel(size = scatter_text_size, max.overlaps = overlap_count) +
    labs(
      title = "PCA (PC2 vs PC3)",
      x = paste0("PC2 (", pc2_var, "%)"),
      y = paste0("PC3 (", pc3_var, "%)"),
      color = "Species"
    ) +
    theme_minimal()
  
  cat("Saving PC2 vs PC3 plots to:", output_dir, "\n")
  # save white background
  ggsave(
    filename = file.path(output_dir, "PC2_PC3/pca_PC2_PC3_white.png"),
    plot = p2,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white",
    create.dir = TRUE
  )
  # save transparent background
  ggsave(
    filename = file.path(output_dir, "PC2_PC3/pca_PC2_PC3_transparent.png"),
    plot = p2,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "transparent"
  )
} else {
  message("Less than 3 PCs available; skipping PC2 vs PC3 plot.")
}

cat("\nPCA plots saved successfully.\n")

if (exists("snakemake") && !is.null(snakemake@log) && length(snakemake@log) > 0) {
  sink(type = "message")  # close message sink
  sink(type = "output")   # close output sink
}