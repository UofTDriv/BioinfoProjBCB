#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

args <- commandArgs(trailingOnly = TRUE)
cfg <- list()
if (length(args) > 0) {
  for (a in args) {
    if (grepl("=", a)) {
      kv <- strsplit(a, "=")[[1]]
      cfg[[kv[1]]] <- kv[2]
    }
  }
}

batch_id <- as.integer(cfg$batch_id %||% stop("Need batch_id!"))
data_dir <- cfg$data_dir %||% "data/processed/stage3_batches"
corr_threshold <- as.numeric(cfg$gene_corr_threshold %||% 0.6)

cat(sprintf("[Worker] Starting Batch %d...\n", batch_id))

assignments <- readRDS(file.path(data_dir, "batch_assignments.rds"))
norm_genes <- readRDS(file.path(data_dir, "norm_genes.rds"))

batch_genes <- assignments$gene[assignments$batch == batch_id]
if(length(batch_genes) == 0) {
  cat("[Worker] No genes for this batch. Exiting.\n")
  q(status=0)
}

# Calculate WGCNA correlation
if (requireNamespace("WGCNA", quietly = TRUE)) {
  cat("[Worker] Calculating WGCNA correlations exclusively for assigned batch...\n")
  correlations <- WGCNA::cor(t(norm_genes[batch_genes, , drop=FALSE]), t(norm_genes), method="spearman", use="pairwise.complete.obs")
} else {
  cat("[Worker] WGCNA not available, falling back to base cor...\n")
  correlations <- cor(t(norm_genes[batch_genes, , drop=FALSE]), t(norm_genes), method="spearman", use="pairwise.complete.obs")
}

correlations[is.na(correlations)] <- 0

# Extract edges efficiently
# Using structural pruning (from < to) based on ROWNAMES and COLNAMES to avoid duplicates across batches
# Find positive or negative correlation edges past threshold
row_names <- rownames(correlations)
col_names <- colnames(correlations)

idx <- which(abs(correlations) > corr_threshold, arr.ind = TRUE)

if (nrow(idx) > 0) {
  from_names <- row_names[idx[,1]]
  to_names <- col_names[idx[,2]]
  values <- correlations[idx]
  
  # Structural structural pruning: only keep edges where from_name < to_name algebraically
  valid_order <- from_names < to_names
  
  edges <- data.table(
    from = from_names[valid_order],
    to = to_names[valid_order],
    weight = abs(values[valid_order]),
    edge_class = ifelse(values[valid_order] > 0, "Positive", "Negative"),
    raw_value = as.numeric(values[valid_order])
  )
} else {
  edges <- data.table(from=character(), to=character(), weight=numeric(), edge_class=character(), raw_value=numeric())
}

out_file <- file.path(data_dir, sprintf("batch_%d_edges.csv", batch_id))
fwrite(edges, out_file)

cat("[Worker] Exported", nrow(edges), "directional edges to", out_file, "\n")
