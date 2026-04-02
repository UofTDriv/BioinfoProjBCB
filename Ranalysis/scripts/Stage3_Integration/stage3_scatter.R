#!/usr/bin/env Rscript

# Skip renv isolation and use system packages
if (file.exists("renv/activate.R")) {
  # Don't activate renv, just ensure key packages are available
  for (pkg in c("dplyr", "tidyr", "stringr", "here")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, quiet = TRUE, verbose = FALSE)
    }
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(here)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

normalize_label <- function(x) {
  x %>% tolower() %>% str_replace_all("[^a-z0-9]", "")
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  kv <- list()
  if (length(args) > 0) {
    for (a in args) {
      if (str_detect(a, "=")) {
        parts <- str_split_fixed(a, "=", 2)
        k <- str_replace(parts[, 1], "^--", "")
        v <- parts[, 2]
        kv[[k]] <- v
      }
    }
  }

  list(
    filtered_gene_counts = kv$filtered_gene_counts %||% here("data", "processed", "filtered_gene_counts.csv"),
    species_abundance = kv$species_abundance, # user can pass bracken or unaligned path
    species_list_nonhuman = kv$species_list_nonhuman, # Optional
    deseq_dir = kv$deseq_dir %||% here("outputs", "DESeq2_results"),
    padj_threshold = as.numeric(kv$padj_threshold %||% 0.05),
    log2fc_threshold = as.numeric(kv$log2fc_threshold %||% 1.0),
    top_n_genes = as.integer(kv$top_n_genes %||% 5000),
    num_batches = as.integer(kv$num_batches %||% 50),
    out_dir = kv$out_dir %||% here("data", "processed", "stage3_batches")
  )
}

cfg <- parse_args()

# Fallback for species_abundance if not provided
if (is.null(cfg$species_abundance)) {
    cfg$species_abundance <- here("data", "processed", "species_list_unaligned.csv")
}

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

load_gene_matrix <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
  if ("Geneid" %in% colnames(x)) x <- x %>% select(-Geneid)
  if ("Neg_Control_W" %in% colnames(x)) x <- x %>% select(-Neg_Control_W)
  x[] <- lapply(x, as.numeric)
  as.matrix(x)
}

load_species_matrix <- function(path, species_list_path = NULL) {
  sp <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  locate_sample_cols <- function(df) {
    bracken_cols <- grep("_bracken_reads$", colnames(df), value = TRUE)
    if (length(bracken_cols) > 0) return(bracken_cols)
    grep("_cladeReads$", colnames(df), value = TRUE)
  }
  sample_cols <- locate_sample_cols(sp)
  
  if (length(sample_cols) == 0) {
    fallbacks <- c(
      here("data", "processed", "brackenRanalysis", "outputs", "nonhuman_bracken260330op.csv"),
      here("data", "processed", "nonhuman_merged.csv"),
      here("data", "processed", "unaligned_merged.csv")
    )
    alt <- fallbacks[file.exists(fallbacks)]
    if (length(alt) > 0) {
      warning("No sample-level abundance cols found. Falling back to: ", alt[[1]])
      sp <- read.csv(alt[[1]], check.names = FALSE, stringsAsFactors = FALSE)
      sample_cols <- locate_sample_cols(sp)
    }
  }

  if (!is.null(species_list_path) && file.exists(species_list_path)) {
    sp_list <- read.csv(species_list_path, check.names = FALSE, stringsAsFactors = FALSE)
    if ("name" %in% colnames(sp_list)) {
      sp <- sp %>% filter(name %in% sp_list$name)
    }
  }

  suffix_pattern <- if (all(grepl("_bracken_reads$", sample_cols))) "_bracken_reads$" else "_cladeReads$"
  
  mat <- sp %>%
    select(name, all_of(sample_cols)) %>%
    rename_with(~ str_replace(.x, suffix_pattern, ""), all_of(sample_cols))

  rownames(mat) <- mat$name
  mat <- mat %>% select(-name)
  mat[is.na(mat)] <- 0
  mat[] <- lapply(mat, as.numeric)

  list(species_info = sp, species_matrix = as.matrix(mat))
}

load_deseq_results <- function(deseq_dir) {
  files <- list.files(deseq_dir, pattern = "_DESeq2_results\\.csv$", full.names = TRUE)
  if (length(files) == 0) stop("No DESeq2 result files found in: ", deseq_dir)

  out <- purrr::map_dfr(files, function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    if (!all(c("gene", "log2FoldChange", "padj") %in% colnames(d))) return(tibble())
    
    species_name <- basename(f) %>%
      str_remove("_DESeq2_results\\.csv$") %>%
      str_replace_all("_", " ")

    d %>%
      transmute(
        species = species_name,
        gene = as.character(gene),
        log2FoldChange = as.numeric(log2FoldChange),
        padj = as.numeric(padj)
      )
  })
  out
}

phase1_normalize <- function(mat) {
  mat <- log2(mat + 1)
  med <- apply(mat, 2, median, na.rm = TRUE)
  sweep(mat, 2, med, FUN = "-")
}

build_species_gene_edges <- function(deseq_df, padj_thr, lfc_thr) {
  sig <- deseq_df %>%
    filter(!is.na(gene), !is.na(log2FoldChange), !is.na(padj)) %>%
    filter(padj < padj_thr, abs(log2FoldChange) > lfc_thr)

  if (nrow(sig) == 0) return(tibble(from=character(), to=character(), weight=numeric(), edge_class=character(), edge_type=character(), raw_value=numeric()))

  sig %>%
    transmute(
      from = species,
      to = gene,
      weight = abs(log2FoldChange),
      edge_class = ifelse(log2FoldChange > 0, "Positive", "Negative"),
      edge_type = "Species_Gene",
      raw_value = log2FoldChange
    ) %>% distinct()
}

cat("[Scatter] Loading gene and species matrices...\n")
gene_mat_raw <- load_gene_matrix(cfg$filtered_gene_counts)
species_loaded <- load_species_matrix(cfg$species_abundance, cfg$species_list_nonhuman)
species_mat_raw <- species_loaded$species_matrix

cat("[Scatter] Loading DESeq2 results...\n")
deseq_df <- load_deseq_results(cfg$deseq_dir)

common_samples <- intersect(colnames(gene_mat_raw), colnames(species_mat_raw))
if (length(common_samples) < 3) stop("Fewer than 3 common samples between host and microbial matrices.")

gene_mat <- gene_mat_raw[, common_samples, drop = FALSE]
species_mat <- species_mat_raw[, common_samples, drop = FALSE]

cat("[Scatter] Phase 1: log2(x+1) + median-centering...\n")
gene_norm <- phase1_normalize(gene_mat)

gene_mad <- apply(gene_norm, 1, mad, na.rm = TRUE)
if (nrow(gene_norm) > cfg$top_n_genes) {
  mad_threshold <- sort(gene_mad, decreasing = TRUE)[cfg$top_n_genes]
  keep_genes <- names(gene_mad)[gene_mad >= mad_threshold]
  filtered_out <- nrow(gene_norm) - length(keep_genes)
  gene_norm <- gene_norm[keep_genes, , drop = FALSE]
  cat(sprintf("[Scatter] MAD filtering removed %d genes; MAD threshold used: %.6f\n", filtered_out, mad_threshold))
} else {
  cat(sprintf("[Scatter] MAD filtering removed 0 genes; MAD threshold used: %.6f\n", min(gene_mad, na.rm = TRUE)))
}

species_norm <- phase1_normalize(species_mat)

cat("[Scatter] Building Species_Gene edges...\n")
species_gene_edges <- build_species_gene_edges(deseq_df, cfg$padj_threshold, cfg$log2fc_threshold)

cat("[Scatter] Generating batch assignments...\n")
gene_use <- intersect(unique(species_gene_edges$to), rownames(gene_norm))
if (length(gene_use) < 100) gene_use <- rownames(gene_norm)

set.seed(42)
k_res <- kmeans(gene_norm[gene_use, ], centers = min(cfg$num_batches, length(gene_use) - 1))
batch_assignments <- data.frame(gene = names(k_res$cluster), batch = as.integer(k_res$cluster), stringsAsFactors = FALSE)

# Assign remaining genes not involved evenly to make sure we have a full gene-gene network if needed
remaining_genes <- setdiff(rownames(gene_norm), gene_use)
if(length(remaining_genes) > 0) {
  rem_assignments <- data.frame(gene = remaining_genes, batch = sample(1:cfg$num_batches, length(remaining_genes), replace=TRUE))
  batch_assignments <- bind_rows(batch_assignments, rem_assignments)
}

saveRDS(gene_norm, file.path(cfg$out_dir, "norm_genes.rds"))
saveRDS(species_norm, file.path(cfg$out_dir, "norm_species.rds"))
saveRDS(species_gene_edges, file.path(cfg$out_dir, "species_gene_edges.rds"))
saveRDS(batch_assignments, file.path(cfg$out_dir, "batch_assignments.rds"))
saveRDS(deseq_df, file.path(cfg$out_dir, "deseq_df.rds"))

cat("[Scatter] Saved objects explicitly for Gather/Workers in:", cfg$out_dir, "\n")
