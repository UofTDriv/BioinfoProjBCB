#!/usr/bin/env Rscript

#### Header ####
# DESeq2 Species-Gene Edge Builder
# Author: Denis Rivard, modified by GitHub Copilot
# Description: Runs per-species DESeq2 and exports species-gene edges

library(here)
library(dplyr)
library(DESeq2)
library(foreach)
library(doParallel)

show_help <- function() {
  cat("Usage: Rscript run_deseq2_species_gene.R [OPTIONS]\n\n")
  cat("Inputs:\n")
  cat("  --filtered-counts <FILE> Path to filtered_gene_counts.csv\n")
  cat("  --unaligned-merged <FILE> Path to unaligned_merged.csv\n")
  cat("  --species-list <FILE>    Path to species_list_unaligned.csv\n")
  cat("Options:\n")
  cat("  --output-dir <DIR>       Output directory for per-species DESeq2 results\n")
  cat("  --edges-output <FILE>    Output CSV for species-gene edges\n")
  cat("  --cores <N>              Number of cores to use (default: detectCores - 1)\n")
  cat("  --use-existing           Skip processing if edges output exists\n")
  cat("  --help, -h               Show this help message\n")
}

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args || "-h" %in% args) {
  show_help()
  stop("Help requested. Execution stopped.")
}

# FILTERED_COUNTS_PATH <- here("data", "processed", "filtered_gene_counts.csv")
FILTERED_COUNTS_PATH <- here("data", "processed", "genes", "dataset_final_filtered_counts.csv")
UNALIGNED_MERGED_PATH <- here("data", "processed", "unaligned_merged.csv")
SPECIES_LIST_PATH <- here("data", "processed", "species_list_unaligned.csv")
OUTPUT_DIR <- here("outputs", "DESeq2_results")
EDGES_OUTPUT_PATH <- here("data", "processed", "genes", "species_gene_edges_extra_filtered.csv")
USE_EXISTING <- FALSE
NUM_CORES <- max(1, parallel::detectCores() - 1)

if (length(args) > 0) {
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--filtered-counts" && i < length(args)) {
      FILTERED_COUNTS_PATH <- normalizePath(args[i + 1], mustWork = FALSE)
      i <- i + 2
    } else if (args[i] == "--unaligned-merged" && i < length(args)) {
      UNALIGNED_MERGED_PATH <- normalizePath(args[i + 1], mustWork = FALSE)
      i <- i + 2
    } else if (args[i] == "--species-list" && i < length(args)) {
      SPECIES_LIST_PATH <- normalizePath(args[i + 1], mustWork = FALSE)
      i <- i + 2
    } else if (args[i] == "--output-dir" && i < length(args)) {
      OUTPUT_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      i <- i + 2
    } else if (args[i] == "--edges-output" && i < length(args)) {
      EDGES_OUTPUT_PATH <- normalizePath(args[i + 1], mustWork = FALSE)
      i <- i + 2
    } else if (args[i] == "--cores" && i < length(args)) {
      NUM_CORES <- max(1, as.numeric(args[i + 1]))
      i <- i + 2
    } else if (args[i] == "--use-existing") {
      USE_EXISTING <- TRUE
      i <- i + 1
    } else {
      i <- i + 1
    }
  }
}

if (USE_EXISTING && file.exists(EDGES_OUTPUT_PATH)) {
  cat("Edges output already exists. Skipping DESeq2 processing:\n")
  cat("  ", EDGES_OUTPUT_PATH, "\n")
  quit(save = "no")
}

if (!file.exists(FILTERED_COUNTS_PATH)) {
  stop("Missing filtered_gene_counts.csv. Run output_processing.R first.")
}
if (!file.exists(UNALIGNED_MERGED_PATH)) {
  stop("Missing unaligned_merged.csv. Run output_processing.R first.")
}
if (!file.exists(SPECIES_LIST_PATH)) {
  stop("Missing species_list_unaligned.csv. Run output_processing.R first.")
}

cat("Loading filtered gene counts...\n")
count_mat <- read.csv(FILTERED_COUNTS_PATH, check.names = FALSE)
if (!"Geneid" %in% colnames(count_mat)) {
  colnames(count_mat)[1] <- "Geneid"
}
rownames(count_mat) <- count_mat$Geneid

# Option 1: Using base R (clean and direct)
colnames(count_mat) <- sub("_count$", "", colnames(count_mat))

if ("Neg_Control_W" %in% colnames(count_mat)) {
  count_mat <- count_mat %>% dplyr::select(-Neg_Control_W)
}

gene_mat <- count_mat %>% dplyr::select(-Geneid)
cat("Final gene matrix:", nrow(gene_mat), "genes x", ncol(gene_mat), "samples\n")

cat("Loading species data...\n")
species_mat_wide_format <- read.csv(UNALIGNED_MERGED_PATH)
species_summary <- read.csv(SPECIES_LIST_PATH)

kraken_mat <- species_mat_wide_format %>%
  dplyr::select(name, ends_with("_cladeReads"))
colnames(kraken_mat) <- gsub("_cladeReads", "", colnames(kraken_mat))
rownames(kraken_mat) <- kraken_mat$name
kraken_mat <- kraken_mat %>% dplyr::select(-name, -starts_with("Neg_Control"))
cat("Kraken matrix:", nrow(kraken_mat), "species x", ncol(kraken_mat), "samples\n")

common_samples <- intersect(colnames(gene_mat), colnames(kraken_mat))
if (length(common_samples) == 0) {
  stop("No overlapping samples between gene and species matrices.")
}

gene_mat <- gene_mat[, common_samples]
kraken_mat <- kraken_mat[, common_samples]
kraken_mat[is.na(kraken_mat)] <- 0
cat("Aligned", length(common_samples), "samples.\n")

if (!all(c("Freq", "cladeReads_mean") %in% colnames(species_summary))) {
  stop("species_list_unaligned.csv is missing required columns: Freq, cladeReads_mean")
}

top_species_df <- species_summary %>%
  filter(Freq > quantile(Freq, 0.90, na.rm = TRUE),
         cladeReads_mean > quantile(cladeReads_mean, 0.70, na.rm = TRUE)) %>%
  arrange(desc(Freq))

target_species <- top_species_df$name
cat("Selected", length(target_species), "species for analysis.\n")

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
}


num_cores <- min(24, NUM_CORES)
cat("Setting up parallel processing with", num_cores, "cores...\n")
cl <- makeCluster(num_cores)
registerDoParallel(cl)

species_gene_edges_list <- foreach(
  i = seq_along(target_species),
  .packages = c("DESeq2", "dplyr", "here"),
  .errorhandling = "remove"
) %dopar% {
  sp <- target_species[i]
  sp_reads <- as.numeric(kraken_mat[sp, ])
  sp_reads[is.na(sp_reads)] <- 0
  n_present <- sum(sp_reads > 0)

  if (n_present < 3) {
    return(list(
      species = sp,
      status = "skipped",
      reason = "<3 samples",
      n_present = n_present,
      log_msg = paste("Species:", sp, "; Status: Skipped (<3 samples present)")
    ))
  }

  if (mean(sp_reads > 0) > 0.20) {
    cutoff <- median(sp_reads[sp_reads > 0])
    condition <- ifelse(sp_reads >= cutoff, "High", "Low")
    split_method <- "median"
  } else {
    condition <- ifelse(sp_reads > 0, "High", "Low")
    split_method <- "presence"
  }
  condition <- factor(condition, levels = c("Low", "High"))

  group_sizes <- table(condition)
  if (min(group_sizes) < 2) {
    return(list(
      species = sp,
      status = "skipped",
      reason = "<2 in group",
      n_present = n_present,
      groups = as.list(group_sizes),
      log_msg = paste("Species:", sp, "; Status: Skipped (<2 samples in one group)")
    ))
  }

  count_matrix <- round(gene_mat)
  colData <- data.frame(condition = condition, row.names = colnames(count_matrix))

  result <- tryCatch({
    dds <- DESeqDataSetFromMatrix(
      countData = count_matrix,
      colData = colData,
      design = ~ condition
    )

    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("condition", "High", "Low"))

    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)

    results_filename <- paste0(gsub(" ", "_", sp), "_DESeq2_results.csv")
    write.csv(res_df, file.path(OUTPUT_DIR, results_filename), row.names = FALSE)

    n_total <- nrow(res_df)
    n_padj_05 <- sum(!is.na(res_df$padj) & res_df$padj < 0.05)
    n_l2fc_2 <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 2)

    sig_genes_df <- res_df %>%
      dplyr::filter(
        !is.na(padj),
        padj < 0.05,
        abs(log2FoldChange) >= 2
      )

    n_final <- nrow(sig_genes_df)

    log_msg <- paste0(
      "Species: ", sp, "; ",
      "Filtering Results: ",
      "Total=", n_total, ", ",
      "padj<0.05=", n_padj_05, " (", round(100 * n_padj_05 / n_total, 1), "%), ",
      "padj<0.05&|L2FC|>=2=", n_l2fc_2, " (", round(100 * n_l2fc_2 / n_total, 1), "%), ",
      "Final=", n_final
    )

    if (n_final > 0) {
      new_edges <- sig_genes_df %>%
        dplyr::transmute(
          from = sp,
          to = gene,
          type = "Species_Gene",
          weight = abs(log2FoldChange),
          direction = ifelse(log2FoldChange > 0, "Positive", "Negative"),
          p_value = pvalue,
          padj = padj,
          logFC = log2FoldChange
        )

      list(
        species = sp,
        status = "success",
        n_present = n_present,
        split_method = split_method,
        groups = as.list(group_sizes),
        n_total = n_total,
        n_padj_05 = n_padj_05,
        n_l2fc_2 = n_l2fc_2,
        n_final = n_final,
        edges = new_edges,
        log_msg = log_msg
      )
    } else {
      list(
        species = sp,
        status = "no_edges",
        n_present = n_present,
        split_method = split_method,
        groups = as.list(group_sizes),
        n_total = n_total,
        n_padj_05 = n_padj_05,
        n_l2fc_2 = n_l2fc_2,
        n_final = 0,
        log_msg = log_msg
      )
    }
  }, error = function(e) {
    list(
      species = sp,
      status = "error",
      error_message = as.character(e$message),
      log_msg = paste("Species:", sp, "; Status: Error -", as.character(e$message))
    )
  })

  return(result)
}

stopCluster(cl)

cat("\n=== DESEQ2 RESULTS ===\n")
for (i in seq_along(species_gene_edges_list)) {
  result <- species_gene_edges_list[[i]]
  if (!is.null(result$log_msg)) {
    cat(result$log_msg, "\n")
  }
}

edge_dfs <- list()
for (i in seq_along(species_gene_edges_list)) {
  result <- species_gene_edges_list[[i]]
  if (!is.null(result$edges)) {
    edge_dfs[[length(edge_dfs) + 1]] <- result$edges
  }
}

if (length(edge_dfs) > 0) {
  species_gene_edges <- do.call(rbind, edge_dfs)
} else {
  species_gene_edges <- data.frame(
    from = character(),
    to = character(),
    type = character(),
    weight = numeric(),
    direction = character(),
    p_value = numeric(),
    padj = numeric(),
    logFC = numeric()
  )
}

cat("Generated", nrow(species_gene_edges), "species-gene edges.\n")
write.csv(species_gene_edges, EDGES_OUTPUT_PATH, row.names = FALSE)
cat("Species-gene edges saved to:", EDGES_OUTPUT_PATH, "\n")
