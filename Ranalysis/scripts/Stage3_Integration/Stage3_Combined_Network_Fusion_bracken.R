#!/usr/bin/env Rscript

# Refactored: This module now acts as a local orchestrator for the high-performance Scatter-Gather workflow.
# The monolithic implementation was deprecated in favor of partitioned execution for high compute efficiency.
# Bracken variables are dynamically pushed to the scattered process.

suppressPackageStartupMessages({
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  # Add balance-mode default arguments to pass down to the scattered pipeline
  args <- c(
    "species_abundance=outputs/brackenRanalysis/outputs/nonhuman_bracken260330op.csv",
    "species_list_nonhuman=outputs/brackenRanalysis/outputs/species_list_nonhuman260330op.csv",
    "annotated_species=outputs/brackenRanalysisAllReports/outputs/annotated_species260330.csv",
    "deseq_dir=outputs/DESeq2_results_bracken",
    "gene_corr_threshold=0.75",
    "species_corr_threshold=0.6",
    "padj_threshold=0.05",
    "log2fc_threshold=1.5",
    "top_n_genes=5000",
    "out_graphml=outputs/Stage3_Network_for_RISK_bracken.graphml",
    "out_homd_graphml=outputs/Stage3_Network_HOMD_Niche_bracken.graphml",
    "out_annotation_csv=outputs/Stage3_RISK_Annotations_bracken.csv",
    "out_validation_csv=outputs/Stage3_Jaccard_Validation_bracken.csv",
    "out_fgsea_csv=outputs/Stage3_fgsea_go_bp_results_bracken.csv",
    "out_species_layer=outputs/Stage3_Species_only_bracken.graphml",
    "out_gene_layer=outputs/Stage3_Gene_only_bracken.graphml",
    "out_dir=data/processed/stage3_balance_gct75_sct60_padj05_log2fc15"
  )
}

script_dir <- "Ranalysis/scripts/Stage3_Integration"

cat("=== Starting Stage3 Scatter-Gather Local Orchestrator (Bracken) ===\n")

# 1. Scatter
scatter_script <- file.path(script_dir, "stage3_scatter.R")
cat("-> Launching Scatter component (Bracken mode)...\n")
system2("Rscript", args = c(scatter_script, args))

# 2. Worker (Parallelized locally)
cat("-> Launching Worker components (Local parallel Execution)...\n")
num_cores <- as.integer(max(1, parallel::detectCores(logical = TRUE) - 1))
num_batches <- 50 # Default matches scatter definition

# Need to parse out_dir to tell workers where to look
out_dir_arg <- grep("^out_dir=", args, value=TRUE)
data_dir <- if(length(out_dir_arg) > 0) sub("^out_dir=", "", out_dir_arg[1]) else "data/processed/stage3_balance_gct75_sct60_padj05_log2fc15"
gene_corr_arg <- grep("^gene_corr_threshold=", args, value=TRUE)
gene_corr_threshold <- if (length(gene_corr_arg) > 0) sub("^gene_corr_threshold=", "", gene_corr_arg[1]) else "0.6"

worker_script <- file.path(script_dir, "stage3_worker.R")
mclapply(1:num_batches, function(b) {
  system2("Rscript", args = c(
    worker_script,
    paste0("batch_id=", b),
    paste0("data_dir=", data_dir),
    paste0("gene_corr_threshold=", gene_corr_threshold)
  ))
}, mc.cores = num_cores)

# 3. Gather
gather_script <- file.path(script_dir, "stage3_gather.R")
cat("-> Launching Gather component...\n")
system2("Rscript", args = c(gather_script, args, paste0("data_dir=", data_dir)))

cat("=== Stage 3 Orchestrator (Bracken) Completed ===\n")
