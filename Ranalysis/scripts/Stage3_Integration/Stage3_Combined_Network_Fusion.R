#!/usr/bin/env Rscript

# Refactored: This module now acts as a local orchestrator for the high-performance Scatter-Gather workflow.
# The monolithic implementation was deprecated in favor of partitioned execution for high compute efficiency.
# You can bypass this orchestrator and run the components directly on a cluster (e.g. SLURM) using the individual scripts.

suppressPackageStartupMessages({
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)

script_dir <- "Ranalysis/scripts/Stage3_Integration"

cat("=== Starting Stage3 Scatter-Gather Local Orchestrator (Unaligned) ===\n")

# 1. Scatter
scatter_script <- file.path(script_dir, "stage3_scatter.R")
cat("-> Launching Scatter component...\n")
system2("Rscript", args = c(scatter_script, args))

# 2. Worker (Parallelized locally)
cat("-> Launching Worker components (Local parallel Execution)...\n")
num_cores <- as.integer(max(1, parallel::detectCores(logical = TRUE) - 1))
num_batches <- 50 # Default matches scatter definition

worker_script <- file.path(script_dir, "stage3_worker.R")
mclapply(1:num_batches, function(b) {
  system2("Rscript", args = c(worker_script, paste0("batch_id=", b)))
}, mc.cores = num_cores)

# 3. Gather
gather_script <- file.path(script_dir, "stage3_gather.R")
cat("-> Launching Gather component...\n")
system2("Rscript", args = c(gather_script, args))

cat("=== Stage 3 Orchestrator Completed ===\n")
