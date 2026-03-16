#!/usr/bin/env Rscript

# main_gene_path_associations.R
# Builds a dense heterogeneous network (Species-Gene, Gene-Gene, Species-Species)
# Updated to utilize pre-filtered data from output_processing.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(igraph)
  library(foreach)
  library(doParallel)
})

# --- ARGUMENT PARSING (Consistent with Stage3 Metagraph) ---
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
    # Inputs from output_processing.R
    filtered_gene_counts = kv$filtered_gene_counts %||% here("data", "processed", "filtered_gene_counts.csv"),
    unaligned_merged = kv$unaligned_merged %||% here("data", "processed", "unaligned_merged.csv"),
    species_list = kv$species_list %||% here("data", "processed", "species_list_unaligned.csv"),
    # Edges from Stage 2
    species_gene_edges = kv$species_gene_edges %||% here("data", "processed", "species_gene_edges.csv"),
    # Outputs
    out_graphml = kv$out_graphml %||% here("outputs", "heterogeneous_network.graphml"),
    num_cores = as.numeric(kv$num_cores %||% 12)
  )
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# --- MAIN EXECUTION ---
main <- function() {
  cfg <- parse_args()
  
  cat("=== PHASE 1: LOADING PRE-FILTERED DATA ===\n")
  # Loading files produced by output_processing.R
  # No manual rRNA/LOC/AS filtering here as it is handled upstream
  gene_mat_raw <- read.csv(cfg$filtered_gene_counts, check.names = FALSE, row.names = 1)
  if("Geneid" %in% colnames(gene_mat_raw)) gene_mat_raw <- gene_mat_raw %>% select(-Geneid)
  
  species_mat_wide <- read.csv(cfg$unaligned_merged, check.names = FALSE)
  species_summary <- read.csv(cfg$species_list)
  species_gene_edges <- read.csv(cfg$species_gene_edges)

  # Align Kraken Matrix
  kraken_mat <- species_mat_wide %>%
    select(name, ends_with("_cladeReads")) %>%
    rename_with(~ gsub("_cladeReads", "", .))
  rownames(kraken_mat) <- kraken_mat$name
  kraken_mat <- kraken_mat %>% select(-name)
  kraken_mat[is.na(kraken_mat)] <- 0

  # Sample Alignment
  common_samples <- intersect(colnames(gene_mat_raw), colnames(kraken_mat))
  gene_mat <- gene_mat_raw[, common_samples]
  kraken_mat <- kraken_mat[, common_samples]
  cat("Aligned", length(common_samples), "samples.\n")

  cat("\n=== PHASE 2: SPECIES SELECTION ===\n")
  top_species_df <- species_summary %>%
    filter(Freq > quantile(Freq, 0.90, na.rm = TRUE), 
           cladeReads_mean > quantile(cladeReads_mean, 0.70, na.rm = TRUE))
  
  target_species <- top_species_df$name
  cat("Selected", length(target_species), "top species.\n")

  cat("\n=== PHASE 3: GENE-GENE EDGES (CLUSTERING) ===\n")
  relevant_genes <- unique(species_gene_edges$to)
  gene_gene_edges <- data.frame()

  if(length(relevant_genes) > 5) {
    sub_gene_mat <- gene_mat[intersect(relevant_genes, rownames(gene_mat)), , drop=FALSE]
    dist_mat <- as.dist(1 - cor(t(sub_gene_mat), method = "spearman"))
    hc <- hclust(dist_mat, method = "ward.D2")
    gene_clusters <- cutree(hc, k = min(15, nrow(sub_gene_mat)-1)) 
    
    cl <- makeCluster(cfg$num_cores)
    registerDoParallel(cl)
    
    gene_gene_edges <- foreach(k = unique(gene_clusters), .combine = rbind, .packages = 'dplyr') %dopar% {
      cluster_genes <- names(gene_clusters)[gene_clusters == k]
      if(length(cluster_genes) < 2) return(NULL)
      
      c_cor <- cor(t(sub_gene_mat[cluster_genes, ]), method="spearman")
      c_cor[lower.tri(c_cor, diag=TRUE)] <- NA 
      
      as.data.frame(as.table(c_cor)) %>%
        filter(!is.na(Freq), abs(Freq) > 0.7) %>% 
        transmute(from = as.character(Var1), to = as.character(Var2), 
                  weight = Freq, type = "Gene_Gene", 
                  direction = ifelse(Freq > 0, "Positive", "Negative"))
    }
    stopCluster(cl)
  }

  cat("\n=== PHASE 4: SPECIES-SPECIES EDGES ===\n")
  sp_cor <- cor(t(kraken_mat[target_species, ]), method="spearman")
  sp_cor[lower.tri(sp_cor, diag=TRUE)] <- NA 
  sp_sp_edges <- as.data.frame(as.table(sp_cor)) %>%
    filter(!is.na(Freq), abs(Freq) > 0.5) %>% 
    transmute(from = as.character(Var1), to = as.character(Var2), 
              weight = Freq, type = "Species_Species", direction = "Positive")

  cat("\n=== PHASE 5: NETWORK ASSEMBLY & METADATA FIX ===\n")
  all_edges <- bind_rows(
    species_gene_edges %>% select(from, to, weight, type, direction),
    gene_gene_edges,
    sp_sp_edges
  )

  # FIX: Ensure 'name' attribute is explicitly set for Cytoscape labels
  node_ids <- unique(c(all_edges$from, all_edges$to))
  nodes <- data.frame(
    id = node_ids,
    name = node_ids, # Cytoscape uses 'name' for labels
    node_class = ifelse(node_ids %in% target_species, "Species", "Gene")
  ) %>%
    left_join(species_summary %>% select(name, Freq, cladeReads_mean), by = "name")

  net <- graph_from_data_frame(d = all_edges, vertices = nodes, directed = FALSE)
  
  cat("Writing Heterogeneous Network to:", cfg$out_graphml, "\n")
  write_graph(net, file = cfg$out_graphml, format = "graphml")
  cat("Done.\n")
}

main()