#!/usr/bin/env Rscript

# Stage3_Integrated_Network_Build_Parallel.R
# Updated to parallelize Part 2: run_species_enrichment

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(igraph)
  library(fgsea)
  library(msigdbr)
  library(foreach)     # Added for parallel loops
  library(doParallel)  # Added for parallel backend
})

# Source utility functions from Stage 1
source(here("Ranalysis", "scripts", "Stage1_Data_acquisition", "utils.R"))

# --- PART 1: BUILD SPECIES-GENE EDGES ---
# Create edges from species to DEGs using DESeq2 results
build_species_gene_edges <- function(deseq_dir, padj_cutoff = 0.05, lfc_cutoff = 1) {
  files <- list.files(deseq_dir, pattern = "_DESeq2_results\\.csv$", full.names = TRUE)
  
  purrr::map_dfr(files, function(f) {
    d <- read.csv(f)
    sp_name <- basename(f) %>% str_remove("_DESeq2_results\\.csv$") %>% str_replace_all("_", " ")
    
    d %>%
      filter(padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff) %>%
      transmute(
        from = sp_name,
        to = as.character(gene),
        weight = log2FoldChange,
        type = "Species_Gene",
        direction = ifelse(log2FoldChange > 0, "Positive", "Negative")
      )
  })
}

# --- PART 2: RUN POSITIVE FGSEA (PARALLELIZED) ---
# Perform gene enrichment and keep positive edges
run_species_enrichment_parallel <- function(deseq_dir, num_cores = 4) {
  cat("Setting up parallel cluster with", num_cores, "cores...\n")
  
  # Load GO:BP gene sets once before entering parallel loop
  msig <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
  pathways <- split(msig$gene_symbol, msig$gs_name)
  
  files <- list.files(deseq_dir, pattern = "_DESeq2_results\\.csv$", full.names = TRUE)
  
  # Initialize Cluster
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Use foreach to iterate through files in parallel
  enrichment_results <- foreach(f = files, .combine = rbind, .packages = c("dplyr", "fgsea", "stringr")) %dopar% {
    d <- read.csv(f)
    sp_name <- basename(f) %>% 
      str_remove("_DESeq2_results\\.csv$") %>% 
      str_replace_all("_", " ")
    
    # Create ranked list for fgsea
    ranks <- d$log2FoldChange
    names(ranks) <- d$gene
    ranks <- sort(ranks, decreasing = TRUE)
    
    # Run fgsea
    fg <- fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500)
    
    # Return positive enrichment only
    fg %>%
      filter(NES > 0, padj < 0.1) %>%
      mutate(species = sp_name)
  }
  
  stopCluster(cl)
  return(as_tibble(enrichment_results))
}

# --- PART 3: SPECIES-SPECIES EDGES ---
# Build edges from Risk Group and Co-occurrence
build_species_species_edges <- function(kraken_mat, epathogen_db, cor_threshold = 0.6) {
  # 1. Co-occurrence (Spearman Correlation)
  sp_cor <- cor(t(kraken_mat), method = "spearman")
  sp_cor[lower.tri(sp_cor, diag = TRUE)] <- NA
  co_occ_edges <- as.data.frame(as.table(sp_cor)) %>%
    filter(!is.na(Freq), Freq > cor_threshold) %>%
    transmute(from = as.character(Var1), to = as.character(Var2), 
              weight = Freq, type = "Species_Species_Correlation")
  
  # 2. Risk Group sharing (from epathogen)
  risk_map <- epathogen_db %>% 
    select(Name, Human.classification) %>% 
    distinct()
  
  shared_risk_edges <- expand.grid(from = risk_map$Name, to = risk_map$Name) %>%
    filter(from < to) %>%
    left_join(risk_map, by = c("from" = "Name")) %>%
    left_join(risk_map, by = c("to" = "Name"), suffix = c(".x", ".y")) %>%
    filter(Human.classification.x == Human.classification.y, 
           Human.classification.x != "NotAnnotated") %>%
    transmute(from, to, weight = 1.0, type = "Species_Species_RiskGroup")
  
  bind_rows(co_occ_edges, shared_risk_edges)
}

# --- PART 4: BUILD GRAPH ---
# Build the final heterogeneous network.graphml

num_cores <- 22 # Adjust based on your environment
deseq_dir <- here("outputs", "DESeq2_results")
epathogen_path <- get_latest_timestamped_file("Ranalysis/databases", "^epathogen.*\\.csv$")
kraken_path <- get_latest_timestamped_file("data/processed", "^unaligned_merged.*\\.csv$")

cat("Part 1: Building Species-Gene edges...\n")
sp_gene_edges <- build_species_gene_edges(deseq_dir)

cat("Part 2: Running parallelized positive fgsea...\n")
enrichment <- run_species_enrichment_parallel(deseq_dir, num_cores = num_cores)

# Link pathways to genes (leading edge) and species to pathways
gene_go_edges <- enrichment %>%
tidyr::unnest(leadingEdge) %>%
transmute(from = pathway, to = leadingEdge, weight = NES, type = "GO_Gene")

sp_go_edges <- enrichment %>%
transmute(from = species, to = pathway, weight = NES, type = "Species_GO")

cat("Part 3: Building Species-Species edges...\n")
kraken_df <- read.csv(kraken_path, row.names = 1, check.names = FALSE)
epathogen_db <- read.csv(epathogen_path)
sp_sp_edges <- build_species_species_edges(kraken_df, epathogen_db)

cat("Part 4: Assembling heterogeneous network...\n")
all_edges <- bind_rows(sp_gene_edges, sp_go_edges, gene_go_edges, sp_sp_edges)

all_nodes <- unique(c(all_edges$from, all_edges$to))
nodes_df <- data.frame(id = all_nodes, name = all_nodes) %>%
mutate(node_class = case_when(
    id %in% enrichment$pathway ~ "GO_Term",
    id %in% sp_gene_edges$from ~ "Species",
    TRUE ~ "Gene"
))

net <- graph_from_data_frame(d = all_edges, vertices = nodes_df, directed = FALSE)

write.csv(enrichment, here("outputs", "stage3_fgsea_results.csv"), row.names = FALSE)
write_graph(net, here("outputs", "Heterogenous_network.graphml"), format = "graphml")

cat("Done. Outputs saved to the 'outputs' folder.\n")
