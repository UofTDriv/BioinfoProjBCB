#!/usr/bin/env Rscript

# Stage3_Integrated_Network_Build.R
# Based on the 4-part process in Stage3_process.pdf

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(igraph)
  library(fgsea)
  library(msigdbr)
})

# Source utility functions from Stage 1
source(here("Ranalysis", "scripts", "Stage1_Data_acquisition", "utils.R"))

# --- PART 1: BUILD SPECIES-GENE EDGES ---
# Create edges from species to DEGs using DESeq2 results [cite: 4, 5]
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

# --- PART 2: RUN POSITIVE FGSEA ---
# Perform GO:BP enrichment and keep positive edges [cite: 23, 24]
run_species_enrichment <- function(deseq_dir) {
  # Load GO:BP gene sets
  msig <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
  pathways <- split(msig$gene_symbol, msig$gs_name)
  
  files <- list.files(deseq_dir, pattern = "_DESeq2_results\\.csv$", full.names = TRUE)
  
  enrichment_results <- purrr::map_dfr(files, function(f) {
    d <- read.csv(f)
    sp_name <- basename(f) %>% str_remove("_DESeq2_results\\.csv$") %>% str_replace_all("_", " ")
    
    # Create ranked list for fgsea
    ranks <- d$log2FoldChange
    names(ranks) <- d$gene
    ranks <- sort(ranks, decreasing = TRUE)
    
    fg <- fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500)
    
    # Keep positive enrichment [cite: 24]
    fg %>%
      filter(NES > 0, padj < 0.1) %>%
      mutate(species = sp_name)
  })
  
  return(enrichment_results)
}

# --- PART 3: SPECIES-SPECIES EDGES ---
# Build edges from Risk Group and Co-occurrence [cite: 49, 51, 54]
build_species_species_edges <- function(kraken_mat, epathogen_db, cor_threshold = 0.6) {
  # 1. Co-occurrence (Spearman Correlation)
  sp_cor <- cor(t(kraken_mat), method = "spearman")
  sp_cor[lower.tri(sp_cor, diag = TRUE)] <- NA
  co_occ_edges <- as.data.frame(as.table(sp_cor)) %>%
    filter(!is.na(Freq), Freq > cor_threshold) %>%
    transmute(from = as.character(Var1), to = as.character(Var2), 
              weight = Freq, type = "Species_Species_Correlation")
  
  # 2. Risk Group sharing (from epathogen results)
  risk_map <- epathogen_db %>% select(Name, Human.classification) %>% distinct()
  
  shared_risk_edges <- expand.grid(from = risk_map$Name, to = risk_map$Name) %>%
    filter(from < to) %>%
    left_join(risk_map, by = c("from" = "Name")) %>%
    left_join(risk_map, by = c("to" = "Name"), suffix = c(".x", ".y")) %>%
    filter(Human.classification.x == Human.classification.y, Human.classification.x != "NotAnnotated") %>%
    transmute(from, to, weight = 1.0, type = "Species_Species_RiskGroup")
  
  bind_rows(co_occ_edges, shared_risk_edges)
}

# --- PART 4: BUILD GRAPH ---
# Integrate all parts into a heterogeneous network 

# Paths (adjust as needed)
deseq_dir <- here("outputs", "DESeq2_results")
epathogen_path <- get_latest_timestamped_file("Ranalysis/databases", "^epathogen.*\\.csv$")
kraken_path <- get_latest_timestamped_file("data/processed", "^unaligned_merged.*\\.csv$")

cat("Part 1: Building Species-Gene edges...\n")
sp_gene_edges <- build_species_gene_edges(deseq_dir)

cat("Part 2: Running positive fgsea...\n")
enrichment <- run_species_enrichment(deseq_dir)

# Link Species to GO Terms and GO Terms to Genes (Leading Edge)
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

# Define nodes with metadata
all_nodes <- unique(c(all_edges$from, all_edges$to))
nodes_df <- data.frame(id = all_nodes, name = all_nodes) %>%
mutate(node_type = case_when(
    id %in% enrichment$pathway ~ "GO_Term",
    id %in% sp_gene_edges$from ~ "Species",
    TRUE ~ "Gene"
))

net <- graph_from_data_frame(d = all_edges, vertices = nodes_df, directed = FALSE)

# Output enrichment results and .graphml [cite: 75, 76]
write.csv(enrichment, here("outputs", "stage3_fgsea_results.csv"), row.names = FALSE)
write_graph(net, here("outputs", "Heterogenous_network.graphml"), format = "graphml")

cat("Workflow complete. Network saved to 'Heterogenous_network.graphml'.\n")
