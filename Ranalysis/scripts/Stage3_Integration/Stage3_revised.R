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
  # Keep only clade read abundances and coerce safely to numeric matrix for correlation
  kraken_num <- kraken_mat %>%
    select(name, ends_with("_cladeReads")) %>%
    rename_with(~ gsub("_cladeReads", "", .)) %>%
    tibble::column_to_rownames("name") %>%
    as.matrix()
  storage.mode(kraken_num) <- "numeric"
  kraken_num[is.na(kraken_num)] <- 0

  # 1. Co-occurrence (Spearman Correlation)
  sp_cor <- suppressWarnings(cor(t(kraken_num), method = "spearman", use = "pairwise.complete.obs"))
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
kraken_df <- read.csv(kraken_path, check.names = FALSE, stringsAsFactors = FALSE)
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
#!/usr/bin/env Rscript

# Stage3_Integrated_Network_Build.R
# Purpose: Full heterogeneous network construction with multi-level parallelization.
# Handles Species-Gene, Species-GO, GO-Gene, and Species-Species layers.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(igraph)
  library(fgsea)
  library(msigdbr)
  library(foreach)
  library(doParallel)
  library(data.table)
  library(tibble)
})

source(here("Ranalysis", "scripts", "Stage1_Data_acquisition", "utils.R"))

# --- GLOBAL OPTIMIZATIONS ---
msig <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
GO_PATHWAYS <- split(msig$gene_symbol, msig$gs_name)

# --- PART 1: SPECIES-GENE EDGES ---
build_species_gene_edges <- function(deseq_dir) {
  cat("Part 1: Building Species-Gene edges...\n")
  files <- list.files(deseq_dir, pattern = "_DESeq2_results\\.csv$", full.names = TRUE)
  purrr::map_dfr(files, function(f) {
    d <- fread(f)
    sp_name <- str_replace_all(str_remove(basename(f), "_DESeq2_results\\.csv$"), "_", " ")
    d %>%
      filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
      transmute(from = sp_name, to = as.character(gene), weight = log2FoldChange, 
                type = "Species_Gene", direction = ifelse(log2FoldChange > 0, "Positive", "Negative"))
  })
}

# --- PART 2: PARALLEL ENRICHMENT ---
run_species_enrichment_parallel <- function(deseq_dir, num_cores) {
  cat("Part 2: Running parallel fgsea enrichment...\n")
  files <- list.files(deseq_dir, pattern = "_DESeq2_results\\.csv$", full.names = TRUE)
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterExport(cl, "GO_PATHWAYS")
  
  results <- foreach(f = files, .combine = rbind, .packages = c("dplyr", "fgsea", "stringr", "data.table")) %dopar% {
    d <- fread(f)
    if(nrow(d) < 10) return(NULL)
    sp_name <- str_replace_all(str_remove(basename(f), "_DESeq2_results\\.csv$"), "_", " ")
    ranks <- sort(setNames(d$log2FoldChange, d$gene), decreasing = TRUE)
    fg <- fgsea(pathways = GO_PATHWAYS, stats = ranks, minSize = 10, maxSize = 500, nproc = 1)
    fg %>% filter(NES > 0, padj < 0.1) %>% mutate(species = sp_name)
  }
  stopCluster(cl)
  return(as_tibble(results))
}

# --- PART 3: PARALLELIZED SPECIES-SPECIES EDGES ---
# Parallelizes both Correlation and Risk Group logic 
build_species_species_parallel <- function(kraken_path, epathogen_db, num_cores, cor_threshold = 0.6) {
  cat("Part 3: Building Species-Species edges in parallel...\n")
  raw_kraken <- fread(kraken_path)
  kraken_mat <- raw_kraken %>%
    select(name, ends_with("_cladeReads")) %>%
    rename_with(~ gsub("_cladeReads", "", .)) %>%
    column_to_rownames("name") %>% 
    as.matrix()
  storage.mode(kraken_mat) <- "numeric"
  kraken_mat[is.na(kraken_mat)] <- 0
  
  # Filter low-abundance to speed up cor()
  keep_sp <- rowSums(kraken_mat > 0, na.rm = TRUE) > (ncol(kraken_mat) * 0.1)
  kraken_mat <- kraken_mat[keep_sp, ]
  
  # 1. Parallel Spearman Correlation
  # We split the matrix rows and calculate correlations in parallel chunks 
  sp_names <- rownames(kraken_mat)
  n_sp <- length(sp_names)

  if (n_sp < 2) {
    warning("Not enough species after filtering to compute correlations; returning only risk-group edges.")
    risk_map <- epathogen_db %>%
      select(Name, Human.classification) %>%
      distinct() %>%
      filter(Human.classification != "NotAnnotated")
    risk_groups <- unique(risk_map$Human.classification)

    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    shared_risk_edges <- foreach(rg = risk_groups, .combine = rbind, .packages = "dplyr") %dopar% {
      members <- risk_map %>% filter(Human.classification == rg) %>% pull(Name)
      if(length(members) < 2) return(NULL)

      expand.grid(from = members, to = members, stringsAsFactors = FALSE) %>%
        filter(from < to) %>%
        mutate(weight = 1.0, type = "Species_Species_RiskGroup")
    }
    stopCluster(cl)
    return(shared_risk_edges)
  }

  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  cat("Calculating correlations using chunks...\n")
  co_occ_edges <- foreach(i = 1:(n_sp-1), .combine = rbind, .packages = "dplyr") %dopar% {
    target_row <- kraken_mat[i, , drop=FALSE]
    other_rows <- kraken_mat[(i+1):n_sp, , drop=FALSE]
    
    # Calculate one-to-many Spearman correlations
    cors <- suppressWarnings(cor(t(target_row), t(other_rows), method = "spearman", use = "pairwise.complete.obs"))
    
    data.frame(
      from = sp_names[i],
      to = sp_names[(i+1):n_sp],
      weight = as.numeric(cors)
    ) %>% filter(!is.na(weight), weight > cor_threshold)
  }
  
  # 2. Parallel Risk Group expansion
  cat("Calculating shared Risk Groups...\n")
  risk_map <- epathogen_db %>% select(Name, Human.classification) %>% distinct() %>% filter(Human.classification != "NotAnnotated")
  risk_groups <- unique(risk_map$Human.classification)
  
  shared_risk_edges <- foreach(rg = risk_groups, .combine = rbind, .packages = "dplyr") %dopar% {
    members <- risk_map %>% filter(Human.classification == rg) %>% pull(Name)
    if(length(members) < 2) return(NULL)
    
    # Efficiently build pairs for this specific group
    expand.grid(from = members, to = members, stringsAsFactors = FALSE) %>%
      filter(from < to) %>%
      mutate(weight = 1.0, type = "Species_Species_RiskGroup")
  }
  
  stopCluster(cl)
  
  co_occ_edges <- co_occ_edges %>% mutate(type = "Species_Species_Correlation")
  return(bind_rows(co_occ_edges, shared_risk_edges))
}

# --- PART 4: MAIN PIPELINE ---
num_cores <- max(1, parallel::detectCores() - 2)
deseq_dir <- here("outputs", "DESeq2_results")
epathogen_path <- get_latest_timestamped_file("Ranalysis/databases", "^epathogen.*\\.csv$")
kraken_path <- get_latest_timestamped_file("data/processed", "^unaligned_merged.*\\.csv$")

# Build Layers
sp_gene_edges <- build_species_gene_edges(deseq_dir)
enrichment <- run_species_enrichment_parallel(deseq_dir, num_cores)

gene_go_edges <- enrichment %>% tidyr::unnest(leadingEdge) %>%
transmute(from = pathway, to = leadingEdge, weight = NES, type = "GO_Gene")

sp_go_edges <- enrichment %>% transmute(from = species, to = pathway, weight = NES, type = "Species_GO")

epathogen_db <- fread(epathogen_path)
sp_sp_edges <- build_species_species_parallel(kraken_path, epathogen_db, num_cores)

# Assembly
cat("Part 4: Assembling final Heterogeneous Network...\n")
all_edges <- bind_rows(sp_gene_edges, sp_go_edges, gene_go_edges, sp_sp_edges)
node_ids <- unique(c(all_edges$from, all_edges$to))
nodes_df <- data.frame(id = node_ids, name = node_ids) %>%
mutate(node_class = case_when(
    id %in% enrichment$pathway ~ "GO_Term",
    id %in% sp_gene_edges$from ~ "Species",
    TRUE ~ "Gene"
))

net <- graph_from_data_frame(d = all_edges, vertices = nodes_df, directed = FALSE)
write_graph(net, here("outputs", "Heterogenous_network.graphml"), format = "graphml")
write.csv(enrichment, here("outputs", "stage3_fgsea_results.csv"), row.names = FALSE)

cat("Workflow complete. Results saved to 'outputs/'.\n")
