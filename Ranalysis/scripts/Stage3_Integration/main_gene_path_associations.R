## Import
# Load necessary libraries
library(here) # For portable path management
library(dplyr)
library(igraph)
library(foreach)
library(doParallel)

# Phase 1 inputs (from output_processing.R)
filtered_counts_path <- here("data", "processed", "filtered_gene_counts.csv")
unaligned_merged_path <- here("data", "processed", "unaligned_merged.csv")
species_list_path <- here("data", "processed", "species_list_unaligned.csv")

if (!file.exists(filtered_counts_path)) {
  stop("Missing filtered_gene_counts.csv. Run output_processing.R first.")
}
if (!file.exists(unaligned_merged_path)) {
  stop("Missing unaligned_merged.csv. Run output_processing.R first.")
}
if (!file.exists(species_list_path)) {
  stop("Missing species_list_unaligned.csv. Run output_processing.R first.")
}

cat("Loading filtered count matrix...\n")
count_mat <- read.csv(filtered_counts_path, check.names = FALSE)
if (!"Geneid" %in% colnames(count_mat)) {
  colnames(count_mat)[1] <- "Geneid"
}
rownames(count_mat) <- count_mat$Geneid
if ("Neg_Control_W" %in% colnames(count_mat)) {
  count_mat <- count_mat %>% select(-Neg_Control_W)
}
gene_mat <- count_mat %>% select(-Geneid)
cat("Final gene matrix:", nrow(gene_mat), "genes x", ncol(gene_mat), "samples\n")

# Load species data
cat("Loading species data...\n")
species_mat_wide_format <- read.csv(unaligned_merged_path)
cat("  Species matrix dimensions:", nrow(species_mat_wide_format), "x", ncol(species_mat_wide_format), "\n")

species_summary <- read.csv(species_list_path)
cat("  Species summary:", nrow(species_summary), "species\n")

# Clean Kraken Columns
cat("Cleaning species matrix...\n")
kraken_mat <- species_mat_wide_format %>%
  select(name, ends_with("_cladeReads")) 
colnames(kraken_mat) <- gsub("_cladeReads", "", colnames(kraken_mat))
rownames(kraken_mat) <- kraken_mat$name
kraken_mat <- kraken_mat %>% select(-name, -starts_with("Neg_Control"))
cat("  Kraken matrix:", nrow(kraken_mat), "species x", ncol(kraken_mat), "samples\n")

# Align Samples
common_samples <- intersect(colnames(gene_mat), colnames(kraken_mat))
gene_mat <- gene_mat[, common_samples]
kraken_mat <- kraken_mat[, common_samples]
kraken_mat[is.na(kraken_mat)] <- 0
cat("Aligned", length(common_samples), "samples.\n")

# Phase 2: Differentiating Species (The "High/Low" Split) ----------------------
cat("\n=== PHASE 2: SPECIES SELECTION ===\n")

# Select Top Species
top_species_df <- species_summary %>%
  filter(Freq > quantile(Freq, 0.90, na.rm = TRUE), 
         cladeReads_mean > quantile(cladeReads_mean, 0.70, na.rm = TRUE)) %>%
  arrange(desc(Freq))

target_species <- top_species_df$name
cat("Selected", length(target_species), "species for analysis.\n")

species_gene_edges_path <- here("data", "processed", "species_gene_edges.csv")
if (!file.exists(species_gene_edges_path)) {
  stop("Missing species_gene_edges.csv. Run Stage2 run_deseq2_species_gene.R first.")
}
species_gene_edges <- read.csv(species_gene_edges_path)
cat("Loaded", nrow(species_gene_edges), "species-gene edges from file.\n")

# 2a. Build and save Species-Gene network before adding Gene-Gene edges
species_gene_nodes <- data.frame(id = unique(c(species_gene_edges$from, species_gene_edges$to)))
species_gene_nodes$node_class <- ifelse(species_gene_nodes$id %in% target_species, "Species", "Gene")

species_gene_net <- graph_from_data_frame(
  d = species_gene_edges[, c("from", "to", "weight", "type", "direction", "p_value")],
  vertices = species_gene_nodes,
  directed = FALSE
)

cat("Writing Species-Gene network to GraphML...\n")
write_graph(species_gene_net, file = "species_gene_network.graphml", format = "graphml")
cat("Species-Gene network saved as species_gene_network.graphml\n")

# Phase 3: Gene-Gene Edges via Clustering (Parallelized Version) ----------------------

relevant_genes <- unique(species_gene_edges$to)

if(length(relevant_genes) > 5) {
  sub_gene_mat <- gene_mat[relevant_genes, , drop=FALSE]
  
  # Calculate distance and cluster
  dist_mat <- as.dist(1 - cor(t(sub_gene_mat), method = "spearman"))
  hc <- hclust(dist_mat, method = "ward.D2")
  
  # Dynamic tree cut: Split into modules
  gene_clusters <- cutree(hc, k = min(15, length(relevant_genes)-1)) 
  unique_clusters <- unique(gene_clusters)
  
  # Setup parallel processing for clustering
  num_cores <- min(24, max(1, parallel::detectCores() - 1))
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Parallelize cluster processing
  gene_gene_edges_list <- foreach(
    k = unique_clusters,
    .packages = c('dplyr'),
    .errorhandling = 'remove'
  ) %dopar% {
    cluster_genes <- names(gene_clusters)[gene_clusters == k]
    if(length(cluster_genes) < 2) {
      return(NULL)
    }
    
    cluster_cor <- cor(t(sub_gene_mat[cluster_genes, , drop=FALSE]), method="spearman")
    
    # Remove lower triangle and diagonal to avoid self-loops and duplicates
    cluster_cor[lower.tri(cluster_cor, diag=TRUE)] <- NA 
    
    # Flatten and convert
    cluster_df <- as.data.frame(as.table(cluster_cor))
    
    # Filter and rename using explicit dplyr:: calls
    cluster_edges <- cluster_df %>%
      dplyr::filter(!is.na(Freq), abs(Freq) > 0.7) %>% 
      dplyr::rename(from = Var1, to = Var2, weight = Freq) %>%
      dplyr::mutate(
        from = as.character(from),
        to = as.character(to),
        type = "Gene_Gene", 
        direction = ifelse(weight > 0, "Positive", "Negative"),
        p_value = NA
      )
    
    return(cluster_edges)
  }
  
  stopCluster(cl)
  
  if(length(gene_gene_edges_list) > 0) {
    gene_gene_edges <- do.call(rbind, gene_gene_edges_list)
  } else {
    gene_gene_edges <- data.frame(
      from = character(),
      to = character(),
      weight = numeric(),
      type = character(),
      direction = character(),
      p_value = character()
    )
  }
} else {
  message("Insufficient genes for clustering.")
  gene_gene_edges <- data.frame(
    from = character(),
    to = character(),
    weight = numeric(),
    type = character(),
    direction = character(),
    p_value = character()
  )
}

# Phase 4: Constructing the Dense Heterogeneous Graph --------------------------
# 1. Species-Species Co-occurrence (The Dense Layer)
sp_cor <- cor(t(kraken_mat[target_species, , drop=FALSE]), method="spearman")
sp_cor[lower.tri(sp_cor, diag=TRUE)] <- NA 

sp_sp_edges <- as.data.frame(as.table(sp_cor)) %>%
  dplyr::filter(!is.na(Freq), abs(Freq) > 0.5) %>% 
  dplyr::rename(from = Var1, to = Var2, weight = Freq) %>%
  dplyr::mutate(
    from = as.character(from),
    to = as.character(to),
    type = "Species_Species", 
    direction = "Positive",
    p_value = NA
  )

# 2b. Combine all layers safely
all_edges <- rbind(
  species_gene_edges[, c("from", "to", "weight", "type", "direction", "p_value")],
  gene_gene_edges[, c("from", "to", "weight", "type", "direction", "p_value")],
  sp_sp_edges[, c("from", "to", "weight", "type", "direction", "p_value")]
)

# 3. Build Node Metadata
nodes <- data.frame(id = unique(c(all_edges$from, all_edges$to)))
nodes$node_class <- ifelse(nodes$id %in% target_species, "Species", "Gene")

# 4. Create igraph object
net <- graph_from_data_frame(d = all_edges, vertices = nodes, directed = FALSE)

print(net)

# Count and print number of gene nodes and species nodes
gene_nodes <- nodes %>% dplyr::filter(node_class == "Gene")
species_nodes <- nodes %>% dplyr::filter(node_class == "Species")

print(paste("Total Gene Nodes:", nrow(gene_nodes)))
print(paste("Total Species Nodes:", nrow(species_nodes)))
print(paste("Total Nodes:", nrow(nodes)))

# Export for Cytoscape
# This creates a GraphML file you can File > Import > Network from File in Cytoscape
cat("Writing Heterogeneous Network to GraphML...\n")
write_graph(net, file = "heterogeneous_network.graphml", format = "graphml")
cat("Graph saved as heterogeneous_network.graphml\n")
