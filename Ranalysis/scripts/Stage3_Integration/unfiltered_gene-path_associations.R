
## Import

# Filename: combine_star_counts.R

# Load necessary libraries
library(here) # For portable path management
library(dplyr)
library(purrr)
library(readr)

# 1. Define the directory containing your STAR count matrices
#    Replace this with the actual path to your 'covid2021_STAR' folder.
star_counts_dir <- here::here("data", "covid2021_STAR")

rm(count_files)

b1 <- read_tsv(here("data","covid2021_STAR","b1countMat.txt"))
b2 <- read_tsv(here("data","covid2021_STAR","b2countMat.txt"))
b3 <- read_tsv(here("data","covid2021_STAR","b3countMat.txt"))
b4 <- read_tsv(here("data","covid2021_STAR","b4countMat.txt"))
b5 <- read_tsv(here("data","covid2021_STAR","b5countMat.txt"))

# 4. Join all the individual sample data frames into a single matrix.
#    The `reduce` function iteratively joins the data frames by the 'gene_id' column.
combined_counts_matrix <- reduce(list(b1, b2, b3, b4, b5), full_join, by = "Geneid")

write.csv(combined_counts_matrix, here("data", "processed", "combined_star_gene_counts.csv"), row.names = FALSE)

count_mat <- read.csv(here("data", "processed", "combined_star_gene_counts.csv"))
nrow(count_mat) # 43388
ncol(count_mat) # 123
colnames(count_mat)
# [1] "Geneid"              "Neg_Control_W_count" "W50504632_count"     "W60804985_count"    
#   [5] "W60805434_count"     "W60805435_count"     "W61006900_count"     "W61006922_count"   ... 

# kraken2 unaligned results
species_mat_long_format <- read.csv("outputs/new_samples/unaligned_merged_long260120op.csv")
colnames(species_mat_long_format)
# [1] "taxID"               "name"                "taxRank"             "taxLineage"         
# [5] "sample"              "cladeReads"          "minimizers"          "distinct_minimizers"
species_mat_wide_format <- read.csv("outputs/new_samples/unaligned_merged260120op.csv")
colnames(species_mat_wide_format)
#  [1] "taxID"                             "name"                             
#  [3] "taxRank"                           "taxLineage"                       
#  [5] "Neg_Control_W_cladeReads"          "Neg_Control_W_minimizers"         
#  [7] "Neg_Control_W_distinct_minimizers" "W50504632_cladeReads"             
#  [9] "W50504632_minimizers"              "W50504632_distinct_minimizers"    
#  [11] "W60804985_cladeReads"              "W60804985_minimizers"             
#  [13] "W60804985_distinct_minimizers" ... 

species_summary <- read.csv("outputs/new_samples/species_list_unaligned260120op.csv")
colnames(species_summary)
# [1] "name"            "taxID"           "taxRank"         "taxLineage"      "Freq"            "cladeReads_mean"
# [7] "cladeReads_max" 
nrow(species_mat_long_format) # 194535
nrow(species_mat_wide_format) # 1965
nrow(species_summary) # 1965

# Phase 1: Data Cleaning & Alignment ------------------------------------
# First, we need to ensure the sample names match perfectly between your gene counts and Kraken2 data.

library(tidyverse)
install.packages("igraph")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
install.packages("pheatmap")
install.packages("foreach")
install.packages("doParallel")
library(igraph)
library(DESeq2)
library(pheatmap)
library(foreach)
library(doParallel)

# 1. Load Data (Assumes variables from your snippet exist)
# Clean Gene Count Columns: Remove "_count" suffix
colnames(count_mat) <- gsub("_count", "", colnames(count_mat))
rownames(count_mat) <- count_mat$Geneid
gene_mat <- count_mat %>% select(-Geneid, -Neg_Control_W) # Remove ID and Control

# Clean Kraken Columns: Remove "_cladeReads" suffix to match gene samples
# We need to ensure we are only matching the 'cladeReads' columns
kraken_mat <- species_mat_wide_format %>%
  select(name, ends_with("_cladeReads")) 
colnames(kraken_mat) <- gsub("_cladeReads", "", colnames(kraken_mat))
rownames(kraken_mat) <- kraken_mat$name
kraken_mat <- kraken_mat %>% select(-name, -starts_with("Neg_Control"))

# Align Samples: Intersect sample names
common_samples <- intersect(colnames(gene_mat), colnames(kraken_mat))
gene_mat <- gene_mat[, common_samples]
kraken_mat <- kraken_mat[, common_samples]

print(paste("Aligned", length(common_samples), "samples."))

# Phase 2: Differentiating Species (The "High/Low" Split) ------------------------------------
# 1. Select Top Species (Filter by Freq and Abundance)
# Adjust quantiles as needed to control "Top" definition
# Phase 2: Differentiating Species (The "High/Low" Split) ----------------------

# 1. Select Top Species
# We filter for high frequency (prevalence) and high average abundance
top_species_df <- species_summary %>%
  filter(Freq > quantile(Freq, 0.90, na.rm = TRUE), 
         cladeReads_mean > quantile(cladeReads_mean, 0.70, na.rm = TRUE)) %>%
  arrange(desc(Freq))

target_species <- top_species_df$name
print(paste("Selected", length(target_species), "species for analysis."))

# 2. Loop through Species to find Gene Correlations
species_gene_edges <- data.frame()

# Progress bar (optional, helpful for long loops)
pb <- txtProgressBar(min = 0, max = length(target_species), style = 3)

for(i in seq_along(target_species)) {
  sp <- target_species[i]
  
  # Fetch reads and handle NAs immediately
  sp_reads <- as.numeric(kraken_mat[sp, ])
  sp_reads[is.na(sp_reads)] <- 0 # CRITICAL FIX: Treat NAs as 0 counts
  
  # Skip if species is effectively absent (present in fewer than 3 samples)
  if(sum(sp_reads > 0) < 3) {
    setTxtProgressBar(pb, i)
    next 
  }
  
  # Define High vs Low/Absent groups
  # If species is widespread (>20% samples), split by median of non-zeros
  # Otherwise, split by Presence (High) vs Absence (Low)
  if(mean(sp_reads > 0) > 0.20) {
    cutoff <- median(sp_reads[sp_reads > 0])
    condition <- ifelse(sp_reads >= cutoff, "High", "Low")
  } else {
    condition <- ifelse(sp_reads > 0, "High", "Low")
  }
  condition <- factor(condition, levels = c("Low", "High"))
  
  # Skip if we don't have at least 2 samples in each group (statistical requirement)
  if(min(table(condition)) < 2) {
    setTxtProgressBar(pb, i)
    next
  }
  
  # Run Wilcoxon Test (Non-parametric t-test equivalent)
  # We suppress warnings to avoid the "ties" spam in console
  p_values <- suppressWarnings(
    apply(gene_mat, 1, function(g) {
      tryCatch(wilcox.test(g ~ condition)$p.value, error = function(e) NA)
    })
  )
  
  # Filter for significant genes (p < 0.05)
  # Ideally, you would use p.adjust(p_values, method="BH") < 0.1 for better stats, 
  # but we'll stick to raw p < 0.05 to ensure you get edges for the graph.
  sig_genes <- names(p_values)[!is.na(p_values) & p_values < 0.01] # Stricter cutoff (0.01) to reduce noise
  
  # Create Edges if genes found
  if(length(sig_genes) > 0){
    
    # Calculate Log2 Fold Change
    # usage of drop=FALSE ensures it stays a matrix even if only 1 gene matches
    sub_mat <- gene_mat[sig_genes, , drop=FALSE] 
    
    # Calculate means for High and Low groups
    means <- t(apply(sub_mat, 1, function(g) tapply(g, condition, mean, na.rm=TRUE)))
    
    # Check dimensions to ensure 'means' is structured correctly
    if(is.null(dim(means))) {
        # If tapply returned a list or vector due to single gene, force structure
        means <- matrix(means, nrow=1)
        colnames(means) <- c("Low", "High")
        rownames(means) <- sig_genes
    }
    
    # LogFC calculation (adding pseudocount +1 to avoid log(0))
    logFC <- log2(means[, "High"] + 1) - log2(means[, "Low"] + 1)
    
    new_edges <- data.frame(
      from = sp,
      to = sig_genes,
      type = "Species_Gene",
      weight = abs(logFC),
      direction = ifelse(logFC > 0, "Positive", "Negative"),
      p_value = p_values[sig_genes]
    )
    species_gene_edges <- rbind(species_gene_edges, new_edges)
  }
  setTxtProgressBar(pb, i)
}
close(pb)

print(paste("Generated", nrow(species_gene_edges), "species-gene edges."))

write.csv(species_gene_edges, here("data", "processed", "species_gene_edges.csv"), row.names = FALSE)
species_gene_edges <- read.csv(here("data", "processed", "species_gene_edges.csv"))

# Phase 3: Gene-Gene Edges via Clustering (Robust Version) ----------------------

relevant_genes <- unique(species_gene_edges$to)

if(length(relevant_genes) > 5) {
  sub_gene_mat <- gene_mat[relevant_genes, , drop=FALSE]
  
  # Calculate distance and cluster
  dist_mat <- as.dist(1 - cor(t(sub_gene_mat), method = "spearman"))
  hc <- hclust(dist_mat, method = "ward.D2")
  
  # Dynamic tree cut: Split into modules
  gene_clusters <- cutree(hc, k = min(15, length(relevant_genes)-1)) 
  
  gene_gene_edges <- data.frame()
  unique_clusters <- unique(gene_clusters)
  
  for(k in unique_clusters) {
    cluster_genes <- names(gene_clusters)[gene_clusters == k]
    if(length(cluster_genes) < 2) next
    
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
    
    gene_gene_edges <- rbind(gene_gene_edges, cluster_edges)
  }
} else {
  message("Insufficient genes for clustering.")
  gene_gene_edges <- data.frame()
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

# 2. Combine all layers safely
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

# Export for Cytoscape (Option A: File Export)
# This creates a GraphML file you can File > Import > Network from File in Cytoscape
write_graph(net, file = "heterogeneous_network.graphml", format = "graphml")
print("Graph saved as heterogeneous_network.graphml")

# Load graph
net <- read_graph("heterogeneous_network.graphml")

# Export for Cytoscape (Option B: Direct Connection via RCy3)
library(RCy3)
createNetworkFromIgraph(net, title = "Dense_Heterogeneous_Net")