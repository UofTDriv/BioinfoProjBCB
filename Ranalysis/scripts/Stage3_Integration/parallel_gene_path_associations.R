# Filename: parallel_gene_path_associations.R
# Description: This script identifies associations between microbial species abundance 
#              and host gene expression, and builds a heterogeneous network.
#              The species-gene correlation analysis is parallelized.

# Load necessary libraries
library(here) # For portable path management
library(dplyr)
library(purrr)
library(readr)
library(foreach)
if (!require("doSNOW", quietly = TRUE)) {
    install.packages("doSNOW")
}
library(doSNOW)
library(tidyverse)
if (!require("igraph", quietly = TRUE)) {
    install.packages("igraph")
}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("DESeq2", quietly = TRUE)) {
    BiocManager::install("DESeq2")
}
if (!require("pheatmap", quietly = TRUE)) {
    install.packages("pheatmap")
}
if (!require("progress", quietly = TRUE)) {
    install.packages("progress")
}
library(igraph)
library(DESeq2)
library(pheatmap)
library(progress)

# 1. Define the directory containing your STAR count matrices
#    Replace this with the actual path to your 'covid2021_STAR' folder.
# star_counts_dir <- here::here("data", "covid2021_STAR")

# b1 <- read_tsv(here("data","covid2021_STAR","b1countMat.txt"))
# b2 <- read_tsv(here("data","covid2021_STAR","b2countMat.txt"))
# b3 <- read_tsv(here("data","covid2021_STAR","b3countMat.txt"))
# b4 <- read_tsv(here("data","covid2021_STAR","b4countMat.txt"))
# b5 <- read_tsv(here("data","covid2021_STAR","b5countMat.txt"))

# # 4. Join all the individual sample data frames into a single matrix.
# #    The `reduce` function iteratively joins the data frames by the 'gene_id' column.
# combined_counts_matrix <- purrr::reduce(list(b1, b2, b3, b4, b5), full_join, by = "Geneid")

# write.csv(combined_counts_matrix, here("data", "processed", "combined_star_gene_counts.csv"), row.names = FALSE)

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


# 1. Load Data (Assumes variables from your snippet exist)
# Clean Gene Count Columns: Remove "_count" suffix
colnames(count_mat) <- gsub("_count", "", colnames(count_mat))
rownames(count_mat) <- count_mat$Geneid
gene_mat <- count_mat %>% select(-"Geneid", -"Neg_Control_W") # Remove ID and Control

# Clean Kraken Columns: Remove "_cladeReads" suffix to match gene samples
# We need to ensure we are only matching the 'cladeReads' columns
kraken_mat <- species_mat_wide_format %>%
  select("name", ends_with("_cladeReads")) 
colnames(kraken_mat) <- gsub("_cladeReads", "", colnames(kraken_mat))
rownames(kraken_mat) <- kraken_mat$name
kraken_mat <- kraken_mat %>% select(-"name", -starts_with("Neg_Control"))

# Align Samples: Intersect sample names
common_samples <- intersect(colnames(gene_mat), colnames(kraken_mat))
gene_mat <- gene_mat[, common_samples]
kraken_mat <- kraken_mat[, common_samples]

print(paste("Aligned", length(common_samples), "samples."))

# Phase 2: Differentiating Species (The "High/Low" Split) ------------------------------------
# 1. Select Top Species (Filter by Freq and Abundance)
# Adjust quantiles as needed to control "Top" definition
top_species_df <- species_summary %>%
  filter(Freq > quantile(Freq, 0.90, na.rm = TRUE), 
         cladeReads_mean > quantile(cladeReads_mean, 0.70, na.rm = TRUE)) %>%
  arrange(desc(Freq))

target_species <- top_species_df$name
print(paste("Selected", length(target_species), "species for analysis."))


# 2. Loop through Species to find Gene Correlations (PARALLELIZED)
# Setup parallel backend
num_cores <- detectCores() # Use all but 2 cores.
cl <- makeCluster(num_cores)
registerDoSNOW(cl)

print(paste("Running gene correlation analysis for", length(target_species), "species on", num_cores, "cores..."))

pb <- progress_bar$new(
    format = "  Processing species [:bar] :percent eta: :eta",
    total = length(target_species), clear = FALSE, width = 60)
opts <- list(progress = function(n) pb$tick())

species_gene_edges <- foreach(i = seq_along(target_species), .combine = 'rbind', .options.snow = opts) %dopar% {
  
  sp <- target_species[i]
  
  # Fetch reads and handle NAs immediately
  sp_reads <- as.numeric(kraken_mat[sp, ])
  sp_reads[is.na(sp_reads)] <- 0 # CRITICAL FIX: Treat NAs as 0 counts
  
  result <- NULL # Initialize result for this iteration
  
  # Skip if species is effectively absent (present in fewer than 3 samples)
  if(sum(sp_reads > 0) >= 3) {
    
    # Define High vs Low/Absent groups
    if(mean(sp_reads > 0) > 0.20) {
      cutoff <- median(sp_reads[sp_reads > 0])
      condition <- ifelse(sp_reads >= cutoff, "High", "Low")
    } else {
      condition <- ifelse(sp_reads > 0, "High", "Low")
    }
    condition <- factor(condition, levels = c("Low", "High"))
    
    # Skip if we don't have at least 2 samples in each group (statistical requirement)
    if(min(table(condition)) >= 2) {
      
      # Run Wilcoxon Test (Non-parametric t-test equivalent)
      p_values <- suppressWarnings(
        apply(gene_mat, 1, function(g) {
          tryCatch(wilcox.test(g ~ condition)$p.value, error = function(e) NA)
        })
      )
      
      # Filter for significant genes (p < 0.01)
      sig_genes <- names(p_values)[!is.na(p_values) & p_values < 0.01]
      
      # Create Edges if genes found
      if(length(sig_genes) > 0){
        
        # Calculate Log2 Fold Change
        sub_mat <- gene_mat[sig_genes, , drop=FALSE] 
        means <- t(apply(sub_mat, 1, function(g) tapply(g, condition, mean, na.rm=TRUE)))
        
        if(is.null(dim(means))) {
            means <- matrix(means, nrow=1)
            colnames(means) <- c("Low", "High")
            rownames(means) <- sig_genes
        }
        
        logFC <- log2(means[, "High"] + 1) - log2(means[, "Low"] + 1)
        
        # This becomes the returned result for this iteration
        result <- data.frame(
          from = sp,
          to = sig_genes,
          type = "Species_Gene",
          weight = abs(logFC),
          direction = ifelse(logFC > 0, "Positive", "Negative"),
          p_value = p_values[sig_genes]
        )
      }
    }
  }
  
  result # Return either NULL or the new_edges data.frame
}

print(paste("Generated", nrow(species_gene_edges), "species-gene edges."))
write.csv(species_gene_edges, here("data", "processed", "species_gene_edges.csv"), row.names = FALSE)


# Phase 3: Gene-Gene Edges via Clustering
print("Phase 3: Gene-Gene Edges via Clustering")
# 1. Subset Gene Matrix
relevant_genes <- unique(species_gene_edges$to)

if(length(relevant_genes) > 5) {
  sub_gene_mat <- gene_mat[relevant_genes, ]
  
  # 2. Hierarchical Clustering
  print("Performing hierarchical clustering of genes...")
  dist_mat <- as.dist(1 - cor(t(sub_gene_mat), method = "spearman"))
  hc <- hclust(dist_mat, method = "ward.D2")
  
  # Dynamic tree cut: Split into modules
  gene_clusters <- cutree(hc, k = min(15, length(relevant_genes)-1)) 
  
  # 3. Create Gene-Gene Edges (Within-Cluster)
  unique_clusters <- unique(gene_clusters)
  
  print("Calculating gene-gene correlations within clusters...")
  pb <- progress_bar$new(
    format = "  Processing cluster [:bar] :percent eta: :eta",
    total = length(unique_clusters), clear = FALSE, width = 60)
  opts <- list(progress = function(n) pb$tick())
  
  gene_gene_edges <- foreach(k = unique_clusters, .combine = 'rbind', .packages = c("dplyr"), .options.snow = opts) %dopar% {
    cluster_genes <- names(gene_clusters)[gene_clusters == k]
    if (length(cluster_genes) < 2) {
      NULL
    } else {
      cluster_cor <- cor(t(sub_gene_mat[cluster_genes, ]), method = "spearman")
      cluster_cor[lower.tri(cluster_cor, diag = TRUE)] <- NA
      
      cluster_edges <- as.data.frame(as.table(cluster_cor)) %>%
        filter(!is.na(.data$Freq), abs(.data$Freq) > 0.7) %>%
        rename(from = .data$Var1, to = .data$Var2, weight = .data$Freq) %>%
        mutate(
          type = "Gene_Gene",
          direction = ifelse(.data$weight > 0, "Positive", "Negative"),
          p_value = NA
        )
      cluster_edges
    }
  }
  if (is.null(gene_gene_edges)) {
    gene_gene_edges <- data.frame()
  }
} else {
  print("Not enough genes for clustering.")
  gene_gene_edges <- data.frame()
}

# Stop the cluster after all parallel tasks are done
stopCluster(cl)

# Phase 4: Constructing the igraph Object
# Create Species-Species Edges
sp_cor <- cor(t(kraken_mat[target_species, ]), method="spearman")
sp_cor[lower.tri(sp_cor, diag=TRUE)] <- NA

sp_sp_edges <- as.data.frame(as.table(sp_cor)) %>%
  filter(!is.na(Freq), abs(Freq) > 0.5) %>%
  rename(from = Var1, to = Var2, weight = Freq) %>%
  mutate(type = "Species_Species", 
         direction = "Positive",
         p_value = NA)

# Combine ALL Edges
all_edges <- rbind(
  species_gene_edges,
  gene_gene_edges,
  sp_sp_edges
)

# Create Graph
nodes <- data.frame(id = unique(c(all_edges$from, all_edges$to)))
nodes$type <- ifelse(nodes$id %in% target_species, "Species", "Gene")

net <- graph_from_data_frame(d = all_edges, vertices = nodes, directed = FALSE)

# Export for Cytoscape
write_graph(net, file = "heterogeneous_network.graphml", format = "graphml")
print("Graph saved as heterogeneous_network.graphml")
