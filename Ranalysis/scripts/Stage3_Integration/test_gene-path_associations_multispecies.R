## Arguments
DESeq2_processed <- FALSE 

## Import
# Load necessary libraries
library(here)
library(dplyr)
library(readr)
library(tidyverse)
library(igraph)
library(DESeq2)

cat("=== LOADING DATA ===\n")

count_mat <- as.data.frame(read.csv(here("data","processed","combined_star_gene_counts.csv")))

# Load species data
cat("Loading species data...\n")
species_mat_wide_format <- read.csv(here("outputs", "new_samples", "unaligned_merged260120op.csv"))
cat("  Species matrix dimensions:", nrow(species_mat_wide_format), "x", ncol(species_mat_wide_format), "\n")

species_summary <- read.csv(here("outputs", "new_samples", "species_list_unaligned260120op.csv"))
cat("  Species summary:", nrow(species_summary), "species\n")

cat("\n=== PHASE 1: DATA CLEANING & FILTERING ===\n")

# Clean Gene Count Columns
cat("Cleaning gene count column names...\n")
colnames(count_mat) <- gsub("_count", "", colnames(count_mat))
cat("  Initial gene count:", nrow(count_mat), "\n")

# Step 1: Remove rRNA genes (if rRNA list exists)
rrna_file <- here("data", "genes", "hgnc_rRNAgenelist.csv")
if(file.exists(rrna_file)) {
  cat("Removing rRNA genes...\n")
  rrnalist <- read.csv(rrna_file)
  count_mat <- subset(count_mat, !(Geneid %in% rrnalist$Approved.symbol))
  cat("  After rRNA removal:", nrow(count_mat), "genes\n")
} else {
  cat("  No rRNA gene list found, skipping rRNA filtering\n")
}

# Step 2: Drop rows which are all 0
cat("Removing zero-count genes...\n")
count_mat <- count_mat %>% 
  rowwise() %>% 
  mutate(total = sum(c_across(where(is.numeric)))) %>%
  filter(total > 0) %>%
  select(-total) %>%
  ungroup()
count_mat <- as.data.frame(count_mat)
cat("  After zero removal:", nrow(count_mat), "genes\n")

# Step 3: Merge alternative splice genes (e.g., GENE-AS1, GENE-AS2)
cat("Merging alternative splice variants (-AS#)...\n")
count_mat$Geneid_trim <- gsub("-AS[0-9]", "", count_mat$Geneid)
count_mat <- count_mat %>% 
  select(-Geneid) %>%
  group_by(Geneid_trim) %>%
  summarize_all(sum) %>%
  ungroup()
count_mat <- as.data.frame(count_mat)
cat("  After AS merging:", nrow(count_mat), "genes\n")

# Step 4: Merge duplicate gene variants (e.g., GENE_1, GENE_2)
cat("Merging duplicate gene variants (_#)...\n")
count_mat$Geneid_trim2 <- gsub("_[0-9]", "", count_mat$Geneid_trim)
count_mat <- count_mat %>%
  select(-Geneid_trim) %>%
  group_by(Geneid_trim2) %>%
  summarize_all(sum) %>%
  ungroup()
count_mat <- as.data.frame(count_mat)
colnames(count_mat)[1] <- "Geneid"
cat("  After duplicate merging:", nrow(count_mat), "genes\n")

# Step 5: Remove LOC genes
cat("Removing LOC genes...\n")
count_mat <- count_mat %>% filter(!grepl("^LOC[0-9]*", Geneid))
cat("  After LOC removal:", nrow(count_mat), "genes\n")

# Step 6: Remove LINC genes
cat("Removing LINC genes...\n")
count_mat <- count_mat %>% filter(!grepl("^LINC[0-9]*", Geneid))
cat("  After LINC removal:", nrow(count_mat), "genes\n")

# Set rownames and create gene matrix
rownames(count_mat) <- count_mat$Geneid
gene_mat <- count_mat %>% select(-Geneid, -Neg_Control_W)
cat("  Final gene matrix:", nrow(gene_mat), "genes x", ncol(gene_mat), "samples\n")

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

# Replace all NAs in kraken_mat with 0
kraken_mat[is.na(kraken_mat)] <- 0
cat("  Replaced NAs in kraken_mat with 0\n")

cat("  Aligned", length(common_samples), "samples between gene and species data\n")
cat("  Sample names:", paste(head(common_samples, 5), collapse=", "), "...\n")

cat("\n=== PHASE 2: SPECIES SELECTION ===\n")

# Select top 10 species for testing
cat("Calculating species quantiles...\n")
freq_90 <- quantile(species_summary$Freq, 0.90, na.rm = TRUE)
reads_70 <- quantile(species_summary$cladeReads_mean, 0.70, na.rm = TRUE)
cat("  Freq 90th percentile:", freq_90, "\n")
cat("  cladeReads 70th percentile:", reads_70, "\n")

top_species_df <- species_summary %>%
  filter(Freq > freq_90, cladeReads_mean > reads_70) %>%
  arrange(desc(Freq)) %>%
  head(10)

cat("\nTop 10 species selected for analysis:\n")
print(top_species_df[, c("name", "taxID", "Freq", "cladeReads_mean")])

target_species_list <- top_species_df$name

# Create output directory
dir.create(here("outputs", "DESeq2_test"), showWarnings = FALSE, recursive = TRUE)

cat("\n=== PHASE 2b: DESeq2 DIFFERENTIAL EXPRESSION (SEQUENTIAL) ===\n")

# Initialize results container
all_results <- list()
species_gene_edges_list <- list()

# Loop through each species
for(species_idx in seq_along(target_species_list)) {
  sp <- target_species_list[species_idx]
  
  cat("SPECIES", species_idx, "OF", length(target_species_list), ":", sp, "\n")
  
  # Fetch reads and handle NAs
  cat("Fetching species reads...\n")
  sp_reads <- as.numeric(kraken_mat[sp, ])
  sp_reads[is.na(sp_reads)] <- 0
  
  n_present <- sum(sp_reads > 0)
  cat("  Samples with species present:", n_present, "/", length(sp_reads), "\n")
  cat("  Min reads:", min(sp_reads), ", Max:", max(sp_reads), ", Mean:", round(mean(sp_reads), 2), "\n")
  
  # Check if sufficient presence
  if(n_present < 3) {
    cat("  SKIPPED: Species present in fewer than 3 samples\n")
    next
  }
  
  # Define High vs Low groups
  cat("Defining sample groups...\n")
  if(mean(sp_reads > 0) > 0.20) {
    cutoff <- median(sp_reads[sp_reads > 0])
    condition <- ifelse(sp_reads >= cutoff, "High", "Low")
    split_method <- "median"
    cat("  Method: Median split (cutoff:", round(cutoff, 2), ")\n")
  } else {
    condition <- ifelse(sp_reads > 0, "High", "Low")
    split_method <- "presence/absence"
    cat("  Method: Presence/Absence split\n")
  }
  condition <- factor(condition, levels = c("Low", "High"))
  
  # Check group sizes
  group_sizes <- table(condition)
  cat("  Group sizes - Low:", group_sizes["Low"], ", High:", group_sizes["High"], "\n")
  
  if(min(group_sizes) < 2) {
    cat("  SKIPPED: At least one group has fewer than 2 samples\n")
    next
  }
  
  # Prepare count matrix
  cat("Running DESeq2 analysis...\n")
  count_matrix <- round(gene_mat)
  colData <- data.frame(condition = condition, row.names = colnames(count_matrix))
  
  # Run DESeq2 with error handling
  tryCatch({
    dds <- DESeqDataSetFromMatrix(
      countData = count_matrix,
      colData = colData,
      design = ~ condition
    )
    
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("condition", "High", "Low"))
    
    # Convert to data frame
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    
    # Save full results
    results_filename <- paste0(gsub(" ", "_", sp), "_DESeq2_results.csv")
    write.csv(res_df, here("outputs", "DESeq2_test", results_filename), row.names = FALSE)
    
    # Count results at different thresholds
    n_total <- nrow(res_df)
    n_no_na <- sum(!is.na(res_df$padj))
    n_padj_05 <- sum(!is.na(res_df$padj) & res_df$padj < 0.05)
    n_padj_01 <- sum(!is.na(res_df$padj) & res_df$padj < 0.01)
    n_l2fc_1 <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1)
    n_l2fc_2 <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 2)
    
    cat("  Filtering results:\n")
    cat("    Total genes:", n_total, "\n")
    cat("    padj < 0.05:", n_padj_05, "(", round(100*n_padj_05/n_total, 1), "%)\n")
    cat("    padj < 0.05 & |L2FC| >= 1:", n_l2fc_1, "(", round(100*n_l2fc_1/n_total, 1), "%)\n")
    cat("    padj < 0.05 & |L2FC| >= 2:", n_l2fc_2, "(", round(100*n_l2fc_2/n_total, 1), "%)\n")
    
    # Apply final filtering (LOC and LINC already removed from count_mat)
    sig_genes_df <- res_df %>%
      filter(
        !is.na(padj),
        padj < 0.05,
        abs(log2FoldChange) >= 2
      )
    
    n_final <- nrow(sig_genes_df)
    cat("    Final significant genes:", n_final, "\n")
    
    # Store result information
    result_info <- list(
      species = sp,
      status = ifelse(n_final > 0, "success", "no_edges"),
      n_present = n_present,
      split_method = split_method,
      n_low = as.numeric(group_sizes["Low"]),
      n_high = as.numeric(group_sizes["High"]),
      n_total = n_total,
      n_padj_05 = n_padj_05,
      n_l2fc_2 = n_l2fc_2,
      n_final = n_final
    )
    all_results[[length(all_results) + 1]] <- result_info
    
    # Create edges if genes found
    if(n_final > 0) {
      species_gene_edges <- sig_genes_df %>%
        transmute(
          from = sp,
          to = gene,
          type = "Species_Gene",
          weight = abs(log2FoldChange),
          direction = ifelse(log2FoldChange > 0, "Positive", "Negative"),
          p_value = pvalue,
          padj = padj,
          logFC = log2FoldChange
        )
      
      species_gene_edges_list[[length(species_gene_edges_list) + 1]] <- species_gene_edges
      cat("  SUCCESS: Created", nrow(species_gene_edges), "edges\n")
    } else {
      cat("  NO EDGES: No significant genes found\n")
    }
    
  }, error = function(e) {
    cat("  ERROR:", as.character(e), "\n")
  })
}

# Combine all edges
cat("=== SUMMARY ===\n")

# Print per-species results
cat("\nPer-species results:\n")
for(i in seq_along(all_results)) {
  r <- all_results[[i]]
  cat("\n", i, ".", r$species, "\n")
  cat("   Status:", r$status, "\n")
  cat("   Samples: Low=", r$n_low, ", High=", r$n_high, "\n")
  cat("   Genes: padj<0.05,|L2FC|>=2:", r$n_l2fc_2, " -> Final edges:", r$n_final, "\n")
}

if(length(species_gene_edges_list) > 0) {
  species_gene_edges_combined <- do.call(rbind, species_gene_edges_list)
  cat("\nTotal species-gene edges:", nrow(species_gene_edges_combined), "\n")
  
  write.csv(species_gene_edges_combined, here("data", "processed", "test_species_gene_edges.csv"), row.names = FALSE)
  cat("Edges saved to: test_species_gene_edges.csv\n")
  
  # Build network
  cat("\n=== BUILDING NETWORK ===\n")
  species_gene_nodes <- data.frame(id = unique(c(species_gene_edges_combined$from, species_gene_edges_combined$to)))
  species_gene_nodes$node_class <- ifelse(species_gene_nodes$id %in% target_species_list, "Species", "Gene")
  
  n_species_nodes <- sum(species_gene_nodes$node_class == "Species")
  n_gene_nodes <- sum(species_gene_nodes$node_class == "Gene")
  
  cat("Total nodes:", nrow(species_gene_nodes), "\n")
  cat("  Species nodes:", n_species_nodes, "\n")
  cat("  Gene nodes:", n_gene_nodes, "\n")
  
  species_gene_net <- graph_from_data_frame(
    d = species_gene_edges_combined[, c("from", "to", "weight", "type", "direction", "p_value")],
    vertices = species_gene_nodes,
    directed = FALSE
  )
  
  write_graph(species_gene_net, file = "test_species_gene_network.graphml", format = "graphml")
  cat("Network saved to: test_species_gene_network.graphml\n")
  
} else {
  cat("\nWARNING: No edges found across all species!\n")
  cat("Suggestion: Check if LOC genes are dominating your results.\n")
  cat("Consider relaxing the L2FC threshold or padj cutoff.\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
