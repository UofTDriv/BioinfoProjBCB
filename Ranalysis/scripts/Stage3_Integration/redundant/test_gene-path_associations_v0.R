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

# Load only batch 1 for testing
cat("Loading batch 1 gene counts...\n")
b1 <- read_tsv(here("data","covid2021_STAR","b1countMat.txt"))
cat("  Batch 1 dimensions:", nrow(b1), "genes x", ncol(b1), "samples\n")

count_mat <- as.data.frame(b1)
rm(b1)

# Load species data
cat("Loading species data...\n")
species_mat_wide_format <- read.csv(here("outputs", "new_samples", "unaligned_merged260120op.csv"))
cat("  Species matrix dimensions:", nrow(species_mat_wide_format), "x", ncol(species_mat_wide_format), "\n")

species_summary <- read.csv(here("outputs", "new_samples", "species_list_unaligned260120op.csv"))
cat("  Species summary:", nrow(species_summary), "species\n")

cat("\n=== PHASE 1: DATA CLEANING & ALIGNMENT ===\n")

# Clean Gene Count Columns
cat("Cleaning gene count column names...\n")
colnames(count_mat) <- gsub("_count", "", colnames(count_mat))
rownames(count_mat) <- count_mat$Geneid
gene_mat <- count_mat %>% select(-Geneid, -Neg_Control_W)
cat("  Gene matrix:", nrow(gene_mat), "genes x", ncol(gene_mat), "samples\n")

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

# Select only top 10 species for testing
cat("Calculating species quantiles...\n")
freq_90 <- quantile(species_summary$Freq, 0.90, na.rm = TRUE)
reads_70 <- quantile(species_summary$cladeReads_mean, 0.70, na.rm = TRUE)
cat("  Freq 90th percentile:", freq_90, "\n")
cat("  cladeReads 70th percentile:", reads_70, "\n")

top_species_df <- species_summary %>%
  filter(Freq > freq_90, cladeReads_mean > reads_70) %>%
  arrange(desc(Freq)) %>%
  head(10)

cat("\nTop 10 species:\n")
print(top_species_df[, c("name", "taxID", "Freq", "cladeReads_mean")])

# Analyze only the TOP 1 species
target_species <- top_species_df$name[1]
cat("\n=== ANALYZING ONLY TOP SPECIES ===\n")
cat("Species:", target_species, "\n")
cat("TaxID:", top_species_df$taxID[1], "\n")
cat("Frequency:", top_species_df$Freq[1], "\n")
cat("Mean cladeReads:", top_species_df$cladeReads_mean[1], "\n")

# Create output directory
dir.create(here("outputs", "DESeq2_test"), showWarnings = FALSE, recursive = TRUE)

cat("\n=== PHASE 2b: DESeq2 DIFFERENTIAL EXPRESSION (NO PARALLELIZATION) ===\n")

sp <- target_species

# Fetch reads and handle NAs
cat("Fetching species reads for:", sp, "\n")
sp_reads <- as.numeric(kraken_mat[sp, ])
cat("  Raw reads (first 10 samples):", paste(head(sp_reads, 10), collapse=", "), "\n")

sp_reads[is.na(sp_reads)] <- 0
cat("  After NA replacement (first 10):", paste(head(sp_reads, 10), collapse=", "), "\n")

n_present <- sum(sp_reads > 0)
cat("  Samples with species present:", n_present, "/", length(sp_reads), "\n")
cat("  Min reads:", min(sp_reads), "\n")
cat("  Max reads:", max(sp_reads), "\n")
cat("  Median reads:", median(sp_reads), "\n")
cat("  Mean reads:", mean(sp_reads), "\n")

# Check if sufficient presence
if(n_present < 3) {
  stop("ERROR: Species present in fewer than 3 samples. Cannot proceed.")
}

# Define High vs Low groups
cat("\nDefining sample groups...\n")
if(mean(sp_reads > 0) > 0.20) {
  cutoff <- median(sp_reads[sp_reads > 0])
  condition <- ifelse(sp_reads >= cutoff, "High", "Low")
  split_method <- "median"
  cat("  Method: Median split (species present in >20% of samples)\n")
  cat("  Cutoff:", cutoff, "\n")
} else {
  condition <- ifelse(sp_reads > 0, "High", "Low")
  split_method <- "presence/absence"
  cat("  Method: Presence/Absence split\n")
}
condition <- factor(condition, levels = c("Low", "High"))

# Check group sizes
group_sizes <- table(condition)
cat("  Group sizes:\n")
cat("    Low:", group_sizes["Low"], "\n")
cat("    High:", group_sizes["High"], "\n")

if(min(group_sizes) < 2) {
  stop("ERROR: At least one group has fewer than 2 samples. Cannot proceed.")
}

# Show which samples are in which group
cat("\n  Sample assignments:\n")
for(i in 1:min(10, length(condition))) {
  cat("    ", common_samples[i], ": ", condition[i], " (reads: ", sp_reads[i], ")\n", sep="")
}
if(length(condition) > 10) {
  cat("    ... and", length(condition) - 10, "more samples\n")
}

# Prepare count matrix
cat("\nPreparing count matrix for DESeq2...\n")
count_matrix <- round(gene_mat)
cat("  Count matrix dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")
cat("  Min count:", min(count_matrix), "\n")
cat("  Max count:", max(count_matrix), "\n")

colData <- data.frame(condition = condition, row.names = colnames(count_matrix))
cat("  Sample metadata created\n")

# Run DESeq2
cat("\nRunning DESeq2 analysis...\n")
cat("  Creating DESeqDataSet...\n")
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = colData,
  design = ~ condition
)
cat("  DESeqDataSet created successfully\n")

cat("  Running DESeq (this may take a few minutes)...\n")
dds <- DESeq(dds, quiet = FALSE)
cat("  DESeq completed\n")

cat("  Extracting results...\n")
res <- results(dds, contrast = c("condition", "High", "Low"))
cat("  Results extracted\n")

# Convert to data frame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Save results
results_filename <- paste0(gsub(" ", "_", sp), "_DESeq2_results.csv")
write.csv(res_df, here("outputs", "DESeq2_test", results_filename), row.names = FALSE)
cat("  Results saved to:", results_filename, "\n")

# Count results at different thresholds
cat("\n=== FILTERING STATISTICS ===\n")
n_total <- nrow(res_df)
cat("Total genes tested:", n_total, "\n")

n_no_na <- sum(!is.na(res_df$padj))
cat("Non-NA padj:", n_no_na, "(", round(100*n_no_na/n_total, 1), "%)\n")

n_padj_05 <- sum(!is.na(res_df$padj) & res_df$padj < 0.05)
cat("padj < 0.05:", n_padj_05, "(", round(100*n_padj_05/n_total, 1), "%)\n")

n_padj_01 <- sum(!is.na(res_df$padj) & res_df$padj < 0.01)
cat("padj < 0.01:", n_padj_01, "(", round(100*n_padj_01/n_total, 1), "%)\n")

n_l2fc_1 <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1)
cat("padj < 0.05 & |L2FC| >= 1:", n_l2fc_1, "(", round(100*n_l2fc_1/n_total, 1), "%)\n")

n_l2fc_2 <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 2)
cat("padj < 0.05 & |L2FC| >= 2:", n_l2fc_2, "(", round(100*n_l2fc_2/n_total, 1), "%)\n")

# Show top significant genes before LOC filter
cat("\nTop 10 significant genes (padj < 0.05, sorted by padj):\n")
top_sig <- res_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj) %>%
  head(10)
if(nrow(top_sig) > 0) {
  print(top_sig[, c("gene", "log2FoldChange", "padj")])
} else {
  cat("  (none found)\n")
}

# Apply strict filtering
cat("\n=== APPLYING FINAL FILTERS ===\n")
cat("Filters:\n")
cat("  1. Remove NA padj\n")
cat("  2. padj < 0.05\n")
cat("  3. |log2FoldChange| >= 2\n")
cat("  4. Remove genes starting with 'LOC'\n")

sig_genes_df <- res_df %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) >= 2,
    !grepl("^LOC", gene)
  )

n_final <- nrow(sig_genes_df)
cat("\nFinal significant genes:", n_final, "\n")

if(n_final > 0) {
  cat("\nSignificant genes:\n")
  print(sig_genes_df[, c("gene", "log2FoldChange", "padj")])
  
  # Create edges
  cat("\n=== CREATING EDGES ===\n")
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
  
  cat("Created", nrow(species_gene_edges), "edges\n")
  
  write.csv(species_gene_edges, here("data", "processed", "test_species_gene_edges_v0.csv"), row.names = FALSE)
  cat("Edges saved to: test_species_gene_edges_v0.csv\n")
  
  # Build network
  cat("\n=== BUILDING NETWORK ===\n")
  species_gene_nodes <- data.frame(id = unique(c(species_gene_edges$from, species_gene_edges$to)))
  species_gene_nodes$node_class <- ifelse(species_gene_nodes$id == sp, "Species", "Gene")
  
  cat("Nodes:", nrow(species_gene_nodes), "(1 species,", nrow(species_gene_nodes)-1, "genes)\n")
  
  species_gene_net <- graph_from_data_frame(
    d = species_gene_edges[, c("from", "to", "weight", "type", "direction", "p_value")],
    vertices = species_gene_nodes,
    directed = FALSE
  )
  
  write_graph(species_gene_net, file = "test_species_gene_network_v0.graphml", format = "graphml")
  cat("Network saved to: test_species_gene_network_v0.graphml\n")
  
} else {
  cat("\nWARNING: No genes passed all filtering criteria!\n")
  cat("Suggestions:\n")
  cat("  - Try relaxing |L2FC| threshold from 2 to 1\n")
  cat("  - Try relaxing padj threshold from 0.05 to 0.1\n")
  cat("  - Check if species groups are well-separated\n")
}

cat("\n=== TEST COMPLETE ===\n")
