## Arguments
DESeq2_processed <- FALSE 

## Import
# Load necessary libraries
library(here)
library(dplyr)
library(purrr)
library(readr)
library(tidyverse)
library(igraph)
library(DESeq2)
library(pheatmap)
library(foreach)
library(doParallel)

cat("=== LOADING DATA ===\n")

# Load combined count matrix
cat("Loading combined gene counts...\n")
count_mat <- as.data.frame(read.csv(here("data","processed","combined_star_gene_counts.csv")))
cat("  Initial dimensions:", nrow(count_mat), "genes x", ncol(count_mat), "samples\n")

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
cat("  Sample names:", paste(head(common_samples, 3), collapse=", "), "...\n")

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
  head(10)  # Only top 10 species

cat("\nTop 10 species selected:\n")
print(top_species_df[, c("name", "taxID", "Freq", "cladeReads_mean")])

target_species <- top_species_df$name
cat("\nSpecies for analysis:", length(target_species), "\n")
for(i in seq_along(target_species)) {
  cat("  ", i, ": ", target_species[i], " (taxID: ", top_species_df$taxID[i], ")\n", sep="")
}

# Create output directory
dir.create(here("outputs", "DESeq2_test"), showWarnings = FALSE, recursive = TRUE)

cat("\n=== PHASE 2b: DESeq2 DIFFERENTIAL EXPRESSION ===\n")

if(!DESeq2_processed) {
  # Setup parallel processing
  num_cores <- 10
  cat("Setting up parallel processing with", num_cores, "cores...\n")
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Export necessary objects to cluster
  clusterExport(cl, c("gene_mat", "kraken_mat", "target_species", "here"))
  
  cat("Running DESeq2 for", length(target_species), "species...\n\n")
  
  species_gene_edges_list <- foreach(
    i = seq_along(target_species),
    .packages = c('DESeq2', 'dplyr', 'here'),
    .errorhandling = 'remove'
  ) %dopar% {
    
    sp <- target_species[i]
    
    # Fetch reads and handle NAs
    sp_reads <- as.numeric(kraken_mat[sp, ])
    sp_reads[is.na(sp_reads)] <- 0
    
    # Count presence
    n_present <- sum(sp_reads > 0)
    
    # Skip if insufficient presence
    if(n_present < 3) {
      return(list(
        species=sp, 
        status="skipped", 
        reason="<3 samples", 
        n_present=n_present,
        log_msg=paste("Species:", sp, "; Status: Skipped (<3 samples present)")
      ))
    }
    
    # Define High vs Low groups
    if(mean(sp_reads > 0) > 0.20) {
      cutoff <- median(sp_reads[sp_reads > 0])
      condition <- ifelse(sp_reads >= cutoff, "High", "Low")
      split_method <- "median"
    } else {
      condition <- ifelse(sp_reads > 0, "High", "Low")
      split_method <- "presence"
    }
    condition <- factor(condition, levels = c("Low", "High"))
    
    # Check group sizes
    group_sizes <- table(condition)
    if(min(group_sizes) < 2) {
      return(list(
        species=sp, 
        status="skipped", 
        reason="<2 in group", 
        n_present=n_present, 
        groups=as.list(group_sizes),
        log_msg=paste("Species:", sp, "; Status: Skipped (<2 samples in one group)")
      ))
    }
    
    # Prepare count matrix
    count_matrix <- round(gene_mat)
    colData <- data.frame(condition = condition, row.names = colnames(count_matrix))
    
    # Run DESeq2
    result <- tryCatch({
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
      
      # Save results
      results_filename <- paste0(gsub(" ", "_", sp), "_DESeq2_results.csv")
      write.csv(res_df, here("outputs", "DESeq2_test", results_filename), row.names = FALSE)
      
      # Count results at different thresholds
      n_total <- nrow(res_df)
      n_no_na <- sum(!is.na(res_df$padj))
      n_padj_05 <- sum(!is.na(res_df$padj) & res_df$padj < 0.05)
      n_l2fc_2 <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 2)
      
      # Apply final filtering (LOC and LINC already removed from count_mat)
      sig_genes_df <- res_df %>%
        dplyr::filter(
          !is.na(padj),
          padj < 0.05,
          abs(log2FoldChange) >= 2
        )
      
      n_final <- nrow(sig_genes_df)
      
      # Create log message
      log_msg <- paste0(
        "Species: ", sp, "; ",
        "Filtering Results: ",
        "Total=", n_total, ", ",
        "padj<0.05=", n_padj_05, " (", round(100*n_padj_05/n_total, 1), "%), ",
        "padj<0.05&|L2FC|>=2=", n_l2fc_2, " (", round(100*n_l2fc_2/n_total, 1), "%), ",
        "Final=", n_final
      )
      
      # Create edges if genes found
      if(n_final > 0) {
        new_edges <- sig_genes_df %>%
          dplyr::transmute(
            from = sp,
            to = gene,
            type = "Species_Gene",
            weight = abs(log2FoldChange),
            direction = ifelse(log2FoldChange > 0, "Positive", "Negative"),
            p_value = pvalue,
            padj = padj,
            logFC = log2FoldChange
          )
        
        list(
          species=sp, 
          status="success", 
          n_present=n_present,
          split_method=split_method,
          groups=as.list(group_sizes),
          n_total=n_total,
          n_padj_05=n_padj_05,
          n_l2fc_2=n_l2fc_2,
          n_final=n_final,
          edges=new_edges,
          log_msg=log_msg
        )
      } else {
        list(
          species=sp, 
          status="no_edges", 
          n_present=n_present,
          split_method=split_method,
          groups=as.list(group_sizes),
          n_total=n_total,
          n_padj_05=n_padj_05,
          n_l2fc_2=n_l2fc_2,
          n_final=0,
          log_msg=log_msg
        )
      }
    }, error = function(e) {
      list(
        species=sp, 
        status="error", 
        error_message=as.character(e$message),
        log_msg=paste("Species:", sp, "; Status: Error -", as.character(e$message))
      )
    })
    
    return(result)
  }
  
  stopCluster(cl)
  
  cat("\n=== DESEQ2 RESULTS ===\n")
  # Print log messages from each species
  for(i in seq_along(species_gene_edges_list)) {
    result <- species_gene_edges_list[[i]]
    if(!is.null(result$log_msg)) {
      cat(result$log_msg, "\n")
    }
  }
  
  # Extract edges from results
  cat("\n=== EXTRACTING EDGES ===\n")
  edge_dfs <- list()
  for(i in seq_along(species_gene_edges_list)) {
    result <- species_gene_edges_list[[i]]
    if(!is.null(result$edges)) {
      edge_dfs[[length(edge_dfs) + 1]] <- result$edges
      cat("Added", nrow(result$edges), "edges from", result$species, "\n")
    }
  }
  
  if(length(edge_dfs) > 0) {
    species_gene_edges <- do.call(rbind, edge_dfs)
  } else {
    species_gene_edges <- data.frame(
      from = character(),
      to = character(),
      type = character(),
      weight = numeric(),
      direction = character(),
      p_value = numeric(),
      padj = numeric(),
      logFC = numeric()
    )
  }
  
  cat("\nTotal species-gene edges:", nrow(species_gene_edges), "\n")
  
  write.csv(species_gene_edges, here("data", "processed", "test_species_gene_edges.csv"), row.names = FALSE)
  
} else {
  cat("Loading pre-computed species-gene edges...\n")
  species_gene_edges <- read.csv(here("data", "processed", "test_species_gene_edges.csv"))
  cat("Loaded", nrow(species_gene_edges), "edges\n")
}

cat("\n=== BUILDING SPECIES-GENE NETWORK ===\n")

if(nrow(species_gene_edges) > 0) {
  species_gene_nodes <- data.frame(id = unique(c(species_gene_edges$from, species_gene_edges$to)))
  species_gene_nodes$node_class <- ifelse(species_gene_nodes$id %in% target_species, "Species", "Gene")
  
  cat("Nodes:", nrow(species_gene_nodes), "(", sum(species_gene_nodes$node_class == "Species"), 
      "species,", sum(species_gene_nodes$node_class == "Gene"), "genes)\n")
  
  species_gene_net <- graph_from_data_frame(
    d = species_gene_edges[, c("from", "to", "weight", "type", "direction", "p_value")],
    vertices = species_gene_nodes,
    directed = FALSE
  )
  
  write_graph(species_gene_net, file = "test_species_gene_network.graphml", format = "graphml")
  cat("Species-Gene network saved\n")
} else {
  cat("No edges found - skipping network creation\n")
}

cat("\n=== TEST COMPLETE ===\n")
