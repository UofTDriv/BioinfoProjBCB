## Arguments
DESeq2_processed <- FALSE 

## Import
# Load necessary libraries
library(here) # For portable path management
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

# Paths
filtered_counts_path <- here("data", "processed", "filtered_gene_counts.csv")
combined_counts_path <- here("data", "processed", "combined_star_gene_counts.csv")
rrna_file <- here("data", "genes", "hgnc_rRNAgenelist.csv")

# Load or create filtered count matrix
if(file.exists(filtered_counts_path)) {
  cat("Loading filtered count matrix...\n")
  count_mat <- read.csv(filtered_counts_path, check.names = FALSE)
  # Handle possible unnamed first column
  if(!"Geneid" %in% colnames(count_mat)) {
    colnames(count_mat)[1] <- "Geneid"
  }
  cat("  Filtered gene count matrix:", nrow(count_mat), "genes x", ncol(count_mat)-1, "samples\n")
} else {
  cat("Filtered count matrix not found, creating from combined counts...\n")
  if(!file.exists(combined_counts_path)) {
    stop("Missing combined_star_gene_counts.csv. Generate it before running this script.")
  }
  count_mat <- as.data.frame(read.csv(combined_counts_path, check.names = FALSE))
  cat("  Combined gene count matrix:", nrow(count_mat), "genes x", ncol(count_mat)-1, "samples\n")

  # Clean Gene Count Columns
  cat("Cleaning gene count column names...\n")
  colnames(count_mat) <- gsub("_count", "", colnames(count_mat))
  cat("  Initial gene count:", nrow(count_mat), "\n")

  # Step 1: Remove rRNA genes (if rRNA list exists)
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

  # Save filtered count matrix
  write.csv(count_mat, filtered_counts_path, row.names = FALSE)
  cat("Filtered count matrix saved to:", filtered_counts_path, "\n")
}

# Set rownames and create gene matrix
rownames(count_mat) <- count_mat$Geneid
gene_mat <- count_mat %>% select(-Geneid, -Neg_Control_W)
cat("Final gene matrix:", nrow(gene_mat), "genes x", ncol(gene_mat), "samples\n")

# Load species data
cat("Loading species data...\n")
species_mat_wide_format <- read.csv(here("outputs", "new_samples", "unaligned_merged260120op.csv"))
cat("  Species matrix dimensions:", nrow(species_mat_wide_format), "x", ncol(species_mat_wide_format), "\n")

species_summary <- read.csv(here("outputs", "new_samples", "species_list_unaligned260120op.csv"))
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

# Create output directory for DESeq2 results if it doesn't exist
dir.create(here("outputs", "DESeq2"), showWarnings = FALSE, recursive = TRUE)

# 2. Loop through Species to find Gene Correlations with DESeq2
if(!DESeq2_processed) {
  # Setup parallel processing
  num_cores <- min(24, max(1, parallel::detectCores() - 1))
  cat("Setting up parallel processing with", num_cores, "cores...\n")
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  species_gene_edges_list <- foreach(
    i = seq_along(target_species),
    .packages = c('DESeq2', 'dplyr', 'here'),
    .errorhandling = 'remove'
  ) %dopar% {
    
    sp <- target_species[i]
    
    # Fetch reads and handle NAs
    sp_reads <- as.numeric(kraken_mat[sp, ])
    sp_reads[is.na(sp_reads)] <- 0
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
    
    # Prepare count matrix for DESeq2
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
      
      # Convert to data frame for filtering
      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)
      
      # Save full DESeq2 results for this species
      results_filename <- paste0(gsub(" ", "_", sp), "_DESeq2_results.csv")
      write.csv(res_df, here("outputs", "DESeq2", results_filename), row.names = FALSE)
      
      # Apply filtering criteria
      n_total <- nrow(res_df)
      n_padj_05 <- sum(!is.na(res_df$padj) & res_df$padj < 0.05)
      n_l2fc_2 <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 2)
      
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
      
      if(n_final > 0){
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
  for(i in seq_along(species_gene_edges_list)) {
    result <- species_gene_edges_list[[i]]
    if(!is.null(result$log_msg)) {
      cat(result$log_msg, "\n")
    }
  }

  # Extract edges from results
  edge_dfs <- list()
  for(i in seq_along(species_gene_edges_list)) {
    result <- species_gene_edges_list[[i]]
    if(!is.null(result$edges)) {
      edge_dfs[[length(edge_dfs) + 1]] <- result$edges
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

  cat("Generated", nrow(species_gene_edges), "species-gene edges.\n")
  write.csv(species_gene_edges, here("data", "processed", "species_gene_edges.csv"), row.names = FALSE)
} else {
  cat("DESeq2_processed is TRUE, loading pre-computed species-gene edges...\n")
  species_gene_edges <- read.csv(here("data", "processed", "species_gene_edges.csv"))
  cat("Loaded", nrow(species_gene_edges), "species-gene edges from file.\n")
}

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
