#!/usr/bin/env Rscript

# Stage3_GO_Microbe_Metagraph.R
# Builds a Cytoscape-ready metagraph by collapsing genes into positively enriched GO:BP metanodes
# and microbes into risk-group or taxonomic metanodes. Includes RISK-compatible annotation export.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(igraph)
  library(fgsea)
  library(msigdbr)
})

# Source utility functions
source(here("Ranalysis", "scripts", "Stage1_Data_acquisition", "utils.R"))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

normalize_label <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("[^a-z0-9]", "")
}

resolve_existing_path <- function(candidates) {
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop("None of the candidate paths exist:\n", paste(candidates, collapse = "\n"))
  }
  existing[[1]]
}

extract_tax_rank <- function(tax_lineage, rank_prefix) {
  pat <- paste0("(?<=\\|", rank_prefix, "_)[^|]+|(?<=^", rank_prefix, "_)[^|]+")
  out <- str_extract(tax_lineage %||% "", pat)
  out[is.na(out) | out == ""] <- paste0("Unknown_", toupper(rank_prefix))
  out
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  kv <- list()
  if (length(args) > 0) {
    for (a in args) {
      if (str_detect(a, "=")) {
        parts <- str_split_fixed(a, "=", 2)
        k <- str_replace(parts[, 1], "^--", "")
        v <- parts[, 2]
        kv[[k]] <- v
      }
    }
  }

  list(
    species_gene_edges = kv$species_gene_edges %||%
      here("data", "processed", "species_gene_edges.csv"),
    filtered_gene_counts = kv$filtered_gene_counts %||%
      here("data", "processed", "filtered_gene_counts.csv"),
    species_list = kv$species_list %||%
      get_latest_timestamped_file(input_dir = "data/processed", pattern = "^species_list_unaligned.*op.*\\.csv$"),
    epathogen = kv$epathogen %||%
      get_latest_timestamped_file(input_dir = "Ranalysis/databases", pattern = "^epathogen.*\\.csv$"),
    deseq_dir = kv$deseq_dir %||% resolve_existing_path(c(
      here("outputs", "DESeq2_test"),
      here("outputs", "DESeq2_results")
    )),
    microbe_grouping = tolower(kv$microbe_grouping %||% "risk"),
    out_graphml = kv$out_graphml %||%
      here("outputs", "stage3_go_microbe_metagraph.graphml"),
    out_overlap = kv$out_overlap %||%
      here("outputs", "stage3_cluster_overlap_jaccard.csv"),
    out_density = kv$out_density %||%
      here("outputs", "stage3_module_density_comparison.csv"),
    out_fgsea = kv$out_fgsea %||%
      here("outputs", "stage3_fgsea_go_bp_results.csv"),
    out_risk_csv = kv$out_risk_csv %||%
      here("outputs", "stage3_risk_go_annotations.csv")
  )
}

read_gene_counts <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)

  if ("Geneid" %in% colnames(x)) {
      x <- x %>% select(-Geneid)
  }

  if ("Neg_Control_W" %in% colnames(x)) {
    x <- x %>% select(-Neg_Control_W)
  }

  x
}

read_deseq_results <- function(deseq_dir) {
  files <- list.files(deseq_dir, pattern = "_DESeq2_results\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No DESeq2 result files found in: ", deseq_dir)
  }

  out <- purrr::map_dfr(files, function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    required <- c("gene", "log2FoldChange", "pvalue")
    if (!all(required %in% colnames(d))) {
      return(tibble())
    }
    sp <- basename(f) %>%
      str_remove("_DESeq2_results\\.csv$") %>%
      str_replace_all("_", " ")

    d %>%
      transmute(
        species = sp,
        species_norm = normalize_label(sp),
        gene = as.character(gene),
        log2FoldChange = as.numeric(log2FoldChange),
        pvalue = as.numeric(pvalue),
        padj = as.numeric(padj)
      )
  })

  if (nrow(out) == 0) {
    stop("DESeq2 files were found but none had required columns: gene, log2FoldChange, pvalue")
  }

  out
}

build_gene_ranking <- function(species_gene_edges, deseq_df) {
  edge_pairs <- species_gene_edges %>%
    transmute(
      species = as.character(from),
      species_norm = normalize_label(from),
      gene = as.character(to)
    ) %>%
    distinct()

  direct_join <- edge_pairs %>%
    left_join(deseq_df, by = c("species", "species_norm", "gene"))

  gene_fallback <- deseq_df %>%
    group_by(gene) %>%
    summarise(
      log2FoldChange_fb = ifelse(
        all(is.na(log2FoldChange)),
        NA_real_,
        median(log2FoldChange, na.rm = TRUE)
      ),
      pvalue_fb = ifelse(
        all(is.na(pvalue)),
        NA_real_,
        min(pvalue, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  scored <- direct_join %>%
    left_join(gene_fallback, by = "gene") %>%
    mutate(
      log2FoldChange_use = ifelse(is.finite(log2FoldChange), log2FoldChange, log2FoldChange_fb),
      pvalue_use = ifelse(is.finite(pvalue), pvalue, pvalue_fb),
      pvalue_use = pmax(pvalue_use, 1e-300),
      rank_score = sign(log2FoldChange_use) * -log10(pvalue_use)
    ) %>%
    filter(is.finite(rank_score), !is.na(gene))

  gene_rank_df <- scored %>%
    group_by(gene) %>%
    summarise(
      rank_score = mean(rank_score, na.rm = TRUE),
      mean_log2FoldChange = mean(log2FoldChange_use, na.rm = TRUE),
      min_pvalue = min(pvalue_use, na.rm = TRUE),
      n_species_links = n_distinct(species),
      .groups = "drop"
    ) %>%
    filter(is.finite(rank_score))

  ranks <- gene_rank_df$rank_score
  names(ranks) <- gene_rank_df$gene
  # Deterministic tie-breaker to reduce arbitrary ordering in fgsea.
  tie_breaker <- rank(names(ranks), ties.method = "first") * 1e-12
  ranks <- ranks + tie_breaker
  ranks <- sort(ranks, decreasing = TRUE)

  list(gene_rank_df = gene_rank_df, ranks = ranks)
}

run_fgsea_go_bp <- function(ranks) {
  msig <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>%
    select(gs_name, gene_symbol) %>%
    distinct()

  pathways <- split(msig$gene_symbol, msig$gs_name)
  pathways <- pathways[lengths(pathways) >= 10]

  if (length(pathways) == 0) {
    stop("No GO:BP pathways available after filtering min size.")
  }

  fg <- fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = 10,
    maxSize = 500
  ) %>%
    as_tibble() %>%
    filter(NES > 0) %>%  # FOCUS ON POSITIVELY UPREGULATED GENE SETS
    arrange(padj, desc(NES))

  if (nrow(fg) == 0) {
    stop("fgsea returned no positive GO terms.")
  }

  fg
}

build_go_gene_modules <- function(fgsea_tbl, genes_of_interest) {
  
  # 1. Filter for significant terms
  go_terms_use <- fgsea_tbl %>%
    dplyr::filter(!is.na(padj), padj < 0.1)

  # Fallback if none are significant
  if (nrow(go_terms_use) == 0) {
    go_terms_use <- fgsea_tbl %>%
      dplyr::slice_head(n = min(30, nrow(fgsea_tbl)))
  }

  # 2. Extract Leading Edge genes
  # Explicitly naming the column in unnest helps avoid ambiguity
  le <- go_terms_use %>%
    dplyr::select(pathway, padj, NES, leadingEdge) %>%
    tidyr::unnest(cols = c(leadingEdge)) %>% 
    dplyr::rename(gene = leadingEdge) %>%
    dplyr::mutate(gene = as.character(gene)) %>%
    dplyr::arrange(padj, dplyr::desc(NES)) %>%
    # If a gene appears in multiple pathways, pick the most significant one
    dplyr::group_by(gene) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(gene, go_term = pathway)

  # 3. Map back to your original gene list
  tibble::tibble(gene = genes_of_interest) %>%
    dplyr::left_join(le, by = "gene") %>%
    dplyr::mutate(go_term = tidyr::replace_na(go_term, "GO_BP_Unassigned"))
}

build_spearman_modules <- function(gene_mat, genes_of_interest) {
  genes <- intersect(genes_of_interest, rownames(gene_mat))
  print(paste("Number of intersecting genes:", length(genes)))
  if (length(genes) < 6) {
    stop("Need at least 6 genes with expression values to build Spearman modules.")
  }

  sub_gene_mat <- gene_mat
  cor_mat <- suppressWarnings(cor(t(sub_gene_mat), method = "spearman", use = "pairwise.complete.obs"))
  cor_mat[is.na(cor_mat)] <- 0

  dist_mat <- as.dist(1 - cor_mat)
  hc <- hclust(dist_mat, method = "ward.D2")
  k <- min(15, length(genes) - 1)
  clusters <- cutree(hc, k = k)

  membership <- tibble(
    gene = names(clusters),
    spearman_cluster = paste0("SPEAR_", formatC(clusters, width = 2, flag = "0"))
  )

  list(
    membership = membership,
    cor_mat = cor_mat
  )
}

module_density <- function(cor_mat, modules_df, module_col, threshold = 0.7) {
  mods <- unique(modules_df[[module_col]])
  purrr::map_dfr(mods, function(mod) {
    g <- modules_df %>% filter(.data[[module_col]] == mod) %>% pull(gene) %>% unique()
    g <- intersect(g, rownames(cor_mat))
    n <- length(g)
    if (n < 2) {
      return(tibble(
        module = mod,
        n_genes = n,
        observed_edges = 0,
        possible_edges = 0,
        density = NA_real_
      ))
    }

    csub <- cor_mat[g, g, drop = FALSE]
    upper <- csub[upper.tri(csub, diag = FALSE)]
    observed <- sum(abs(upper) > threshold, na.rm = TRUE)
    possible <- choose(n, 2)

    tibble(
      module = mod,
      n_genes = n,
      observed_edges = observed,
      possible_edges = possible,
      density = observed / possible
    )
  })
}

compute_jaccard_table <- function(spearman_df, go_df) {
  s_sets <- split(spearman_df$gene, spearman_df$spearman_cluster)
  g_sets <- split(go_df$gene, go_df$go_term)

  grid <- expand.grid(
    spearman_cluster = names(s_sets),
    go_term = names(g_sets),
    stringsAsFactors = FALSE
  )

  grid %>%
    rowwise() %>%
    mutate(
      intersection = length(intersect(s_sets[[spearman_cluster]], g_sets[[go_term]])),
      union = length(union(s_sets[[spearman_cluster]], g_sets[[go_term]])),
      jaccard = ifelse(union == 0, 0, intersection / union)
    ) %>%
    ungroup() %>%
    arrange(desc(jaccard), desc(intersection))
}

annotate_microbes <- function(species_df, epathogen_df) {
  if (!"name" %in% colnames(species_df)) {
    stop("species_list file must contain a 'name' column")
  }

  if (!"taxLineage" %in% colnames(species_df)) {
    species_df$taxLineage <- NA_character_
  }
  if (!"taxID" %in% colnames(species_df)) {
    species_df$taxID <- NA
  }

  species_df <- species_df %>%
    mutate(
      species_norm = normalize_label(name),
      taxID = suppressWarnings(as.numeric(taxID)),
      Family = extract_tax_rank(taxLineage, "f"),
      Genus = extract_tax_rank(taxLineage, "g")
    )

  epi <- epathogen_df %>%
    transmute(
      taxID = suppressWarnings(as.numeric(taxID)),
      Name = as.character(Name),
      Name_norm = normalize_label(Name),
      RiskGroup = as.character(Human.classification)
    )

  by_taxid <- epi %>%
    filter(!is.na(taxID), !is.na(RiskGroup), RiskGroup != "") %>%
    group_by(taxID) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(taxID, RiskGroup_taxid = RiskGroup)

  by_name <- epi %>%
    filter(!is.na(Name_norm), Name_norm != "", !is.na(RiskGroup), RiskGroup != "") %>%
    group_by(Name_norm) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(Name_norm, RiskGroup_name = RiskGroup)

  species_df %>%
    left_join(by_taxid, by = "taxID") %>%
    left_join(by_name, by = c("species_norm" = "Name_norm")) %>%
    mutate(
      RiskGroup = coalesce(RiskGroup_taxid, RiskGroup_name, "NotAnnotated")
    ) %>%
    select(-RiskGroup_taxid, -RiskGroup_name)
}

build_metagraph <- function(species_gene_edges, gene_to_go, microbe_annot, microbe_grouping) {
  microbe_grouping <- tolower(microbe_grouping)
  if (!microbe_grouping %in% c("risk", "family", "genus")) {
    stop("microbe_grouping must be one of: risk, family, genus")
  }

  if (microbe_grouping == "risk") {
    microbe_annot <- microbe_annot %>% mutate(microbe_meta = RiskGroup)
    meta_type <- "RiskGroup"
  } else if (microbe_grouping == "family") {
    microbe_annot <- microbe_annot %>% mutate(microbe_meta = Family)
    meta_type <- "Family"
  } else {
    microbe_annot <- microbe_annot %>% mutate(microbe_meta = Genus)
    meta_type <- "Genus"
  }

  edge_tbl <- species_gene_edges %>%
    transmute(
      species = as.character(from),
      gene = as.character(to),
      edge_weight = as.numeric(weight),
      direction = as.character(direction)
    ) %>%
    left_join(microbe_annot %>% select(name, microbe_meta), by = c("species" = "name")) %>%
    left_join(gene_to_go, by = c("gene" = "gene")) %>%
    mutate(
      microbe_meta = coalesce(microbe_meta, "Unknown_Microbe_Group"),
      go_term = coalesce(go_term, "GO_BP_Unassigned")
    )

  collapsed_edges <- edge_tbl %>%
    group_by(microbe_meta, go_term) %>%
    summarise(
      weight = mean(edge_weight, na.rm = TRUE),
      edge_count = n(),
      species_count = n_distinct(species),
      gene_count = n_distinct(gene),
      positive_fraction = mean(direction == "Positive", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(from = microbe_meta, to = go_term)

  # Added 'name' attribute explicitly to fix Cytoscape integer ID issue
  microbe_nodes <- edge_tbl %>%
    distinct(microbe_meta, species) %>%
    group_by(microbe_meta) %>%
    summarise(member_count = n_distinct(species), .groups = "drop") %>%
    transmute(
      id = microbe_meta,
      name = microbe_meta,
      node_class = "MicrobeMeta",
      grouping_method = meta_type,
      member_count = member_count
    )

  # Added 'name' attribute explicitly to fix Cytoscape integer ID issue
  gene_nodes <- edge_tbl %>%
    distinct(go_term, gene) %>%
    group_by(go_term) %>%
    summarise(member_count = n_distinct(gene), .groups = "drop") %>%
    transmute(
      id = go_term,
      name = go_term,
      node_class = "GOMeta",
      grouping_method = "GO_BP",
      member_count = member_count
    )

  nodes <- bind_rows(microbe_nodes, gene_nodes) %>% distinct(id, .keep_all = TRUE)

  graph_from_data_frame(d = collapsed_edges, vertices = nodes, directed = FALSE)
}

cfg <- parse_args()

cat("Loading inputs...\n")
species_gene_edges <- read.csv(cfg$species_gene_edges, stringsAsFactors = FALSE)
gene_mat <- read_gene_counts(cfg$filtered_gene_counts)
species_list <- read.csv(cfg$species_list, stringsAsFactors = FALSE)
epathogen_db <- read.csv(cfg$epathogen, stringsAsFactors = FALSE)
deseq_df <- read_deseq_results(cfg$deseq_dir)

required_edge_cols <- c("from", "to", "weight")
if (!all(required_edge_cols %in% colnames(species_gene_edges))) {
  stop("species_gene_edges must contain columns: from, to, weight")
}

genes_of_interest <- unique(as.character(species_gene_edges$to))
print(paste("Number of genes of interest:", length(genes_of_interest)))
print(paste("Number of genes in gene_mat:", nrow(gene_mat)))
print("Head of genes_of_interest:")
print(head(genes_of_interest, 10))
print("Head of gene_mat rownames:")
print(head(rownames(gene_mat), 10))

cat("Building fgsea gene ranking from DESeq2 outputs...\n")
ranking <- build_gene_ranking(species_gene_edges, deseq_df)
if (length(ranking$ranks) < 20) {
  stop("Too few ranked genes for stable GO enrichment. Ranked genes: ", length(ranking$ranks))
}

cat("Running GO:BP enrichment with fgsea...\n")
fgsea_tbl <- run_fgsea_go_bp(ranking$ranks)

# Format leading edges for exports
fgsea_tbl_out <- fgsea_tbl %>%
  mutate(leadingEdge = purrr::map_chr(leadingEdge, ~ paste(.x, collapse = ";")))

write.csv(fgsea_tbl_out, cfg$out_fgsea, row.names = FALSE)

cat("Exporting GO annotations for RISK...\n")
# Filter significant annotations and format for RISK load_annotation_csv()
risk_df <- fgsea_tbl_out %>%
  filter(!is.na(padj), padj < 0.1) %>%
  tidyr::unnest(leadingEdge) %>%    # expand list-column of genes
  dplyr::rename(label = pathway, nodes = leadingEdge) %>%
  dplyr::select(label, nodes)
write.csv(risk_df, cfg$out_risk_csv, row.names = FALSE)

cat("Assigning genes to GO metaclusters from leading edges...\n")
gene_to_go <- build_go_gene_modules(fgsea_tbl, genes_of_interest)

cat("Building Spearman modules for cluster comparison...\n")
spearman <- build_spearman_modules(gene_mat, genes_of_interest)

go_membership <- gene_to_go %>%
  filter(gene %in% spearman$membership$gene)
spearman_membership <- spearman$membership %>%
  filter(gene %in% go_membership$gene)

cat("Computing Jaccard overlap and edge-density differences...\n")
jaccard_tbl <- compute_jaccard_table(spearman_membership, go_membership)

spearman_density <- module_density(
  cor_mat = spearman$cor_mat,
  modules_df = spearman_membership,
  module_col = "spearman_cluster",
  threshold = 0.7
) %>%
  rename(spearman_cluster = module)

go_density <- module_density(
  cor_mat = spearman$cor_mat,
  modules_df = go_membership,
  module_col = "go_term",
  threshold = 0.7
) %>%
  rename(go_term = module)

best_match <- jaccard_tbl %>%
  group_by(spearman_cluster) %>%
  slice_head(n = 1) %>%
  ungroup()

density_comp <- best_match %>%
  left_join(spearman_density, by = "spearman_cluster") %>%
  left_join(go_density, by = "go_term", suffix = c("_spearman", "_go")) %>%
  mutate(density_diff_go_minus_spearman = density_go - density_spearman)

write.csv(jaccard_tbl, cfg$out_overlap, row.names = FALSE)
write.csv(density_comp, cfg$out_density, row.names = FALSE)

cat("Annotating microbes and creating meta-group mapping...\n")
microbe_annot <- annotate_microbes(species_list, epathogen_db)

cat("Constructing collapsed metagraph...\n")
metagraph <- build_metagraph(
  species_gene_edges = species_gene_edges,
  gene_to_go = gene_to_go,
  microbe_annot = microbe_annot,
  microbe_grouping = cfg$microbe_grouping
)

cat("Writing GraphML for Cytoscape...\n")
write_graph(metagraph, file = cfg$out_graphml, format = "graphml")

cat("Done.\n")
cat("GraphML: ", cfg$out_graphml, "\n", sep = "")
cat("fgsea table: ", cfg$out_fgsea, "\n", sep = "")
cat("RISK Annotations CSV: ", cfg$out_risk_csv, "\n", sep = "")
cat("Jaccard table: ", cfg$out_overlap, "\n", sep = "")
cat("Density comparison: ", cfg$out_density, "\n", sep = "")
