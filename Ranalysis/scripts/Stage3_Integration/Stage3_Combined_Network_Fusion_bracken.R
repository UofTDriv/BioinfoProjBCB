#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(igraph)
  library(fgsea)
  library(msigdbr)
  library(foreach)
  library(doParallel)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

parse_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0) return(default)
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

normalize_label <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("[^a-z0-9]", "")
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
    filtered_gene_counts = kv$filtered_gene_counts %||%
      here("data", "processed", "filtered_gene_counts.csv"),
    species_abundance = kv$species_abundance %||%
      here("data", "processed", "brackenRanalysis", "outputs", "nonhuman_bracken260330op.csv"),
    species_list_nonhuman = kv$species_list_nonhuman %||%
      here("data", "processed", "brackenRanalysis", "outputs", "species_list_nonhuman260330op.csv"),
    epathogen = kv$epathogen %||%
      here("Ranalysis", "databases", "epathogen-2025-07-15-result.csv"),
    deseq_dir = kv$deseq_dir %||% here("outputs", "DESeq2_results_bracken"),

    padj_threshold = as.numeric(kv$padj_threshold %||% 0.05),
    log2fc_threshold = as.numeric(kv$log2fc_threshold %||% 1.0),
    species_corr_threshold = as.numeric(kv$species_corr_threshold %||% 0.5),
    gene_corr_threshold = as.numeric(kv$gene_corr_threshold %||% 0.6),

    num_cores = as.integer(kv$num_cores %||% max(1L, parallel::detectCores(logical = TRUE) - 1L)),
    export_phases = parse_bool(kv$export_phases, default = FALSE),

    out_graphml = kv$out_graphml %||%
      here("outputs", "Stage3_Network_for_RISK_bracken.graphml"),
    out_annotation_csv = kv$out_annotation_csv %||%
      here("outputs", "Stage3_RISK_Annotations_bracken.csv"),
    out_validation_csv = kv$out_validation_csv %||%
      here("outputs", "Stage3_Jaccard_Validation_bracken.csv"),
    out_fgsea_csv = kv$out_fgsea_csv %||%
      here("outputs", "Stage3_fgsea_go_bp_results_bracken.csv"),
    out_species_layer = kv$out_species_layer %||%
      here("outputs", "Stage3_Species_only_bracken.graphml"),
    out_gene_layer = kv$out_gene_layer %||%
      here("outputs", "Stage3_Gene_only_bracken.graphml")
  )
}

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

load_gene_matrix <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)

  if ("Geneid" %in% colnames(x)) x <- x %>% select(-Geneid)
  if ("Neg_Control_W" %in% colnames(x)) x <- x %>% select(-Neg_Control_W)

  x[] <- lapply(x, as.numeric)
  as.matrix(x)
}

load_species_matrix <- function(path, species_list_path = NULL) {
  sp <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  locate_sample_cols <- function(df) {
    bracken_cols <- grep("_bracken_reads$", colnames(df), value = TRUE)
    if (length(bracken_cols) > 0) return(bracken_cols)
    grep("_cladeReads$", colnames(df), value = TRUE)
  }

  sample_cols <- locate_sample_cols(sp)
  if (length(sample_cols) == 0) {
    fallbacks <- c(
      here("data", "processed", "brackenRanalysis", "outputs", "nonhuman_bracken260330op.csv"),
      here("data", "processed", "nonhuman_merged.csv"),
      here("data", "processed", "unaligned_merged.csv")
    )
    alt <- fallbacks[file.exists(fallbacks)]
    if (length(alt) > 0) {
      warning("No sample-level abundance columns found in species file. Falling back to: ", alt[[1]])
      sp <- read.csv(alt[[1]], check.names = FALSE, stringsAsFactors = FALSE)
      sample_cols <- locate_sample_cols(sp)
    }
  }

  if (!"name" %in% colnames(sp)) {
    stop("Species abundance file must contain a 'name' column.")
  }

  if (length(sample_cols) == 0) {
    stop("Species abundance file must contain sample columns ending with '_bracken_reads' or '_cladeReads'.")
  }

  if (!is.null(species_list_path) && file.exists(species_list_path)) {
    sp_list <- read.csv(species_list_path, check.names = FALSE, stringsAsFactors = FALSE)
    if ("name" %in% colnames(sp_list)) {
      sp <- sp %>% filter(name %in% sp_list$name)
    } else {
      warning("species_list_nonhuman file does not contain 'name'; skipping species filtering.")
    }
  }

  suffix_pattern <- if (all(grepl("_bracken_reads$", sample_cols))) {
    "_bracken_reads$"
  } else {
    "_cladeReads$"
  }

  mat <- sp %>%
    select(name, all_of(sample_cols)) %>%
    rename_with(~ str_replace(.x, suffix_pattern, ""), all_of(sample_cols))

  rownames(mat) <- mat$name
  mat <- mat %>% select(-name)
  mat[is.na(mat)] <- 0
  mat[] <- lapply(mat, as.numeric)

  list(
    species_info = sp,
    species_matrix = as.matrix(mat)
  )
}

load_deseq_results <- function(deseq_dir) {
  files <- list.files(deseq_dir, pattern = "_DESeq2_results\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No DESeq2 result files found in: ", deseq_dir)
  }

  out <- purrr::map_dfr(files, function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    need <- c("gene", "log2FoldChange", "padj")
    if (!all(need %in% colnames(d))) {
      return(tibble())
    }

    species_name <- basename(f) %>%
      str_remove("_DESeq2_results\\.csv$") %>%
      str_replace_all("_", " ")

    d %>%
      transmute(
        species = species_name,
        species_norm = normalize_label(species_name),
        gene = as.character(gene),
        log2FoldChange = as.numeric(log2FoldChange),
        padj = as.numeric(padj),
        pvalue = as.numeric(pvalue)
      )
  })

  if (nrow(out) == 0) {
    stop("No valid DESeq2 files found with columns gene/log2FoldChange/padj.")
  }

  out
}

load_epathogen_map <- function(path) {
  epi <- read.csv(path, stringsAsFactors = FALSE)

  if (!all(c("Name", "Human.classification") %in% colnames(epi))) {
    warning("ePathogen file missing expected columns; species risk_group will be set to Unknown.")
    return(tibble(name_norm = character(), risk_group = character(), taxID = numeric()))
  }

  epi %>%
    transmute(
      taxID = suppressWarnings(as.numeric(taxID)),
      name_norm = normalize_label(Name),
      risk_group = as.character(Human.classification)
    ) %>%
    filter(!is.na(risk_group), risk_group != "") %>%
    distinct(taxID, name_norm, .keep_all = TRUE)
}

phase1_normalize <- function(mat) {
  mat <- log2(mat + 1)
  med <- apply(mat, 2, median, na.rm = TRUE)
  sweep(mat, 2, med, FUN = "-")
}

build_species_gene_edges <- function(deseq_df, padj_thr, lfc_thr) {
  sig <- deseq_df %>%
    filter(!is.na(gene), !is.na(log2FoldChange), !is.na(padj)) %>%
    filter(padj < padj_thr, abs(log2FoldChange) > lfc_thr)

  if (nrow(sig) == 0) {
    warning("No significant Species_Gene edges found with current thresholds.")
    return(tibble(
      from = character(),
      to = character(),
      weight = numeric(),
      edge_class = character(),
      edge_type = character(),
      raw_value = numeric()
    ))
  }

  sig %>%
    transmute(
      from = species,
      to = gene,
      weight = abs(log2FoldChange),
      edge_class = ifelse(log2FoldChange > 0, "Positive", "Negative"),
      edge_type = "Species_Gene",
      raw_value = log2FoldChange
    ) %>%
    distinct()
}

make_parallel_backend <- function(num_cores) {
  cl <- parallel::makeCluster(max(1L, num_cores))
  doParallel::registerDoParallel(cl)
  cl
}

parallel_extract_edges_from_cor <- function(cor_mat, threshold, edge_type) {
  n <- nrow(cor_mat)
  if (n < 2) {
    return(tibble(
      from = character(),
      to = character(),
      weight = numeric(),
      edge_class = character(),
      edge_type = character(),
      raw_value = numeric()
    ))
  }

  idx <- seq_len(n - 1)
  edges <- foreach(i = idx, .combine = bind_rows, .packages = "dplyr") %dopar% {
    vals <- cor_mat[i, (i + 1):n]
    keep <- which(!is.na(vals) & abs(vals) > threshold)
    if (length(keep) == 0) return(tibble())

    j <- keep + i
    tibble(
      from = rownames(cor_mat)[i],
      to = rownames(cor_mat)[j],
      weight = abs(vals[keep]),
      edge_class = ifelse(vals[keep] >= 0, "Positive", "Negative"),
      edge_type = edge_type,
      raw_value = as.numeric(vals[keep])
    )
  }

  edges
}

run_positive_fgsea <- function(deseq_df, padj_thr, lfc_thr) {
  positive <- deseq_df %>%
    filter(!is.na(gene), !is.na(log2FoldChange), !is.na(padj)) %>%
    filter(padj < padj_thr, log2FoldChange > lfc_thr)

  if (nrow(positive) == 0) {
    warning("No positively upregulated genes passed thresholds for fgsea.")
    return(list(fgsea = tibble(), annotation = tibble(label = character(), nodes = character())))
  }

  gene_score <- positive %>%
    group_by(gene) %>%
    summarise(score = max(log2FoldChange, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(score), !is.na(score))

  if (nrow(gene_score) < 10) {
    warning("Too few positive genes for stable fgsea; returning empty annotation.")
    return(list(fgsea = tibble(), annotation = tibble(label = character(), nodes = character())))
  }

  stats <- gene_score$score
  names(stats) <- gene_score$gene
  stats <- sort(stats, decreasing = TRUE)

  msig <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>%
    select(gs_name, gene_symbol) %>%
    distinct()

  pathways <- split(msig$gene_symbol, msig$gs_name)
  pathways <- pathways[lengths(pathways) >= 10]

  if (length(pathways) == 0) {
    warning("No GO:BP pathways available after filtering.")
    return(list(fgsea = tibble(), annotation = tibble(label = character(), nodes = character())))
  }

  fg <- fgsea(
    pathways = pathways,
    stats = stats,
    minSize = 10,
    maxSize = 500
  ) %>%
    as_tibble() %>%
    filter(!is.na(padj)) %>%
    arrange(padj, desc(NES))

  if (nrow(fg) == 0 || !any(fg$padj < 0.05)) {
    warning("No significant GO terms remained after fgsea; writing empty annotation CSV.")
    return(list(fgsea = fg, annotation = tibble(label = character(), nodes = character())))
  }

  fg_sig <- fg %>%
    filter(padj < 0.05)

  gene_to_top_term <- fg_sig %>%
    select(pathway, leadingEdge) %>%
    unnest(cols = c(leadingEdge)) %>%
    mutate(leadingEdge = as.character(leadingEdge)) %>%
    group_by(pathway) %>%
    summarise(nodes = paste(leadingEdge, collapse = ";"), .groups = "drop") %>%
    transmute(label = pathway, nodes = nodes)

  list(fgsea = fg, annotation = gene_to_top_term)
}

build_nodes <- function(edges, species_set, gene_set, species_meta = NULL) {
  node_ids <- unique(c(edges$from, edges$to))

  nodes <- tibble(
    id = node_ids,
    name = node_ids,
    node_class = case_when(
      node_ids %in% species_set ~ "Species",
      node_ids %in% gene_set ~ "Gene",
      TRUE ~ "Unknown"
    )
  )

  if (!is.null(species_meta) && nrow(species_meta) > 0) {
    nodes <- nodes %>%
      left_join(species_meta, by = c("name" = "species_name"))
  }

  nodes
}

build_graph <- function(edges, nodes, out_path) {
  ensure_parent_dir(out_path)
  g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  write_graph(g, out_path, format = "graphml")
  g
}

compute_coexpression_clusters <- function(gene_gene_edges, gene_set) {
  if (nrow(gene_gene_edges) == 0) {
    return(tibble(gene = gene_set, statistical_cluster = paste0("stat_", seq_along(gene_set))))
  }

  gg <- graph_from_data_frame(
    d = gene_gene_edges %>% select(from, to),
    directed = FALSE,
    vertices = tibble(name = gene_set)
  )

  comps <- components(gg)$membership
  tibble(
    gene = names(comps),
    statistical_cluster = paste0("stat_", as.integer(comps))
  )
}

compute_jaccard_validation <- function(stat_clusters, annotation_df) {
  if (nrow(annotation_df) == 0 || nrow(stat_clusters) == 0) {
    return(tibble(
      statistical_cluster = character(),
      go_module = character(),
      intersection = integer(),
      union = integer(),
      jaccard = numeric()
    ))
  }

  stat_sets <- split(stat_clusters$gene, stat_clusters$statistical_cluster)
  go_sets <- split(annotation_df$nodes, annotation_df$label)

  expand.grid(
    statistical_cluster = names(stat_sets),
    go_module = names(go_sets),
    stringsAsFactors = FALSE
  ) %>%
    rowwise() %>%
    mutate(
      intersection = length(intersect(stat_sets[[statistical_cluster]], go_sets[[go_module]])),
      union = length(union(stat_sets[[statistical_cluster]], go_sets[[go_module]])),
      jaccard = ifelse(union == 0, 0, intersection / union)
    ) %>%
    ungroup() %>%
    arrange(desc(jaccard), desc(intersection))
}

# Main execution
cfg <- parse_args()

cat("[Stage3-Bracken] Loading host and microbial inputs...\n")
gene_mat_raw <- load_gene_matrix(cfg$filtered_gene_counts)
species_loaded <- load_species_matrix(cfg$species_abundance, cfg$species_list_nonhuman)
species_mat_raw <- species_loaded$species_matrix
species_info <- species_loaded$species_info
epathogen_map <- load_epathogen_map(cfg$epathogen)

deseq_df <- load_deseq_results(cfg$deseq_dir)

common_samples <- intersect(colnames(gene_mat_raw), colnames(species_mat_raw))
if (length(common_samples) < 3) {
  stop("Fewer than 3 common samples between host and microbial matrices.")
}

gene_mat <- gene_mat_raw[, common_samples, drop = FALSE]
species_mat <- species_mat_raw[, common_samples, drop = FALSE]

cat("[Stage3-Bracken] Phase 1: log2(x+1) + median-centering...\n")
gene_norm <- phase1_normalize(gene_mat)
species_norm <- phase1_normalize(species_mat)

cat("[Stage3-Bracken] Building Species_Gene edges from DESeq2 significance...\n")
species_gene_edges <- build_species_gene_edges(
  deseq_df = deseq_df,
  padj_thr = cfg$padj_threshold,
  lfc_thr = cfg$log2fc_threshold
)

if (nrow(species_gene_edges) == 0) {
  warning("No Species_Gene edges found. Graph outputs will contain only correlation layers if available.")
}

species_in_edges <- unique(species_gene_edges$from)
if (length(species_in_edges) == 0) {
  species_in_edges <- rownames(species_norm)
}

gene_in_edges <- unique(species_gene_edges$to)
if (length(gene_in_edges) == 0) {
  gene_in_edges <- rownames(gene_norm)
}

species_use <- intersect(species_in_edges, rownames(species_norm))
gene_use <- intersect(gene_in_edges, rownames(gene_norm))

if (length(species_use) < 2) {
  warning("Too few overlapping species for species-species correlations.")
  species_species_edges <- tibble(
    from = character(), to = character(), weight = numeric(),
    edge_class = character(), edge_type = character(), raw_value = numeric()
  )
} else {
  cat("[Stage3-Bracken] Phase 2: species-species Pearson correlations...\n")
  sp_cor <- suppressWarnings(cor(t(species_norm[species_use, , drop = FALSE]), method = "pearson", use = "pairwise.complete.obs"))
  sp_cor[is.na(sp_cor)] <- 0

  cl <- make_parallel_backend(cfg$num_cores)
  species_species_edges <- parallel_extract_edges_from_cor(
    cor_mat = sp_cor,
    threshold = cfg$species_corr_threshold,
    edge_type = "Species_Species"
  )
  parallel::stopCluster(cl)
}

if (length(gene_use) < 2) {
  warning("Too few overlapping genes for gene-gene correlations.")
  gene_gene_edges <- tibble(
    from = character(), to = character(), weight = numeric(),
    edge_class = character(), edge_type = character(), raw_value = numeric()
  )
} else {
  cat("[Stage3-Bracken] Building gene-gene Spearman correlations...\n")
  gg_cor <- suppressWarnings(cor(t(gene_norm[gene_use, , drop = FALSE]), method = "spearman", use = "pairwise.complete.obs"))
  gg_cor[is.na(gg_cor)] <- 0

  cl <- make_parallel_backend(cfg$num_cores)
  gene_gene_edges <- parallel_extract_edges_from_cor(
    cor_mat = gg_cor,
    threshold = cfg$gene_corr_threshold,
    edge_type = "Gene_Gene"
  )
  parallel::stopCluster(cl)
}

cat("[Stage3-Bracken] Phase 3: positive-only GO:BP fgsea for RISK annotation...\n")
fgsea_out <- run_positive_fgsea(
  deseq_df = deseq_df,
  padj_thr = cfg$padj_threshold,
  lfc_thr = cfg$log2fc_threshold
)

annotation_df <- fgsea_out$annotation %>% distinct(label, nodes)

ensure_parent_dir(cfg$out_annotation_csv)
write.csv(annotation_df, cfg$out_annotation_csv, row.names = FALSE)

ensure_parent_dir(cfg$out_fgsea_csv)
if (nrow(fgsea_out$fgsea) > 0) {
  fgsea_export <- fgsea_out$fgsea %>%
    mutate(leadingEdge = map_chr(leadingEdge, ~ paste(.x, collapse = ";")))
  write.csv(fgsea_export, cfg$out_fgsea_csv, row.names = FALSE)
} else {
  write.csv(tibble(), cfg$out_fgsea_csv, row.names = FALSE)
}

cat("[Stage3-Bracken] Phase 4: network assembly and attribute mapping...\n")
all_edges <- bind_rows(
  species_gene_edges,
  gene_gene_edges,
  species_species_edges
) %>%
  mutate(
    from = as.character(from),
    to = as.character(to),
    weight = as.numeric(weight),
    edge_class = as.character(edge_class),
    edge_type = as.character(edge_type)
  )

if (nrow(all_edges) == 0) {
  stop("No edges available to build network. Adjust thresholds or input data.")
}

species_set <- unique(c(rownames(species_norm), species_gene_edges$from))
gene_set <- unique(c(rownames(gene_norm), species_gene_edges$to))

species_meta <- species_info %>%
  transmute(
    species_name = as.character(name),
    taxID = suppressWarnings(as.numeric(taxID)),
    name_norm = normalize_label(species_name)
  ) %>%
  left_join(
    epathogen_map %>% select(taxID, name_norm, risk_group),
    by = c("taxID", "name_norm")
  ) %>%
  mutate(risk_group = ifelse(is.na(risk_group) | risk_group == "", "Unknown", risk_group)) %>%
  select(species_name, taxID, risk_group)

nodes <- build_nodes(
  all_edges,
  species_set = species_set,
  gene_set = gene_set,
  species_meta = species_meta
)

build_graph(all_edges, nodes, cfg$out_graphml)

stat_clusters <- compute_coexpression_clusters(
  gene_gene_edges = gene_gene_edges,
  gene_set = intersect(gene_use, unique(c(all_edges$from, all_edges$to)))
)

jaccard_tbl <- compute_jaccard_validation(stat_clusters, annotation_df)
ensure_parent_dir(cfg$out_validation_csv)
write.csv(jaccard_tbl, cfg$out_validation_csv, row.names = FALSE)

if (cfg$export_phases) {
  if (nrow(species_species_edges) > 0) {
    sp_nodes <- build_nodes(
      species_species_edges,
      species_set = species_set,
      gene_set = gene_set,
      species_meta = species_meta
    )
    build_graph(species_species_edges, sp_nodes, cfg$out_species_layer)
  }
  if (nrow(gene_gene_edges) > 0) {
    gg_nodes <- build_nodes(
      gene_gene_edges,
      species_set = species_set,
      gene_set = gene_set,
      species_meta = species_meta
    )
    build_graph(gene_gene_edges, gg_nodes, cfg$out_gene_layer)
  }
}

cat("[Stage3-Bracken] Done.\n")
cat("GraphML:", cfg$out_graphml, "\n")
cat("Annotation CSV:", cfg$out_annotation_csv, "\n")
cat("Validation CSV:", cfg$out_validation_csv, "\n")
