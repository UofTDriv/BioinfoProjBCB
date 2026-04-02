#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(igraph)
  library(fgsea)
  library(msigdbr)
  library(purrr)
  library(here)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

parse_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0) return(default)
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

normalize_label <- function(x) {
  x %>% tolower() %>% str_replace_all("[^a-z0-9]", "")
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
    data_dir = kv$data_dir %||% here("data", "processed", "stage3_batches"),
    annotated_species = kv$annotated_species %||% here("outputs", "annotated_species260330.csv"),
    species_corr_threshold = as.numeric(kv$species_corr_threshold %||% 0.5),
    padj_threshold = as.numeric(kv$padj_threshold %||% 0.05),
    log2fc_threshold = as.numeric(kv$log2fc_threshold %||% 1.0),
      export_dir = kv$export_dir %||% here("outputs", "stage3_balance_exports"),
    out_graphml = kv$out_graphml %||% here("outputs", "Stage3_Network_for_RISK.graphml"),
    out_homd_graphml = kv$out_homd_graphml %||% here("outputs", "Stage3_Network_HOMD_Niche.graphml"),
    out_annotation_csv = kv$out_annotation_csv %||% here("outputs", "Stage3_RISK_Annotations.csv"),
    out_validation_csv = kv$out_validation_csv %||% here("outputs", "Stage3_Jaccard_Validation.csv"),
    out_fgsea_csv = kv$out_fgsea_csv %||% here("outputs", "Stage3_fgsea_go_bp_results.csv"),
    export_phases = parse_bool(kv$export_phases, default = FALSE),
    out_species_layer = kv$out_species_layer %||% here("outputs", "Stage3_Species_only.graphml"),
    out_gene_layer = kv$out_gene_layer %||% here("outputs", "Stage3_Gene_only.graphml")
  )
}

cfg <- parse_args()

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

dir.create(cfg$export_dir, recursive = TRUE, showWarnings = FALSE)

cat("[Gather] Reconstructing Species-Species layer...\n")
norm_species <- readRDS(file.path(cfg$data_dir, "norm_species.rds"))
if (requireNamespace("WGCNA", quietly = TRUE)) {
  sp_cor <- WGCNA::cor(t(norm_species), method = "pearson", use = "pairwise.complete.obs")
} else {
  sp_cor <- cor(t(norm_species), method = "pearson", use = "pairwise.complete.obs")
}
sp_cor[is.na(sp_cor)] <- 0

# Extract upper triangle for species to avoid duplicates
sp_idx <- which(abs(sp_cor) > cfg$species_corr_threshold, arr.ind = TRUE)
sp_edges <- data.table(from=character(), to=character(), weight=numeric(), edge_class=character(), raw_value=numeric())
if(nrow(sp_idx) > 0) {
    rn <- rownames(sp_cor)[sp_idx[,1]]
    cn <- colnames(sp_cor)[sp_idx[,2]]
    valid <- rn < cn
    if(any(valid)) {
        val <- sp_cor[sp_idx][valid]
        sp_edges <- data.table(
            from = rn[valid],
            to = cn[valid],
            weight = abs(val),
            edge_class = ifelse(val > 0, "Positive", "Negative"),
            raw_value = as.numeric(val)
        )
    }
}
sp_edges[, edge_type := "Species_Species"]

cat("[Gather] Ingesting all Gene-Gene batch_*.csv concurrently via data.table...\n")
batch_files <- list.files(cfg$data_dir, pattern = "batch_[0-9]+_edges\\.csv$", full.names = TRUE)
if (length(batch_files) > 0) {
  gene_gene_edges <- rbindlist(lapply(batch_files, fread))
  gene_gene_edges[, edge_type := "Gene_Gene"]
} else {
  warning("No batch files found for genes!")
  gene_gene_edges <- data.table(from=character(), to=character(), weight=numeric(), edge_class=character(), raw_value=numeric(), edge_type=character())
}

cat("[Gather] Combining networks...\n")
species_gene_edges <- readRDS(file.path(cfg$data_dir, "species_gene_edges.rds"))
all_edges <- rbind(sp_edges, gene_gene_edges, fill=TRUE)
all_edges <- rbind(all_edges, species_gene_edges, fill=TRUE) %>% as_tibble()

correlation_edge_types <- c("Species_Species", "Gene_Gene")
all_edges <- all_edges %>%
  mutate(
    weight = ifelse(edge_type %in% correlation_edge_types, 1 - as.numeric(weight), as.numeric(weight))
  )

cat("[Gather] Loading attributes... \n")
load_annotated_species <- function(path) {
  empty_result <- tibble(
    species_name = character(),
    tax_rank = character(),
    risk_group = character(),
    homd_category = character(),
    homd = character()
  )

  if(!file.exists(path)) return(empty_result)
  ann <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("name", "RiskGroup", "HOMD.Category", "taxRank") %in% colnames(ann))) return(empty_result)
  
  ann %>%
    transmute(
      taxID = suppressWarnings(as.numeric(taxID)),
      species_name = as.character(name),
      tax_rank = as.character(taxRank),
      risk_group = as.character(RiskGroup),
      homd_category = as.character(HOMD.Category),
      homd = as.character(HOMD)
    ) %>%
    filter(!is.na(species_name), species_name != "") %>%
    distinct(species_name, .keep_all = TRUE)
}

species_meta <- load_annotated_species(cfg$annotated_species) %>%
  mutate(
    risk_group = ifelse(is.na(risk_group) | risk_group == "", "Unknown", risk_group),
    homd_category = ifelse(is.na(homd_category) | homd_category == "", "NotAnnotated", homd_category),
    tax_rank = ifelse(is.na(tax_rank) | tax_rank == "", "Unknown", tax_rank)
  )

build_nodes <- function(edges, species_set, gene_set, species_meta = NULL) {
  node_ids <- unique(c(edges$from, edges$to))
  species_nodes_from_edges <- unique(edges$from[edges$edge_type %in% c("Species_Gene", "Species_Species")])
  gene_nodes_from_edges <- unique(edges$to[edges$edge_type %in% c("Species_Gene", "Gene_Gene")])
  gene_nodes_from_gg <- unique(c(edges$from[edges$edge_type == "Gene_Gene"], edges$to[edges$edge_type == "Gene_Gene"]))

  nodes <- tibble(
    id = node_ids,
    name = node_ids,
    node_class = case_when(
      node_ids %in% species_set | node_ids %in% species_nodes_from_edges ~ "Species",
      node_ids %in% gene_set | node_ids %in% gene_nodes_from_gg ~ "Gene",
      TRUE ~ "Unknown"
    )
  )

  if (!is.null(species_meta) && nrow(species_meta) > 0) {
    nodes <- nodes %>% left_join(species_meta, by = c("name" = "species_name"))
  }
  nodes
}

build_graph <- function(edges, nodes, out_path, export_csv = FALSE, export_dir = cfg$export_dir) {
  ensure_parent_dir(out_path)
  g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  write_graph(g, out_path, format = "graphml")
  if (isTRUE(export_csv)) {
    write.csv(nodes, file.path(export_dir, "Stage3_Nodes.csv"), row.names = FALSE)
    write.csv(
      edges %>% mutate(distance = weight),
      file.path(export_dir, "Stage3_Edges.csv"),
      row.names = FALSE
    )
  }
  g
}

# Run FGSEA
cat("[Gather] Running positive-only GO:BP fgsea for RISK annotation...\n")
deseq_df <- readRDS(file.path(cfg$data_dir, "deseq_df.rds"))

run_positive_fgsea <- function(deseq_df, padj_thr, lfc_thr) {
  positive <- deseq_df %>%
    filter(!is.na(gene), !is.na(log2FoldChange), !is.na(padj)) %>%
    filter(padj < padj_thr, log2FoldChange > lfc_thr)

  if (nrow(positive) == 0) return(list(fgsea = tibble(), annotation = tibble(label = character(), nodes = character())))

  gene_score <- positive %>%
    group_by(gene) %>%
    summarise(score = max(log2FoldChange, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(score), !is.na(score))

  if (nrow(gene_score) < 10) return(list(fgsea = tibble(), annotation = tibble(label = character(), nodes = character())))

  stats <- gene_score$score
  names(stats) <- gene_score$gene
  stats <- sort(stats, decreasing = TRUE)

  msig <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>%
    select(gs_name, gene_symbol) %>%
    distinct()

  pathways <- split(msig$gene_symbol, msig$gs_name)
  pathways <- pathways[lengths(pathways) >= 10]

  if (length(pathways) == 0) return(list(fgsea = tibble(), annotation = tibble(label = character(), nodes = character())))

  fg <- fgsea(pathways = pathways, stats = stats, minSize = 10, maxSize = 500) %>%
    as_tibble() %>% filter(!is.na(padj)) %>% arrange(padj, desc(NES))

  if (nrow(fg) == 0 || !any(fg$padj < 0.05)) return(list(fgsea = fg, annotation = tibble(label = character(), nodes = character())))

  gene_to_top_term <- fg %>% filter(padj < 0.05) %>%
    select(pathway, leadingEdge) %>% unnest(cols = c(leadingEdge)) %>%
    mutate(leadingEdge = as.character(leadingEdge)) %>%
    group_by(pathway) %>% summarise(nodes = paste(leadingEdge, collapse = ";"), .groups = "drop") %>%
    transmute(label = pathway, nodes = nodes)

  list(fgsea = fg, annotation = gene_to_top_term)
}

fgsea_out <- run_positive_fgsea(deseq_df, cfg$padj_threshold, cfg$log2fc_threshold)

annotation_df <- fgsea_out$annotation %>% distinct(label, nodes)
ensure_parent_dir(cfg$out_annotation_csv)
write.csv(annotation_df, cfg$out_annotation_csv, row.names = FALSE)

ensure_parent_dir(cfg$out_fgsea_csv)
if (nrow(fgsea_out$fgsea) > 0) {
  write.csv(fgsea_out$fgsea %>% mutate(leadingEdge = map_chr(leadingEdge, ~ paste(.x, collapse = ";"))), cfg$out_fgsea_csv, row.names = FALSE)
} else {
  write.csv(tibble(), cfg$out_fgsea_csv, row.names = FALSE)
}

species_set <- unique(c(rownames(norm_species), species_gene_edges$from))
gene_set <- unique(c(readRDS(file.path(cfg$data_dir, "norm_genes.rds")) %>% rownames(), species_gene_edges$to))

nodes <- build_nodes(all_edges, species_set, gene_set, species_meta)
build_graph(all_edges, nodes, cfg$out_graphml, export_csv = TRUE)

# HOMD Niche Network
build_homd_edges <- function(species_meta) {
  annotated_species <- species_meta %>% filter(!is.na(homd_category), homd_category != "NotAnnotated", homd_category != "")
  if (nrow(annotated_species) < 2) return(tibble(from = character(), to = character(), weight = numeric(), edge_class = character(), edge_type = character(), raw_value = character()))

  edges_list <- list()
  homd_groups <- split(annotated_species$species_name, annotated_species$homd_category)
  for (homd_cat in names(homd_groups)) {
    species_in_cat <- homd_groups[[homd_cat]]
    if (length(species_in_cat) >= 2) {
      for (i in 1:(length(species_in_cat) - 1)) {
        for (j in (i + 1):length(species_in_cat)) {
          edges_list[[length(edges_list) + 1]] <- tibble(from = species_in_cat[i], to = species_in_cat[j], weight = 1.0, edge_class = "HOMD_Shared", edge_type = "HOMD_Niche", raw_value = homd_cat)
        }
      }
    }
  }
  if (length(edges_list) == 0) return(tibble())
  bind_rows(edges_list)
}

homd_edges <- build_homd_edges(species_meta)
if (nrow(homd_edges) > 0) {
  homd_nodes <- build_nodes(homd_edges, species_set, gene_set, species_meta)
  build_graph(homd_edges, homd_nodes, cfg$out_homd_graphml)
  cat("HOMD niche network built with", nrow(homd_edges), "edges\n")
}

# Validation logic
compute_coexpression_clusters <- function(gene_gene_edges, gene_set) {
  if (nrow(gene_gene_edges) == 0) return(tibble(gene = gene_set, statistical_cluster = paste0("stat_", seq_along(gene_set))))
  gg <- graph_from_data_frame(d = gene_gene_edges %>% select(from, to), directed = FALSE, vertices = tibble(name = gene_set))
  comps <- components(gg)$membership
  tibble(gene = names(comps), statistical_cluster = paste0("stat_", as.integer(comps)))
}

stat_clusters <- compute_coexpression_clusters(
  gene_gene_edges = all_edges %>% filter(edge_type=="Gene_Gene"),
  gene_set = intersect(gene_set, unique(c(all_edges$from, all_edges$to)))
)

if (nrow(annotation_df) > 0 && nrow(stat_clusters) > 0) {
  stat_sets <- split(stat_clusters$gene, stat_clusters$statistical_cluster)
  go_sets <- split(annotation_df$nodes, annotation_df$label)
  jaccard_tbl <- expand.grid(statistical_cluster = names(stat_sets), go_module = names(go_sets), stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(
      intersection = length(intersect(stat_sets[[statistical_cluster]], go_sets[[go_module]])),
      union = length(union(stat_sets[[statistical_cluster]], go_sets[[go_module]])),
      jaccard = ifelse(union == 0, 0, intersection / union)
    ) %>%
    ungroup() %>% arrange(desc(jaccard), desc(intersection))
  
  write.csv(jaccard_tbl, cfg$out_validation_csv, row.names = FALSE)
}

if (cfg$export_phases) {
  if (nrow(sp_edges) > 0) {
    build_graph(sp_edges, build_nodes(sp_edges, species_set, gene_set, species_meta), cfg$out_species_layer)
  }
  gg_edges <- all_edges %>% filter(edge_type=="Gene_Gene")
  if (nrow(gg_edges) > 0) {
    build_graph(gg_edges, build_nodes(gg_edges, species_set, gene_set, species_meta), cfg$out_gene_layer)
  }
}

cat("[Gather] Finished aggregation onto single nodeset. Stage 3 workflow complete.\n")
cat("GraphML:", cfg$out_graphml, "\n")
