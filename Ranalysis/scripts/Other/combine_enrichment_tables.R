# Combine enrichment tables and format for RISK import
# RISK expects: label (annotation term) and nodes (gene names separated by ;)

library(dplyr)
library(fs)
library(purrr)

# Get all enrichment CSV files from the enrichment folder
enrich_files <- dir_ls("data/other/enrichment", glob = "*.csv", recurse = FALSE)

# Read and combine all enrichment files
combined_enrichment <- enrich_files |>
  map_df(read.csv) |>
  distinct()

# Format for RISK: create label and nodes columns
risk_annotation <- combined_enrichment |>
  # Create label combining term ID and name
  mutate(
    label = paste(term.id, term.name, sep = ": "),
    # nodes column is already in the correct format (pipe-delimited genes)
    nodes = intersecting.genes
  ) |>
  # Select only the columns RISK needs
  select(label, nodes) |>
  # Remove any duplicates
  distinct()

# Write to CSV for RISK import
write.csv(
  risk_annotation,
  "data/other/enrichment/risk_annotation.csv",
  row.names = FALSE,
  quote = TRUE
)

risk_annotation