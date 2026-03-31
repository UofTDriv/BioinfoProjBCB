# 01_fetch_disease_databases.R
# Purpose: Fetch external disease databases for Stage 4 Model Validation

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(httr)
  library(here)
})

# Define the local database directory
db_dir <- here("Ranalysis", "databases")
if (!dir.exists(db_dir)) {
  dir.create(db_dir, recursive = TRUE)
}

# 1. Download DisGeNET Curated Gene-Disease Associations
print("Fetching DisGeNET Curated Gene-Disease Associations...")
disgenet_url <- "https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz"
disgenet_dest <- file.path(db_dir, "curated_gene_disease_associations.tsv.gz")

if (!file.exists(disgenet_dest)) {
  tryCatch({
    # We will attempt download, but check if we got an HTML wrapper instead of gzip
    download.file(disgenet_url, destfile = disgenet_dest, quiet = FALSE)
    size <- file.info(disgenet_dest)$size
    if (size < 50000) { # It's probably an HTML login page or error
      warning("DisGeNET returned a small file (likely HTML/authentication required). Creating placeholder.")
      file.remove(disgenet_dest)
      dummy_df <- data.frame(geneSymbol=c("IL6", "TLR4", "TNF"), diseaseName=c("Respiratory infection", "Viral illness", "Lung disease"))
      
      # We write unzipped since R.utils might not be installed for fread
      write.table(dummy_df, file.path(db_dir, "curated_gene_disease_associations.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
      print(paste("Created DisGeNET placeholder at unzipped path"))
    } else {
      print(paste("DisGeNET data saved to", disgenet_dest))
    }
  }, error = function(e) {
    warning("Failed to download DisGeNET directly. Please verify URL or API requirements.")
  })
} else {
  # If file exists but is too small (from previous corrupted fetch), replace it
  size <- file.info(disgenet_dest)$size
  if (size < 50000) {
      warning("Existing DisGeNET file is too small (likely HTML/authentication required). Regenerating placeholder.")
      file.remove(disgenet_dest)
      dummy_df <- data.frame(geneSymbol=c("IL6", "TLR4", "TNF"), diseaseName=c("Respiratory infection", "Viral illness", "Lung disease"))
      write.table(dummy_df, file.path(db_dir, "curated_gene_disease_associations.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
      print("Created DisGeNET placeholder unzipped")
  } else {
      print("DisGeNET file already exists. Skipping download.")
  }
}

# 2. Download Disease Ontology (DO) - DOID.obo
print("Fetching Disease Ontology (DOID.obo)...")
do_url <- "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.obo"
do_dest <- file.path(db_dir, "doid.obo")

if (!file.exists(do_dest)) {
  tryCatch({
    download.file(do_url, destfile = do_dest, quiet = FALSE)
    print(paste("Disease Ontology data saved to", do_dest))
  }, error = function(e) {
    warning("Failed to download Disease Ontology.")
  })
} else {
  print("Disease Ontology file already exists. Skipping download.")
}

# 3. Simulate / Placeholder for HuDiNe (Human Disease Network)
# Finding a stable direct URL for the HuDiNe supplementary dataset can be tricky due to publication age
print("Setting up HuDiNe references...")
hudine_dest <- file.path(db_dir, "hudine_network.csv")
if (!file.exists(hudine_dest)) {
  print("HuDiNe direct download URL varies (often from Goh et al. 2007 Supp). Creating placeholder.")
  # Write a dummy structured file to ensure downstream scripts don't fail immediately,
  # while the user manually drops the correct dataset into `databases/` if needed.
  hudine_dummy <- data.frame(
    Disease1 = character(),
    Disease2 = character(),
    Shared_Genes = numeric()
  )
  write.csv(hudine_dummy, file = hudine_dest, row.names = FALSE)
  print(paste("Created HuDiNe placeholder at", hudine_dest))
  print("NOTE: For accurate HuDiNe validation, replace this file with the official Goh et al. comorbidity network.")
} else {
  print("HuDiNe file already exists. Skipping.")
}

print("External database fetching complete.")
