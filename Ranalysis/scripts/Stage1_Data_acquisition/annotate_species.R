#!/usr/bin/env Rscript

# Species Annotation Script for Risk Groups and HOMD Classification
# Author: Denis Rivard, modified by GitHub Copilot
# Description: Annotates species with risk groups from epathogen database and HOMD classifications

# Load required libraries
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(stringr)
})

# Utility function from utils.R
get_latest_timestamped_file <- function(input_dir = "inputs", pattern, database_type = NULL) {
  # Special handling for database files
  if (!is.null(database_type)) {
    db_dir <- here("databases")
    if (!dir.exists(db_dir)) {
      stop("Databases directory not found: ", db_dir)
    }
    
    if (database_type == "epathogen") {
      pattern <- "epathogen.*\\.csv$"
    } else if (database_type == "HOMD") {
      pattern <- "HOMD_taxon_table.*\\.(xls|xlsx|tsv|csv)$"
    }
    
    files <- list.files(path = db_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  } else {
    files <- list.files(path = here(input_dir), pattern = pattern, full.names = TRUE)
  }
  
  if (length(files) == 0) {
    if (!is.null(database_type)) {
      stop("No matching database files found for type: ", database_type)
    } else {
      stop("No matching files found in the directory.")
    }
  }
  
  # Extract various timestamp formats from filenames:
  basenames <- basename(files)
  iso_ts <- stringr::str_extract(basenames, "\\d{4}-\\d{2}-\\d{2}")
  ymd8_ts <- stringr::str_extract(basenames, "(?<!\\d)\\d{8}(?!\\d)")
  ymd6_ts <- stringr::str_extract(basenames, "(?<!\\d)\\d{6}(?!\\d)")

  dates <- rep(as.Date(NA), length(files))
  # parse in order of specificity: YYYY-MM-DD, YYYYMMDD, YYMMDD (e.g., runtime_info250819op.csv -> 250819)
  has_iso <- !is.na(iso_ts)
  if (any(has_iso)) dates[has_iso] <- as.Date(iso_ts[has_iso], format = "%Y-%m-%d")

  has_8 <- !is.na(ymd8_ts)
  if (any(has_8)) dates[has_8 & is.na(dates)] <- as.Date(ymd8_ts[has_8 & is.na(dates)], format = "%Y%m%d")

  has_6 <- !is.na(ymd6_ts)
  if (any(has_6)) dates[has_6 & is.na(dates)] <- as.Date(ymd6_ts[has_6 & is.na(dates)], format = "%y%m%d")

  if (all(is.na(dates))) {
    # Fallback to file modification time if no timestamps in filename
    file_info <- file.info(files)
    latest_file <- rownames(file_info)[which.max(file_info$mtime)]
  } else {
    # select the most recent among parsed dates
    valid_idx <- which(!is.na(dates))
    best_idx <- valid_idx[which.max(dates[valid_idx])]
    latest_file <- files[best_idx]
  }
  # Return just the file path
  latest_file
}

# Load databases function (will be called with configurable path)
load_databases <- function(databases_dir = here("Ranalysis", "databases")) {
  cat("Loading databases from:", databases_dir, "\n")
  
  # Get latest epathogen database
  epathogen_file <- get_latest_timestamped_file(database_type = "epathogen")
  cat("Using epathogen database:", basename(epathogen_file), "\n")
  epathogen_db <- read.csv(epathogen_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Get latest HOMD database
  HOMD_file <- get_latest_timestamped_file(database_type = "HOMD")
  cat("Using HOMD database:", basename(HOMD_file), "\n")
  HOMD_db <- read.table(HOMD_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t", fill = TRUE)
  
  return(list(epathogen_db = epathogen_db, HOMD_db = HOMD_db))
}

#' Annotate species with Risk Groups from epathogen database
#' 
#' @param df Input dataframe with 'name' and 'taxID' columns
#' @param epathogen_db Epathogen database dataframe
#' @return Dataframe with additional 'RiskGroup' column
annotate_risk_groups <- function(df, epathogen_db) {
  
  # Validate required columns
  required_cols <- c("name", "taxID")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  cat("Annotating risk groups using epathogen database...\n")
  cat("Epathogen database contains", nrow(epathogen_db), "entries\n")
  
  # Initialize RiskGroup column
  df$RiskGroup <- "NotAnnotated"
  
  # Step 1: Match by taxID (when both taxID and database taxID are not NA)
  taxid_matches <- df$taxID %in% epathogen_db$taxID[!is.na(epathogen_db$taxID)] & !is.na(df$taxID)
  if (any(taxid_matches)) {
    # Create lookup for taxID matches
    taxid_lookup <- epathogen_db %>%
      filter(!is.na(taxID)) %>%
      select(taxID, Human.classification) %>%
      distinct(taxID, .keep_all = TRUE)
    
    # Match by taxID
    df <- df %>%
      left_join(taxid_lookup, by = "taxID", suffix = c("", "_temp")) %>%
      mutate(RiskGroup = ifelse(!is.na(Human.classification), Human.classification, RiskGroup)) %>%
      select(-Human.classification)
    
    taxid_matched <- sum(!is.na(df$RiskGroup) & df$RiskGroup != "NotAnnotated")
    cat("Matched", taxid_matched, "species by taxID\n")
  }
  
  # Step 2: Match remaining by Name
  unmatched_by_name <- df$RiskGroup == "NotAnnotated" & !is.na(df$name)
  if (any(unmatched_by_name)) {
    # Create lookup for Name matches
    name_lookup <- epathogen_db %>%
      filter(!is.na(Name)) %>%
      select(Name, Human.classification) %>%
      distinct(Name, .keep_all = TRUE)
    
    # Match by Name (case-insensitive)
    df_unmatched <- df[unmatched_by_name, ]
    df_unmatched <- df_unmatched %>%
      mutate(name_lower = tolower(name)) %>%
      left_join(name_lookup %>% mutate(name_lower = tolower(Name)), 
                by = "name_lower", suffix = c("", "_temp")) %>%
      mutate(RiskGroup = ifelse(!is.na(Human.classification), Human.classification, RiskGroup)) %>%
      select(-Human.classification, -name_lower)
    
    # Ensure column consistency before assignment
    df_unmatched <- df_unmatched %>% select(names(df))
    df[unmatched_by_name, ] <- df_unmatched
    
    name_matched <- sum(df$RiskGroup != "NotAnnotated") - taxid_matched
    cat("Matched", name_matched, "additional species by Name\n")
  }
  
  # Summary
  total_matched <- sum(df$RiskGroup != "NotAnnotated")
  total_species <- nrow(df)
  cat("Total species annotated with risk groups:", total_matched, "out of", total_species, 
      "(", round(100 * total_matched / total_species, 1), "%)\n")
  
  return(df)
}

#' Annotate species with HOMD classifications
#' 
#' @param df Input dataframe with 'taxID' column and existing annotations
#' @param HOMD_db HOMD database dataframe
#' @return Dataframe with additional 'HOMD' and 'HOMD.Category' columns
annotate_homd <- function(df, HOMD_db) {
  
  cat("Annotating HOMD classifications...\n")
  cat("HOMD database contains", nrow(HOMD_db), "entries\n")
  
  # Initialize HOMD columns
  df$HOMD <- "NotAnnotated"
  df$HOMD.Category <- "NotAnnotated"
  
  # Match by taxID to NCBI.Taxon.ID
  homd_lookup <- HOMD_db %>%
    filter(!is.na(NCBI.Taxon.ID) & NCBI.Taxon.ID != "") %>%
    select(NCBI.Taxon.ID, Body.Site.s.) %>%
    distinct(NCBI.Taxon.ID, .keep_all = TRUE) %>%
    rename(taxID = NCBI.Taxon.ID, body_sites = Body.Site.s.)
  
  # Convert taxID to numeric if it isn't already
  if (!is.numeric(df$taxID)) {
    df$taxID <- as.numeric(df$taxID)
  }
  if (!is.numeric(homd_lookup$taxID)) {
    homd_lookup$taxID <- as.numeric(homd_lookup$taxID)
  }
  
  # Join with HOMD data
  df <- df %>%
    left_join(homd_lookup, by = "taxID") %>%
    mutate(
      HOMD = ifelse(!is.na(body_sites), body_sites, HOMD),
      # Categorize based on body sites content
      HOMD.Category = case_when(
        is.na(body_sites) | body_sites == "" ~ "NotAnnotated",
        grepl("Opportunistic Pathogen", body_sites, ignore.case = TRUE) ~ "Opportunist",
        grepl("Pathogen", body_sites, ignore.case = TRUE) & 
          !grepl("Opportunistic", body_sites, ignore.case = TRUE) ~ "Pathogen",
        !is.na(body_sites) ~ "Microbiome",
        TRUE ~ "NotAnnotated"
      )
    ) %>%
    select(-body_sites)
  
  # Summary
  homd_matched <- sum(df$HOMD != "NotAnnotated")
  total_species <- nrow(df)
  cat("Total species annotated with HOMD:", homd_matched, "out of", total_species, 
      "(", round(100 * homd_matched / total_species, 1), "%)\n")
  
  # HOMD category summary
  homd_summary <- df %>%
    count(HOMD.Category) %>%
    arrange(desc(n))
  cat("HOMD Category summary:\n")
  print(homd_summary)
  
  return(df)
}

#' Process dataframe with complete species annotation
#' 
#' @param df Input dataframe with 'name' and 'taxID' columns
#' @param databases_dir Path to databases directory
#' @return Annotated dataframe with RiskGroup, HOMD, and HOMD.Category columns
process_species_annotation <- function(df, databases_dir = here("Ranalysis", "databases")) {
  
  cat("Starting species annotation process...\n")
  cat("Input dataframe has", nrow(df), "rows\n")
  
  # Load databases
  dbs <- load_databases(databases_dir)
  
  # Annotate risk groups
  df_annotated <- annotate_risk_groups(df, dbs$epathogen_db)
  
  # Annotate HOMD
  df_annotated <- annotate_homd(df_annotated, dbs$HOMD_db)
  
  cat("Species annotation completed successfully!\n")
  
  return(df_annotated)
}

#' Parse command line arguments
#' 
#' @return List of parsed arguments
parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Default configuration
  config <- list(
    df = NULL,
    output = NULL,
    databases_dir = here("Ranalysis", "databases"),
    help = FALSE
  )
  
  # Parse command line arguments
  if (length(args) > 0) {
    i <- 1
    while (i <= length(args)) {
      arg <- args[i]
      
      if (arg == "--help" || arg == "-h") {
        config$help <- TRUE
        break
      } else if (arg == "--df") {
        if (i + 1 <= length(args)) {
          i <- i + 1
          config$df <- args[i]
        } else {
          stop("--df requires an argument")
        }
      } else if (arg == "--output") {
        if (i + 1 <= length(args)) {
          i <- i + 1
          config$output <- args[i]
        } else {
          stop("--output requires an argument")
        }
      } else if (arg == "--databases-dir") {
        if (i + 1 <= length(args)) {
          i <- i + 1
          config$databases_dir <- args[i]
        } else {
          stop("--databases-dir requires an argument")
        }
      } else {
        warning(paste("Unknown argument:", arg))
      }
      
      i <- i + 1
    }
  }
  
  return(config)
}

#' Display help message
show_help <- function() {
  cat("Usage: Rscript annotate_species.R [OPTIONS]\n\n")
  cat("Options:\n")
  cat("  --df <object>             R object name containing the dataframe (required)\n")
  cat("                            Must have 'name' and 'taxID' columns\n")
  cat("  --output <object>         R object name to store the annotated result\n")
  cat("                            (default: annotated_species)\n")
  cat("  --databases-dir <DIR>     Path to databases directory\n")
  cat("                            (default: Ranalysis/databases)\n")
  cat("  --help, -h               Show this help message\n\n")
  cat("Description:\n")
  cat("  This script annotates species with:\n")
  cat("  - RiskGroup: From epathogen database (RG1-RG4 or NotAnnotated)\n")
  cat("  - HOMD: Body site information from HOMD database\n")
  cat("  - HOMD.Category: Pathogen/Opportunist/Microbiome/NotAnnotated\n\n")
  cat("Input Requirements:\n")
  cat("  The input dataframe must contain:\n")
  cat("  - 'name': Species name column\n")
  cat("  - 'taxID': Taxonomy ID column\n\n")
  cat("Database Files:\n")
  cat("  The script automatically finds the latest timestamped database files:\n")
  cat("  - epathogen*.csv files for risk group annotation\n")
  cat("  - HOMD_taxon_table*.xls files for HOMD classification\n\n")
  cat("Examples:\n")
  cat("  run_as_script(here('scripts', 'annotate_species.R'), '--df', 'my_species_data')\n")
  cat("  run_as_script(here('scripts', 'annotate_species.R'), '--df', 'species_list', '--output', 'annotated_data')\n\n")
}

#' Main function to orchestrate species annotation
#' 
#' @param args Parsed command line arguments
main <- function(args = NULL) {
  
  # If no args provided, parse from command line
  if (is.null(args)) {
    args <- parse_arguments()
  }
  
  # Check for help request
  if (args$help) {
    show_help()
    return(invisible())
  }

  # Validate required arguments
  if (is.null(args$df)) {
    cat("Error: --df argument is required.\n\n")
    show_help()
    stop("Missing required arguments")
  }
  
  # Get dataframe from environment
  cat("Loading input dataframe from object:", args$df, "\n")
  if (!exists(args$df, envir = .GlobalEnv)) {
    stop(paste("Object", args$df, "not found in environment"))
  }
  df <- get(args$df, envir = .GlobalEnv)
  
  # Process annotation
  cat("Processing species annotation...\n")
  annotated_df <- process_species_annotation(df, args$databases_dir)
  
  # Store result in global environment
  output_name <- if (!is.null(args$output)) args$output else "annotated_species"
  assign(output_name, annotated_df, envir = .GlobalEnv)
  cat("Annotated dataframe stored as:", output_name, "\n")
  
  return(annotated_df)
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# Check for help request first
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args || "-h" %in% args) {
  show_help()
  stop("Help requested. Execution stopped.")
}

# Run main function
main()

rm(args)
