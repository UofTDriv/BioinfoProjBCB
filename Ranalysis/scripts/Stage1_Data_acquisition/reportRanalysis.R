#!/usr/bin/env Rscript

#### Header ####
# Kraken/Bracken Analysis Orchestration Script
# Author: Denis Rivard, modified by GitHub Copilot
# Last updated: August 26th 2025
# Description: Orchestrates downstream analysis of kraken2/bracken results including
#              processing kreports, species annotation with risk groups and HOMD data,
#              and correlation analysis between datasets
# Usage: Rscript reportRanalysis.R --input-dir <path> --proj-name <name> [options]

library(here)
library(parallel)
library(doParallel)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

# Utility functions from utils.R
run_as_script <- function() {
  return(!interactive() && length(sys.calls()) == 0)
}

# Utility function from output_processing.R
output_csv_file <- function(x, file_name, dir, script_name) {
  timestamp <- format(Sys.time(), "%y%m%d")
  # Construct the full file path
  file_path <- file.path(here(dir), paste0(file_name, timestamp, script_name, ".csv"))
  write.csv(x = x, file = file_path, row.names = F)
}

get_latest_timestamped_file <- function(pattern, path = ".") {
  files <- list.files(path = path, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    cat("No files found matching pattern:", pattern, "in path:", path, "\n")
    return(NULL)
  }
  
  # Extract timestamps from the files
  # Assuming format: prefix_YYMMDD_suffix
  timestamp_pattern <- "([0-9]{6})"
  timestamps <- str_extract(files, timestamp_pattern)
  
  if (all(is.na(timestamps))) {
    cat("Warning: No timestamps found in file names. Returning first file.\n")
    return(files[1])
  }
  
  # Handle files without timestamps
  valid_files <- files[!is.na(timestamps)]
  valid_timestamps <- timestamps[!is.na(timestamps)]
  
  if (length(valid_files) == 0) {
    cat("Warning: No files with valid timestamps found. Returning first file.\n")
    return(files[1])
  }
  
  # Find the most recent file
  latest_index <- which.max(valid_timestamps)
  latest_file <- valid_files[latest_index]
  
  cat("Latest file found:", basename(latest_file), "with timestamp:", valid_timestamps[latest_index], "\n")
  
  return(latest_file)
}

# Function to run scripts with command-line arguments
run_script_with_args <- function(script_path, args_vector) {
  # Handle vector of arguments
  assign("commandArgs", function(trailingOnly = TRUE) args_vector, envir = .GlobalEnv)
  # Run the script
  source(script_path)
}

# Try to load openxlsx for Excel export functionality
excel_available <- FALSE
tryCatch({
  library(openxlsx)
  excel_available <- TRUE
}, error = function(e) {
  cat("Warning: openxlsx not available, Excel export will be skipped\n")
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to display help
show_help <- function() {
  cat("Usage: Rscript reportRanalysis.R [OPTIONS]\n\n")
  cat("Description:\n")
  cat("  Orchestrates comprehensive analysis of kraken2/bracken outputs including:\n")
  cat("  - Processing kreport and bracken files from multiple datasets\n")
  cat("  - Species annotation with risk groups and HOMD classifications\n")
  cat("  - Correlation analysis between unaligned and non-human datasets\n")
  cat("  - Generation of summary reports and visualizations\n\n")
  cat("Required Options:\n")
  cat("  --input-dir <DIR>         Path to input directory containing kreport/bracken folders\n")
  cat("  --proj-name <NAME>        Project name for output directory naming\n\n")
  cat("Report Options:\n")
  cat("  --all-reports             Include all report Analyses\n")
  cat("  --rt-stats                Include Runtime Analysis\n")
  cat("  --rr-stats                Include Read/Alignment Report Statistics\n")
  cat("  --kr-stats                Include Kraken Report Statistics\n")
  cat("  --sa-stats                Include Species Annotation Statistics\n")
  cat("  --pd-report               Pathogen Detection Analysis (RG3&4 calls) \n")
  cat("  --split-report            Split each per sample report into its own pdf document\n")
  cat("  --no-per-sample           Skip analysis and report generation for each sample\n")
  cat("  --no-corr                 Skip correlation analysis and report generation\n")
  cat("  --no-params               Skip adding parameter columns\n")
  cat("  --subspecies              Include subspecies in analysis\n")
  cat("  --top-n-plot <N>          Number of top species to show in plots (default: 50)\n\n")
  cat("Optional Options:\n")
  cat("  --output-base-dir <DIR>   Base output directory (default: outputs)\n")
  cat("  --databases-dir <DIR>     Path to databases directory (default: Ranalysis/databases)\n")
  cat("  --top-n-freq <N>          Number of top species for frequency analysis (default: 25)\n")
  cat("Filtering:\n")
  cat("  --min-reads <N>           Minimum clade reads threshold (default: 0)\n")
  cat("  --exclude-taxid <ID>      Exclude species with this taxonomy ID\n")
  cat("  --minimizer-ratio <R>     Filter by minimum ratio of distinct_minimizers/cladeReads\n")
  cat("  --minimizer-threshold <N> Filter by minimum distinct_minimizers threshold\n")
  cat("  --no-protein-coding       Skip protein-coding-only gene filtering\n")
  cat("  --no-remove-mt            Skip mitochondrial gene removal\n")
  cat("  --no-go-annotation        Skip GO annotation gene filtering\n")
  cat("  --cores <N>               Number of cores to use (default: detectCores())\n")
  cat("  --help, -h                Show this help message\n\n")
  cat("Optional Manual Output Processing Inputs:\n")
  cat("  --unaligned-kreports <DIR> Path to unaligned kreports folder (default: NULL)\n")
  cat("  --unaligned-bracken <DIR>  Path to unaligned bracken folder (default: NULL)\n")
  cat("  --nonhuman-kreports <DIR>  Path to nonhuman kreports folder (default: NULL)\n")
  cat("  --nonhuman-bracken <DIR>   Path to nonhuman bracken folder (default: NULL)\n")
  cat("  --rpm-dir <DIR>           Path to RPM count matrix folder (default: data/METASEQ_pipeline/RPM)\n")
  cat("  --counts-output-dir <DIR> Path to output filtered counts (default: data/processed)\n")
  cat("  --rrna-file <FILE>        Path to rRNA gene list CSV (default: Ranalysis/databases/hgnc_rRNAgenelist.csv)\n")
  cat("  --skip-counts             Skip RPM count matrix merging and gene filtering\n")
  cat("  --runtime-dir <DIR>       Path to runtime folder for parsing runtime files (default: NULL)\n")
  cat("  --rrstats-dir <DIR>       Path to read statistics folder for parsing rrstats files (default: NULL)\n")
  cat("  --metadata-dir <DIR>      Path to directory containing sample metadata CSV files (default: NULL)\n")
  cat("Examples:\n")
  cat("  Rscript new_reportRanalysis.R --input-dir /path/to/data --proj-name my_project\n")
  cat("  Rscript new_reportRanalysis.R --input-dir /path/to/data --proj-name my_project --top-n 50 --cores 8\n")
  cat("  Rscript new_reportRanalysis.R --input-dir /path/to/data --proj-name my_project --no-corr --subspecies\n\n")
}
# Check for help request
if ("--help" %in% args || "-h" %in% args) {
  show_help()
  stop("Help requested. Execution stopped.")
} else {rm(show_help)}

# Configuration structure
create_config <- function(args = commandArgs(trailingOnly = TRUE)) {
  config <- list(
    INPUT_DIR = here("data","METASEQ_pipeline"),
    PROJ_NAME = NULL,
    RANALYSIS_DIR = here(),  # Default to current here() location
    OUTPUT_BASE_DIR = here("data","processed"),
    DATABASES_DIR = here("Ranalysis","databases"),
    TOP_N_FREQ = 25,
    TOP_N_PLOT = 50,
    MIN_READS = 10,
    EXCLUDE_TAXID = NULL,
    PERFORM_CORRELATION = TRUE,
    ADD_PARAMS = TRUE,
    INCLUDE_SUBSPECIES = FALSE,
    MINIMIZER_RATIO = NULL,
    MINIMIZER_THRESHOLD = 500,
    CORES = detectCores(),
    RPM_DIR = here("data", "METASEQ_pipeline", "RPM"),
    COUNTS_OUTPUT_DIR = here("data", "processed"),
    RRNA_FILE = here("Ranalysis", "databases", "hgnc_rRNAgenelist.csv"),
    PROCESS_COUNTS = TRUE,
    PROTEIN_CODING_FILTER = TRUE,
    REMOVE_MT_GENES = TRUE,
    GO_ANNOTATION_FILTER = TRUE,
    # Report options
    RT_STATS = FALSE,
    RR_STATS = FALSE,
    KR_STATS = FALSE,
    SA_STATS = FALSE,
    PD_STATS = FALSE,
    SPLIT_REPORT = FALSE,
    NO_PER_SAMPLE = FALSE,
    # Manual input paths
    UNALIGNED_KREPORTS_DIR = NULL,
    UNALIGNED_BRACKEN_DIR = NULL,
    NONHUMAN_KREPORTS_DIR = NULL,
    NONHUMAN_BRACKEN_DIR = NULL,
    RUNTIME_DIR = NULL,
    RRSTATS_DIR = NULL,
    METADATA_DIR = NULL
  )
  
  # Parse command line arguments
  if (length(args) > 0) {
    i <- 1
    while (i <= length(args)) {
      if (args[i] == "--input-dir" && i < length(args)) {
        config$INPUT_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--proj-name" && i < length(args)) {
        config$PROJ_NAME <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--output-base-dir" && i < length(args)) {
        config$OUTPUT_BASE_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--ranalysis-dir" && i < length(args)) {
        config$RANALYSIS_DIR <- args[i + 1]
        config$DATABASES_DIR <- file.path(args[i + 1], "databases")
        i <- i + 2
      } else if (args[i] == "--databases-dir" && i < length(args)) {
        config$DATABASES_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--top-n-freq" && i < length(args)) {
        config$TOP_N_FREQ <- as.numeric(args[i + 1])
        i <- i + 2
      } else if (args[i] == "--top-n-plot" && i < length(args)) {
        config$TOP_N_PLOT <- as.numeric(args[i + 1])
        i <- i + 2
      } else if (args[i] == "--min-reads" && i < length(args)) {
        config$MIN_READS <- as.numeric(args[i + 1])
        i <- i + 2
      } else if (args[i] == "--exclude-taxid" && i < length(args)) {
        config$EXCLUDE_TAXID <- as.numeric(args[i + 1])
        i <- i + 2
      } else if (args[i] == "--no-corr") {
        config$PERFORM_CORRELATION <- FALSE
        i <- i + 1
      } else if (args[i] == "--no-params") {
        config$ADD_PARAMS <- FALSE
        i <- i + 1
      } else if (args[i] == "--subspecies") {
        config$INCLUDE_SUBSPECIES <- TRUE
        i <- i + 1
      } else if (args[i] == "--minimizer-ratio" && i < length(args)) {
        config$MINIMIZER_RATIO <- as.numeric(args[i + 1])
        i <- i + 2
      } else if (args[i] == "--minimizer-threshold" && i < length(args)) {
        config$MINIMIZER_THRESHOLD <- as.numeric(args[i + 1])
        i <- i + 2
      } else if (args[i] == "--cores" && i < length(args)) {
        config$CORES <- as.numeric(args[i + 1])
        i <- i + 2
      # Report options
      } else if (args[i] == "--all-reports") {
        config$RT_STATS <- TRUE
        config$RR_STATS <- TRUE
        config$KR_STATS <- TRUE
        config$SA_STATS <- TRUE
        config$PD_STATS <- TRUE
        i <- i + 1
      } else if (args[i] == "--rt-stats") {
        config$RT_STATS <- TRUE
        i <- i + 1
      } else if (args[i] == "--rr-stats") {
        config$RR_STATS <- TRUE
        i <- i + 1
      } else if (args[i] == "--kr-stats") {
        config$KR_STATS <- TRUE
        i <- i + 1
      } else if (args[i] == "--sa-stats") {
        config$SA_STATS <- TRUE
        i <- i + 1
      } else if (args[i] == "--pd-stats") {
        config$PD_STATS <- TRUE
        i <- i + 1
      } else if (args[i] == "--pd-report") {
        config$PD_STATS <- TRUE
        i <- i + 1
      } else if (args[i] == "--split-report") {
        config$SPLIT_REPORT <- TRUE
        i <- i + 1
      } else if (args[i] == "--no-per-sample") {
        config$NO_PER_SAMPLE <- TRUE
        i <- i + 1
      # Manual input paths
      } else if (args[i] == "--unaligned-kreports" && i < length(args)) {
        config$UNALIGNED_KREPORTS_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--unaligned-bracken" && i < length(args)) {
        config$UNALIGNED_BRACKEN_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--nonhuman-kreports" && i < length(args)) {
        config$NONHUMAN_KREPORTS_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--nonhuman-bracken" && i < length(args)) {
        config$NONHUMAN_BRACKEN_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--runtime-dir" && i < length(args)) {
        config$RUNTIME_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--rrstats-dir" && i < length(args)) {
        config$RRSTATS_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--metadata-dir" && i < length(args)) {
        config$METADATA_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--rpm-dir" && i < length(args)) {
        config$RPM_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--counts-output-dir" && i < length(args)) {
        config$COUNTS_OUTPUT_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--rrna-file" && i < length(args)) {
        config$RRNA_FILE <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--skip-counts") {
        config$PROCESS_COUNTS <- FALSE
        i <- i + 1
      } else if (args[i] == "--no-protein-coding") {
        config$PROTEIN_CODING_FILTER <- FALSE
        i <- i + 1
      } else if (args[i] == "--no-remove-mt") {
        config$REMOVE_MT_GENES <- FALSE
        i <- i + 1
      } else if (args[i] == "--no-go-annotation") {
        config$GO_ANNOTATION_FILTER <- FALSE
        i <- i + 1
      } else if (args[i] == "--help" || args[i] == "-h") {
        show_help()
        stop("Help requested. Execution stopped.")
      } else {
        i <- i + 1
      }
    }
  }
  
  # Validate required arguments
  if (is.null(config$INPUT_DIR)) {
    stop("ERROR: --input-dir is required.")
  }
  if (is.null(config$PROJ_NAME)) {
    stop("ERROR: --proj-name is required.")
  }
  if (!dir.exists(config$INPUT_DIR)) {
    stop("ERROR: Input directory does not exist: ", config$INPUT_DIR)
  }
  if (config$SA_STATS && !dir.exists(config$DATABASES_DIR)) {
    stop("ERROR: Databases directory does not exist: ", config$DATABASES_DIR)
  }

  # Align RPM directory with input directory if available
  if (!is.null(config$INPUT_DIR)) {
    rpm_candidate <- file.path(config$INPUT_DIR, "RPM")
    if (dir.exists(rpm_candidate)) {
      config$RPM_DIR <- rpm_candidate
    }
  }
  
  # Create output paths
  config$output_dir <- file.path(config$OUTPUT_BASE_DIR, config$PROJ_NAME, "outputs")
  config$reports_dir <- file.path(config$OUTPUT_BASE_DIR, config$PROJ_NAME, "reports")
  
  return(config)
}

# Setup parallel processing
setup_parallel <- function(cores) {
  cat("Setting up parallel processing with", cores, "cores\n")
  # Reuse existing cluster if present and valid
  if (exists("cl", envir = .GlobalEnv)) {
    cl_existing <- tryCatch(get("cl", envir = .GlobalEnv), error = function(e) NULL)
    if (!is.null(cl_existing) && inherits(cl_existing, "cluster")) {
      workers <- tryCatch(getDoParWorkers(), error = function(e) NA)
      if (!is.na(workers) && workers > 0) {
        cat("Existing cluster detected with", workers, "workers. Reusing it.\n")
        registerDoParallel(cl_existing)
        return(cl_existing)
      } else {
        cat("Existing cluster object found but not registered. Registering it now.\n")
        registerDoParallel(cl_existing)
        return(cl_existing)
      }
    }
  }
  # No valid existing cluster found -> create a new one
  start_time <- Sys.time()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  assign("cl", cl, envir = .GlobalEnv)
  end_time <- Sys.time()
  workers <- tryCatch(getDoParWorkers(), error = function(e) NA)
  cat("Parallel cluster initialized with", workers, "workers, in", 
      round(difftime(end_time, start_time, units = "secs")), "seconds\n")
  return(cl)
}

# Build arguments for output_processing.R (updated to work with new version)
build_output_processing_args <- function(config) {
  # Use manual paths if provided, otherwise auto-detect
  unaligned_kreports_dir <- config$UNALIGNED_KREPORTS_DIR
  nonhuman_kreports_dir <- config$NONHUMAN_KREPORTS_DIR
  unaligned_bracken_dir <- config$UNALIGNED_BRACKEN_DIR
  nonhuman_bracken_dir <- config$NONHUMAN_BRACKEN_DIR
  
  # Auto-detect directories only if manual paths are not provided
  if (is.null(unaligned_kreports_dir) || is.null(nonhuman_kreports_dir)) {
    # Try different expected directory structures
    if (dir.exists(file.path(config$INPUT_DIR, "kraken2"))) {
      # Standard pipeline structure - look for subdirectories
      kraken2_dirs <- list.dirs(file.path(config$INPUT_DIR, "kraken2"), full.names = TRUE, recursive = FALSE)
      unaligned_dirs <- grep("unaligned", kraken2_dirs, value = TRUE)
      nonhuman_dirs <- grep("nonhuman", kraken2_dirs, value = TRUE)
      
      if (is.null(unaligned_kreports_dir) && length(unaligned_dirs) > 0) unaligned_kreports_dir <- unaligned_dirs[1]
      if (is.null(nonhuman_kreports_dir) && length(nonhuman_dirs) > 0) nonhuman_kreports_dir <- nonhuman_dirs[1]
    } else if (dir.exists(file.path(config$INPUT_DIR, "unaligned_kreports"))) {
      if (is.null(unaligned_kreports_dir)) unaligned_kreports_dir <- file.path(config$INPUT_DIR, "unaligned_kreports")
    } else if (dir.exists(file.path(config$INPUT_DIR, "kreports"))) {
      # Use general kreports directory for both
      if (is.null(unaligned_kreports_dir)) unaligned_kreports_dir <- file.path(config$INPUT_DIR, "kreports")
    } else {
      # Check if INPUT_DIR itself contains kreport files (direct directory case)
      kreport_files <- list.files(config$INPUT_DIR, pattern = "\\.(kreport|k2report)$", recursive = TRUE)
      if (length(kreport_files) > 0) {
        # Check subdirectories for unaligned and nonhuman patterns
        subdirs <- list.dirs(config$INPUT_DIR, full.names = TRUE, recursive = FALSE)
        unaligned_subdirs <- grep("unaligned", subdirs, value = TRUE)
        nonhuman_subdirs <- grep("nonhuman", subdirs, value = TRUE)
        
        if (is.null(unaligned_kreports_dir) && length(unaligned_subdirs) > 0) {
          unaligned_kreports_dir <- unaligned_subdirs[1]
        }
        if (is.null(nonhuman_kreports_dir) && length(nonhuman_subdirs) > 0) {
          nonhuman_kreports_dir <- nonhuman_subdirs[1]
        }
      }
    }
    
    if (is.null(nonhuman_kreports_dir)) {
      if (dir.exists(file.path(config$INPUT_DIR, "kraken2"))) {
        # Look for nonhuman subdirectories within the main kraken2 folder
        kraken2_dirs <- list.dirs(file.path(config$INPUT_DIR, "kraken2"), full.names = TRUE, recursive = FALSE)
        nonhuman_dirs <- grep("nonhuman", kraken2_dirs, value = TRUE)
        if (length(nonhuman_dirs) > 0) nonhuman_kreports_dir <- nonhuman_dirs[1]
      } else if (dir.exists(file.path(config$INPUT_DIR, "kraken2_nonhuman"))) {
        # Fallback: check for legacy structure (kraken2_nonhuman folder) - deprecated
        cat("Warning: Using legacy kraken2_nonhuman directory structure. Consider updating to new kraken2/{param_label}_nonhuman structure.\n")
        kraken2_nonhuman_dirs <- list.dirs(file.path(config$INPUT_DIR, "kraken2_nonhuman"), full.names = TRUE, recursive = FALSE)
        nonhuman_dirs <- grep("nonhuman", kraken2_nonhuman_dirs, value = TRUE)
        if (length(nonhuman_dirs) > 0) nonhuman_kreports_dir <- nonhuman_dirs[1]
      } else if (dir.exists(file.path(config$INPUT_DIR, "nonhuman_kreports"))) {
        nonhuman_kreports_dir <- file.path(config$INPUT_DIR, "nonhuman_kreports")
      } else if (dir.exists(file.path(config$INPUT_DIR, "kreports"))) {
        nonhuman_kreports_dir <- file.path(config$INPUT_DIR, "kreports")
      }
    }
  }
  
  # Check if we found valid kreport directories
  if (is.null(unaligned_kreports_dir)) {
    stop("ERROR: Could not find unaligned kreports directory in: ", config$INPUT_DIR)
  }
  
  # Count samples in unaligned directory
  unaligned_samples <- 0
  if (dir.exists(unaligned_kreports_dir)) {
    unaligned_files <- list.files(unaligned_kreports_dir, pattern = "\\.(kreport|k2report)$")
    unaligned_samples <- length(unaligned_files)
  }
  
  # Check nonhuman directory and validate sample count match
  nonhuman_valid <- FALSE
  if (!is.null(nonhuman_kreports_dir) && dir.exists(nonhuman_kreports_dir)) {
    nonhuman_files <- list.files(nonhuman_kreports_dir, pattern = "\\.(kreport|k2report)$")
    nonhuman_samples <- length(nonhuman_files)
    
    if (nonhuman_samples == unaligned_samples && nonhuman_samples > 0) {
      nonhuman_valid <- TRUE
      cat("Found matching nonhuman kreports directory with", nonhuman_samples, "samples\n")
    } else {
      cat("Warning: Nonhuman kreports directory found but sample count mismatch (unaligned:", unaligned_samples, ", nonhuman:", nonhuman_samples, ")\n")
      nonhuman_kreports_dir <- NULL
    }
  } else {
    cat("Warning: No valid nonhuman kreports directory found. Analysis will use only unaligned data.\n")
    nonhuman_kreports_dir <- NULL
  }
  
  args <- c(
    "--unaligned-kreports", unaligned_kreports_dir,
    "--output-dir", config$output_dir,
    "--top-n", as.character(config$TOP_N_FREQ),
    "--min-reads", as.character(config$MIN_READS),
    "--output"
  )

  if (!is.null(config$RPM_DIR)) {
    args <- c(args, "--rpm-dir", config$RPM_DIR)
  }
  if (!is.null(config$COUNTS_OUTPUT_DIR)) {
    args <- c(args, "--counts-output-dir", config$COUNTS_OUTPUT_DIR)
  }
  if (!is.null(config$RRNA_FILE)) {
    args <- c(args, "--rrna-file", config$RRNA_FILE)
  }
  if (!config$PROCESS_COUNTS) {
    args <- c(args, "--skip-counts")
  }
  if (!config$PROTEIN_CODING_FILTER) {
    args <- c(args, "--no-protein-coding")
  }
  if (!config$REMOVE_MT_GENES) {
    args <- c(args, "--no-remove-mt")
  }
  if (!config$GO_ANNOTATION_FILTER) {
    args <- c(args, "--no-go-annotation")
  }
  
  if (config$CORES > 1) {
    args <- c(args, "--parallel")
  } else {
    cat("Note: Running in single-core mode. For faster processing, consider using multiple cores with --cores <N> > 1.\n")
  }
  
  # Only add nonhuman kreports if valid
  if (nonhuman_valid && !is.null(nonhuman_kreports_dir)) {
    args <- c(args, "--nonhuman-kreports", nonhuman_kreports_dir)
  }
  
  # Add bracken directories (use manual paths if provided, otherwise auto-detect)
  if (!is.null(unaligned_bracken_dir)) {
    args <- c(args, "--unaligned-bracken", unaligned_bracken_dir)
  } else if (dir.exists(file.path(config$INPUT_DIR, "bracken"))) {
    # Standard pipeline structure - look for subdirectories
    bracken_dirs <- list.dirs(file.path(config$INPUT_DIR, "bracken"), full.names = TRUE, recursive = FALSE)
    unaligned_bracken_dirs <- grep("unaligned", bracken_dirs, value = TRUE)
    
    if (length(unaligned_bracken_dirs) > 0) {
      args <- c(args, "--unaligned-bracken", unaligned_bracken_dirs[1])
    }
  } else {
    # Try alternative directory structures
    if (dir.exists(file.path(config$INPUT_DIR, "unaligned_bracken"))) {
      args <- c(args, "--unaligned-bracken", file.path(config$INPUT_DIR, "unaligned_bracken"))
    }
  }
  
  # Only add nonhuman bracken if nonhuman kreports are valid
  if (nonhuman_valid) {
    if (!is.null(nonhuman_bracken_dir)) {
      args <- c(args, "--nonhuman-bracken", nonhuman_bracken_dir)
    } else if (dir.exists(file.path(config$INPUT_DIR, "bracken"))) {
      # Standard pipeline structure - look for subdirectories
      bracken_dirs <- list.dirs(file.path(config$INPUT_DIR, "bracken"), full.names = TRUE, recursive = FALSE)
      nonhuman_bracken_dirs <- grep("nonhuman", bracken_dirs, value = TRUE)
      
      if (length(nonhuman_bracken_dirs) > 0) {
        args <- c(args, "--nonhuman-bracken", nonhuman_bracken_dirs[1])
      }
    } else {
      # Try alternative directory structures
      if (dir.exists(file.path(config$INPUT_DIR, "nonhuman_bracken"))) {
        args <- c(args, "--nonhuman-bracken", file.path(config$INPUT_DIR, "nonhuman_bracken"))
      }
      
      # Also try legacy bracken directory structure 
      if (dir.exists(file.path(config$INPUT_DIR, "bracken"))) {
        bracken_dir <- file.path(config$INPUT_DIR, "bracken")
        # Note: param_label detection would need to be implemented for this fallback
        unaligned_bracken_files <- list.files(bracken_dir, pattern = "_unaligned_.*\\.bracken$", full.names = TRUE)
        nonhuman_bracken_files <- list.files(bracken_dir, pattern = "_nonhuman_.*\\.bracken$", full.names = TRUE)
        
        if (length(unaligned_bracken_files) > 0) {
          args <- c(args, "--unaligned-bracken", bracken_dir)
        }
        if (length(nonhuman_bracken_files) > 0) {
          args <- c(args, "--nonhuman-bracken", bracken_dir)
        }
      }
    }
  }
  
  # Add metadata arguments (use manual paths if provided, otherwise auto-detect)
  if (!is.null(config$RUNTIME_DIR)) {
    args <- c(args, "--runtime-dir", config$RUNTIME_DIR)
  } else if (dir.exists(config$INPUT_DIR) && file.exists(file.path(config$INPUT_DIR, "runtime"))) {
    args <- c(args, "--runtime-dir", config$INPUT_DIR)
  }
  
  if (!is.null(config$RRSTATS_DIR)) {
    args <- c(args, "--rrstats-dir", config$RRSTATS_DIR)
  } else if (dir.exists(config$INPUT_DIR) && dir.exists(file.path(config$INPUT_DIR, "QC"))) {
    args <- c(args, "--rrstats-dir", file.path(config$INPUT_DIR, "QC"))
  }
  
  if (!is.null(config$METADATA_DIR)) {
    args <- c(args, "--metadata-dir", config$METADATA_DIR)
  } else if (dir.exists(config$INPUT_DIR)) {
    metadata_files <- list.files(config$INPUT_DIR, pattern = "sample_metadata.*\\.csv$", full.names = TRUE, recursive = FALSE)
    if (length(metadata_files) > 0) {
      args <- c(args, "--metadata-dir", config$INPUT_DIR)
    }
  }
  
  # Add optional flags
  if (config$ADD_PARAMS) {
    args <- c(args, "--params")
  }
  if (config$INCLUDE_SUBSPECIES) {
    args <- c(args, "--subspecies")
  }
  if (!is.null(config$EXCLUDE_TAXID)) {
    args <- c(args, "--exclude-taxid", as.character(config$EXCLUDE_TAXID))
  }
  if (!is.null(config$MINIMIZER_RATIO)) {
    args <- c(args, "--minimizer-ratio", as.character(config$MINIMIZER_RATIO))
  }
  if (!is.null(config$MINIMIZER_THRESHOLD)) {
    args <- c(args, "--minimizer-threshold", as.character(config$MINIMIZER_THRESHOLD))
  }
  
  return(args)
}

# ============================================================================
# Analysis Reports
# ============================================================================

# Function to create correlation analysis
create_correlation_analysis <- function(data1, data2, data1_name, data2_name, data1_col_suffix, data2_col_suffix, output_dir, top_n_label = config$TOP_N_PLOT) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  top_n_label <- max(1, floor(config$TOP_N_PLOT / 2))
  
  # Get sample columns
  cols1 <- grep(paste0("_", data1_col_suffix, "$"), colnames(data1), value = TRUE)
  cols2 <- grep(paste0("_", data2_col_suffix, "$"), colnames(data2), value = TRUE)
  
  samples1 <- str_remove(cols1, paste0("_", data1_col_suffix, "$"))
  samples2 <- str_remove(cols2, paste0("_", data2_col_suffix, "$"))
  
  common_samples <- intersect(samples1, samples2)
  
  if (length(common_samples) == 0) {
    warning("No common samples found between datasets")
    return(NULL)
  }
  
  cat("Found", length(common_samples), "common samples for correlation analysis\n")
  
  start_time <- Sys.time()
  
  # Process samples with parallel processing when available
  if (exists("cl", envir = .GlobalEnv) && !is.null(get("cl", envir = .GlobalEnv))) {
    cat("Processing correlation analysis in parallel using", length(get("cl", envir = .GlobalEnv)), "cores...\n")
    
    # Set up parallel backend
    cluster <- get("cl", envir = .GlobalEnv)
    
    # Export required variables to cluster nodes
    clusterExport(cluster, c("data1", "data2", "data1_col_suffix", "data2_col_suffix", 
                             "data1_name", "data2_name"), envir = environment())
    
    # Parallel processing using foreach
    sample_results <- foreach(sample_id = common_samples, 
                              .packages = c("dplyr", "ggplot2", "ggrepel"),
                              .combine = 'c') %dopar% {
      col1 <- paste0(sample_id, "_", data1_col_suffix)
      col2 <- paste0(sample_id, "_", data2_col_suffix)
      
      if (col1 %in% colnames(data1) && col2 %in% colnames(data2)) {
        # Select data from both datasets
        data1_subset <- data1 %>% select(taxID, name, all_of(col1))
        data2_subset <- data2 %>% select(taxID, name, all_of(col2))
        
        # Rename the columns to generic names before joining
        colnames(data1_subset)[3] <- "reads1"
        colnames(data2_subset)[3] <- "reads2"
        
        sample_comparison <- inner_join(
          data1_subset,
          data2_subset,
          by = c("taxID", "name")
        ) %>%
          filter(reads1 > 0 & reads2 > 0)
        
        cor_coef <- cor(log10(sample_comparison$reads1 + 1), 
                        log10(sample_comparison$reads2 + 1), 
                        use = "complete.obs")
        
        # Create labels for only the top_n_label most different points
        sample_comparison <- sample_comparison %>%
          mutate(
            log_reads1 = log10(reads1 + 1),
            log_reads2 = log10(reads2 + 1),
            diff = abs(log_reads1 - log_reads2)
          ) %>%
          arrange(desc(diff)) %>%
          mutate(
            rank = row_number(),
            label = ifelse(rank <= top_n_label, name, NA_character_)
          ) %>%
          select(-diff, -rank)
        
        # Create the correlation plot
        gg_correlation <- ggplot(sample_comparison, aes(x = log_reads1, y = log_reads2)) +
          geom_point(alpha = 0.6, size = 1.5) + 
          geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted") +
          geom_text_repel(aes(label = label), na.rm = TRUE, size = 2.5, max.overlaps = top_n_label) +
          labs(
            x = paste0("Log10(", data1_name, " ", data1_col_suffix, " + 1)"),
            y = paste0("Log10(", data2_name, " ", data2_col_suffix, " + 1)"),
            title = paste0("Correlation: ", data1_name, " vs ", data2_name, " (", sample_id, ")"),
            subtitle = paste0("Pearson r = ", round(cor_coef, 5), " (n = ", nrow(sample_comparison), " species)")
          ) + coord_fixed()
        
        list(list(
          sample = sample_id,
          pearson_r = cor_coef,
          n_species = nrow(sample_comparison),
          plot = gg_correlation
        ))
      } else {
        list(NULL)
      }
    }
    
    # Remove NULL entries
    sample_results <- sample_results[!sapply(sample_results, is.null)]
    
  } else {
    # Fallback to sequential processing
    cat("Processing correlation analysis sequentially...\n")
    sample_results <- list()
    for (sample_id in common_samples) {
      col1 <- paste0(sample_id, "_", data1_col_suffix)
      col2 <- paste0(sample_id, "_", data2_col_suffix)
      
      if (col1 %in% colnames(data1) && col2 %in% colnames(data2)) {
        # Select data from both datasets
        data1_subset <- data1 %>% select(taxID, name, all_of(col1))
        data2_subset <- data2 %>% select(taxID, name, all_of(col2))
        
        # Rename the columns to generic names before joining
        colnames(data1_subset)[3] <- "reads1"
        colnames(data2_subset)[3] <- "reads2"
        
        sample_comparison <- inner_join(
          data1_subset,
          data2_subset,
          by = c("taxID", "name")
        ) %>%
          filter(reads1 > 0 & reads2 > 0)
        
        cor_coef <- cor(log10(sample_comparison$reads1 + 1), 
                        log10(sample_comparison$reads2 + 1), 
                        use = "complete.obs")
        
        # Create labels for only the top_n_label most different points
        sample_comparison <- sample_comparison %>%
          mutate(
            log_reads1 = log10(reads1 + 1),
            log_reads2 = log10(reads2 + 1),
            diff = abs(log_reads1 - log_reads2)
          ) %>%
          arrange(desc(diff)) %>%
          mutate(
            rank = row_number(),
            label = ifelse(rank <= top_n_label, name, NA_character_)
          ) %>%
          select(-diff, -rank)
        
        # Create the correlation plot
        gg_correlation <- ggplot(sample_comparison, aes(x = log_reads1, y = log_reads2)) +
            geom_point(alpha = 0.6, size = 1.5) + 
            geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted") +
            geom_text_repel(aes(label = label), na.rm = TRUE, size = 2.5, max.overlaps = top_n_label) +
            labs(
              x = paste0("Log10(", data1_name, " ", data1_col_suffix, " + 1)"),
              y = paste0("Log10(", data2_name, " ", data2_col_suffix, " + 1)"),
              title = paste0("Correlation: ", data1_name, " vs ", data2_name, " (", sample_id, ")"),
              subtitle = paste0("Pearson r = ", round(cor_coef, 5), " (n = ", nrow(sample_comparison), " species)")
            ) + coord_fixed() 
        
        sample_results[[length(sample_results) + 1]] <- list(
          sample = sample_id,
          pearson_r = cor_coef,
          n_species = nrow(sample_comparison),
          plot = gg_correlation
        )
      }
    }
  }
  
  end_time <- Sys.time()
  cat("Finished correlation analysis in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
  
  # Extract correlation results and plots
  correlation_results <- tibble(
    sample = character(),
    pearson_r = numeric(),
    n_species = integer()
  )
  
  plots_list <- list()
  last_plot <- NULL
  
  # Process results (compatible with both parallel and sequential)
  for (result in sample_results) {
    if (!is.null(result)) {
      correlation_results <- bind_rows(
        correlation_results,
        tibble(
          sample = result$sample,
          pearson_r = result$pearson_r,
          n_species = result$n_species
        )
      )
      plots_list[[result$sample]] <- result$plot
      last_plot <- result$plot
    }
  }
  
  cat("Saving plots to PDF\n")
  # Save plots to PDF
  pdf(file.path(output_dir, paste0("correlation_", data1_name, "_vs_", data2_name, "_by_sample.pdf")), 
      width = 12, height = 8)
  for (plot in plots_list) {
    print(plot)
  }
  dev.off()
  cat("Correlation plots saved to PDF\n")
  
  # Print correlation results table
  cat("\n=== Correlation Analysis Results ===\n")
  correlation_results <- correlation_results %>%
    mutate(pearson_r = round(pearson_r, 5))
  print(correlation_results)
  cat("=====================================\n")
  
  return(list(common_samples = common_samples, correlation_results = correlation_results))
}

runread_stats_report <- function(runread_data, report_dir) {
  pdf(file.path(report_dir, "run_read_stats_report.pdf"), width = 8, height = 6)
  
  title_text <- "Part 2: Batch RNAseq Statistics"
  
  # Create reference text instead of displaying the table
  ref_text <- paste0(
    "RNAseq statistics have been compiled and saved as CSV files.\n\n",
    "The RNAseq statistics data includes:\n",
    "- Sample names\n",
    "- Number of input reads\n",
    "- Average input read length\n",
    "- Uniquely mapped reads\n",
    "- Average mapped length\n", 
    "- Reads assigned to genes\n",
    "- Percent mapped reads\n",
    "- Uniquely mapped percent\n",
    "- Flag for samples with <10M unique reads\n\n",
    "Data source: Parsed from run_read_stats*.txt files in rrstats/ folder\n",
    "Output location: Check outputs/ directory for CSV files containing this data\n\n",
    "Key findings:\n",
    "- Total samples with read statistics: ", nrow(runread_data), "\n",
    "- Samples with <10M unique reads: ", sum(runread_data$less_than_10M_unique, na.rm = TRUE), "\n",
    "- Average input reads: ", format(round(mean(runread_data$num_input_reads, na.rm = TRUE)), big.mark = ","), "\n",
    "- Average mapping percentage: ", round(mean(runread_data$percent_mapped_reads, na.rm = TRUE), 1), "%\n",
    "- Average unique mapping percentage: ", round(mean(runread_data$uniquely_mapped_percent, na.rm = TRUE), 1), "%"
  )
  
  summary_plot <- ggplot() + 
    theme_void() + 
    labs(title = title_text, subtitle = ref_text)
  print(summary_plot)
  
  dev.off()
  
  # Export to Excel with conditional formatting if openxlsx is available
  tryCatch({
    if (excel_available) {
      excel_file <- file.path(report_dir, "run_read_stats_report.xlsx")
      
      # Create workbook
      wb <- createWorkbook()
      addWorksheet(wb, "Read_Stats")

      # write.csv(runread_data, file = file.path(report_dir, "run_read_stats_report.csv"), row.names = FALSE)
      # runread_data <- read.table(file.path(report_dir, "run_read_stats_report.csv"), header = TRUE, sep = ",")
      
      # Create headers without underscores
      headers <- colnames(runread_data)
      clean_headers <- gsub("_", " ", headers)
      clean_headers <- tools::toTitleCase(clean_headers)
      
      # Prepare data for Excel export - ensure proper data types
      excel_data <- runread_data
      
      # Convert numeric columns to ensure they're properly formatted
      numeric_cols <- c("num_input_reads", "avg_input_read_length", "uniquely_mapped_reads", 
                       "avg_mapped_length", "reads_assigned_to_genes", "percent_mapped_reads", 
                       "uniquely_mapped_percent")
      
      for (col in numeric_cols) {
        if (col %in% colnames(excel_data)) {
          excel_data[[col]] <- as.numeric(excel_data[[col]])
        }
      }
      
      # Ensure logical column is properly formatted
      if ("less_than_10M_unique" %in% colnames(excel_data)) {
        excel_data[["less_than_10M_unique"]] <- as.logical(excel_data[["less_than_10M_unique"]])
      }
      
      # Ensure sample column is character
      if ("sample" %in% colnames(excel_data)) {
        excel_data[["sample"]] <- as.character(excel_data[["sample"]])
      }
      
      # Write headers and data separately
      writeData(wb, "Read_Stats", t(clean_headers), startRow = 1, colNames = FALSE)
      writeData(wb, "Read_Stats", excel_data, startRow = 2, colNames = FALSE)
      
      # Add conditional formatting
      data_rows <- 2:(nrow(excel_data) + 1)
      # Data bars for num_input_reads
      if ("num_input_reads" %in% headers) {
        num_reads_col <- which(headers == "num_input_reads")
        conditionalFormatting(wb, "Read_Stats",
                            cols = num_reads_col,
                            rows = data_rows,
                            type = "dataBar",
                            style = "#0070C0")
        cat("Applied data bars to num_input_reads column\n")
      }
      
      # THESE BREAK THE EXCEL FILE - DISABLED FOR NOW
      # # Green data bars for uniquely_mapped_reads
      # if ("uniquely_mapped_reads" %in% headers) {
      #   unique_col <- which(headers == "uniquely_mapped_reads")
      #   conditionalFormatting(wb, "Read_Stats", 
      #                       cols = unique_col, 
      #                       rows = data_rows,
      #                       type = "dataBar",
      #                       style = "#00B050")
      #   cat("Applied green data bars to uniquely_mapped_reads column\n")
      # }
      # # Orange data bars for reads_assigned_to_genes
      # if ("reads_assigned_to_genes" %in% headers) {
      #   assigned_col <- which(headers == "reads_assigned_to_genes")
      #   conditionalFormatting(wb, "Read_Stats",
      #                       cols = assigned_col,
      #                       rows = data_rows,
      #                       type = "dataBar",
      #                       style = "#FF8C00")
      #   cat("Applied orange data bars to reads_assigned_to_genes column\n")
      # }
    
      # Red-White-Blue color scale for avg_mapped_length
      if ("avg_mapped_length" %in% headers) {
        avg_length_col <- which(headers == "avg_mapped_length")
        conditionalFormatting(wb, "Read_Stats",
                            cols = avg_length_col,
                            rows = data_rows,
                            type = "colorScale",
                            style = c("#FF0000", "#FFFFFF", "#0000FF"))
        # Add borders to color scale cells
        border_style <- createStyle(border = "TopBottomLeftRight", borderStyle = "thin")
        addStyle(wb, "Read_Stats", border_style, rows = data_rows, cols = avg_length_col)
        cat("Applied red-white-blue color scale to avg_mapped_length column\n")
      }
      
      # Green-Yellow color scale for percent_mapped_reads
      if ("percent_mapped_reads" %in% headers) {
        percent_mapped_col <- which(headers == "percent_mapped_reads")
        conditionalFormatting(wb, "Read_Stats",
                            cols = percent_mapped_col,
                            rows = data_rows,
                            type = "colorScale",
                            style = c("#00FF00", "#FFFF00"))
        # Add borders to color scale cells
        border_style <- createStyle(border = "TopBottomLeftRight", borderStyle = "thin")
        addStyle(wb, "Read_Stats", border_style, rows = data_rows, cols = percent_mapped_col)
        cat("Applied green-yellow color scale to percent_mapped_reads column\n")
      }

      # Green-Yellow color scale for uniquely_mapped_percent
      if ("uniquely_mapped_percent" %in% headers) {
        unique_percent_col <- which(headers == "uniquely_mapped_percent")
        conditionalFormatting(wb, "Read_Stats",
                            cols = unique_percent_col,
                            rows = data_rows,
                            type = "colorScale",
                            style = c("#00FF00", "#FFFF00"))
        # Add borders to color scale cells
        border_style <- createStyle(border = "TopBottomLeftRight", borderStyle = "thin")
        addStyle(wb, "Read_Stats", border_style, rows = data_rows, cols = unique_percent_col)
        cat("Applied green-yellow color scale to uniquely_mapped_percent column\n")
      }
      
      # Auto-size columns
      setColWidths(wb, "Read_Stats", cols = 1:length(headers), widths = "auto")
      
      # Add header formatting
      header_style <- createStyle(
        textDecoration = "bold",
        fgFill = "#E6E6FA",
        border = "TopBottomLeftRight",
        borderStyle = "thin"
      )
      addStyle(wb, "Read_Stats", header_style, rows = 1, cols = 1:length(headers))
      
      # Save workbook
      saveWorkbook(wb, excel_file, overwrite = TRUE)
      cat("Read statistics exported to Excel:", excel_file, "\n")
    } else {
      cat("Excel export skipped - openxlsx not available\n")
    }
  }, error = function(e) {
    cat("Warning: Failed to export read statistics to Excel:", e$message, "\n")
    print(e)
  })
}

kraken_stats_report <- function(unaligned_results, nonhuman_results = NULL, combined_report_data, config, report_dir, batch) {
  if (batch) {
    pdf(file.path(report_dir, "kraken_stats_report.pdf"), width = 8, height = 6)
    
    title_text <- "Part 3: Kraken2 Statistics"
    
    # Check if nonhuman data is available
    has_nonhuman <- !is.null(nonhuman_results) && !is.null(nonhuman_results$merged)
    
    # Create reference text instead of displaying the table
    ref_text <- paste0(
      "Kraken2 classification results have been compiled and saved as CSV files.\n\n",
      "The Kraken2 results data includes:\n",
      "- Sample names\n",
      "- Species classifications\n",
      "- Unaligned clade reads\n"
    )
    
    if (has_nonhuman) {
      ref_text <- paste0(ref_text, "- Non-human clade reads\n")
      ref_text <- paste0(ref_text, "- Dataset information (unaligned/nonhuman)\n")
    } else {
      ref_text <- paste0(ref_text, "- Dataset information (unaligned only - no nonhuman data available)\n")
    }
    
    ref_text <- paste0(ref_text,
      "- Merged input read counts from RRstats\n",
      "- Additional optional columns (confidence, minimum hit groups, etc.)\n\n",
      "Data source: Combined from unaligned_kreports/", if (has_nonhuman) " and nonhuman_kreports/" else "", "\n",
      "Processing: output_processing.R with species filtering and merging\n",
      "Merging: Combined with runtime and read statistics by sample name\n",
      "Output location: Check outputs/ directory for:\n",
      "  - sample_report_data_with_metadata*op.csv (combined data)\n",
      "  - sample_report_data_[dataset]_*op.csv (dataset-specific files)\n\n",
      "Key findings:\n",
      "- Average number of species found per sample before secondary filtering: ", mean(combined_report_data$species_count, na.rm = TRUE), "\n",
      "- Average number of species found per sample after secondary filtering: ", mean(combined_report_data$filtered_species_count, na.rm = TRUE), "\n",
      "- Number of unique species across all unaligned samples: ", nrow(unaligned_results$species_list), "\n"
    )
    
    if (config$INCLUDE_SUBSPECIES) {
      ref_text <- paste(ref_text, "- Average subspecies across all samples before secondary filtering: ", mean(combined_report_data$subspecies_count), "\n", sep = "")
      ref_text <- paste(ref_text, "- Average subspecies across all samples after secondary filtering: ", mean(combined_report_data$subspecies_count), "\n", sep = "")
    }
    
    top_clades <- unaligned_results$merged %>% 
      arrange(desc(cladeReads_mean)) %>%
      slice_head(n = 10)
    ref_text <- paste(ref_text, paste(top_clades$name, top_clades$cladeReads_mean, sep = ": ", collapse = "\n"), sep = "\n")
    
    # Create summary plot using ggplot
    summary_plot <- ggplot() + 
      theme_void() + 
      labs(title = title_text, subtitle = ref_text)
    print(summary_plot)
    
    if (has_nonhuman) {
      cat("Creating Kraken2 statistics report using non-human sample results...\n")
      species_list <- nonhuman_results$species_list
    } else {
      species_list <- unaligned_results$species_list
    }
    
    # Assign to global environment for the plot scripts to access
    assign("species_list", species_list, envir = .GlobalEnv)
    
    # Source the function-based plotting (disable main execution)
    assign("..main_execution_disabled..", TRUE, envir = .GlobalEnv)
    source(file.path(config$RANALYSIS_DIR, "scripts", "Stage1_Data_acquisition", "annotate_species_plots.R"))
    
    # Kingdom basic plot
    kingdom_basic <- generate_annotation_plot(
      df = species_list,
      plot_type = "kingdom"
    )
    print(kingdom_basic)
    
    # Kingdom detailed plots
    detailed_kingdom_plot <- generate_annotation_plot(
      df = species_list,
      plot_type = "kingdom",
      detailed = "Freq",
      log = TRUE
    )
    print(detailed_kingdom_plot)
    
    # Check if bracken data is available for advanced plotting
    if (has_nonhuman) {
      bracken_cols <- grep("bracken_reads_mean", colnames(nonhuman_results$species_list), value = TRUE)
    } else {
      # For unaligned results, check the same way
      bracken_cols <- grep("bracken_reads_mean", colnames(unaligned_results$species_list), value = TRUE)
    }
    
    if (length(bracken_cols) > 0) {
      # First: Show log mean bracken reads kingdom plot
      bracken_kingdom_plot <- generate_annotation_plot(
        df = species_list,
        plot_type = "kingdom",
        detailed = "bracken_reads_mean",
        log = TRUE
      )
      print(bracken_kingdom_plot)
      
      # Second: Create side-by-side comparison plot (clade reads vs bracken reads)
      cat("Creating side-by-side clade vs bracken reads comparison plot...\n")
      
      # Calculate dynamic N based on config$TOP_N_PLOT
      dynamic_n <- max(1, floor(config$TOP_N_PLOT / 2))
      
      # Prepare data for side-by-side comparison (ranked left-to-right by species)
      top_species_comparison <- unaligned_results$species_list %>%
        filter(!is.na(cladeReads_mean) & !is.na(bracken_reads_mean)) %>%
        arrange(desc(cladeReads_mean)) %>%
        slice_head(n = dynamic_n) %>%
        mutate(species_rank = row_number()) %>%
        select(name, cladeReads_mean, bracken_reads_mean, species_rank)
      
      if (nrow(top_species_comparison) > 0) {
        # Transform data for side-by-side plotting
        comparison_long <- top_species_comparison %>%
          pivot_longer(cols = c("cladeReads_mean", "bracken_reads_mean"), 
                       names_to = "read_type", values_to = "reads") %>%
          mutate(read_type = case_when(
            read_type == "cladeReads_mean" ~ "Kraken2 Clade",
            read_type == "bracken_reads_mean" ~ "Bracken"
          )) %>%
          mutate(read_type = factor(read_type, levels = c("Kraken2 Clade", "Bracken")))
        
        # Create side-by-side bar plot with improved ranking and log scale
        comparison_plot <- ggplot(comparison_long, 
                                 aes(x = reorder(name, species_rank), y = log10(reads + 1), fill = read_type)) +
          geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
          scale_fill_manual(values = c("Kraken2 Clade" = "steelblue", "Bracken" = "orange")) +
          scale_y_continuous(trans = "log10", 
                             labels = function(x) format(10^x, scientific = FALSE)) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
                legend.position = "bottom") +
          labs(
            title = paste0("Kraken2 Clade Reads vs Bracken Reads (Top ", dynamic_n, " Species)"),
            subtitle = paste("Ranked by Kraken2 clade reads (log10 scale)"),
            x = "Species (ranked by clade reads)",
            y = "Log10(Reads + 1)",
            fill = "Read Type"
          )
        
        print(comparison_plot)
      } else {
        plot.new()
        text(0.5, 0.5, "Insufficient data for Bracken comparison plot", 
             adj = c(0.5, 0.5), cex = 1, col = "gray50")
      }
    } else {
      # No bracken data - show the original clade reads mean plot
      mean_kingdom_plot <- generate_annotation_plot(
        df = species_list,
        plot_type = "kingdom",
        detailed = "cladeReads_mean",
        log = TRUE
      )
      print(mean_kingdom_plot)
      
      # Create placeholder plot for missing bracken data
      placeholder_plot <- ggplot() + 
        theme_void() + 
        labs(title = "Bracken Data", subtitle = "No Bracken data available for comparison")
      print(placeholder_plot)
    }
    
    # Add percentage plots using calculate_species_percentages
    cat("Creating percentage-based plots...\n")
    tryCatch({
      # Determine which merged_long data to use for percentage calculations
      merged_long_data <- if (has_nonhuman && !is.null(nonhuman_results$merged_long)) {
        nonhuman_results$merged_long
      } else {
        unaligned_results$merged_long
      }
      
      # Calculate species percentages
      species_percentages <- calculate_species_percentages(merged_long_data)
      
      if (nrow(species_percentages) > 0) {
        # Create percentage-based kingdom plot
        # Add percentage column to species_list for plotting
        species_list_with_pct <- species_list %>%
          left_join(
            species_percentages %>%
              group_by(name) %>%
              summarise(avg_percentage = mean(percentage, na.rm = TRUE), .groups = 'drop'),
            by = "name"
          ) %>%
          mutate(avg_percentage = ifelse(is.na(avg_percentage), 0, avg_percentage))
        
        percentage_kingdom_plot <- generate_annotation_plot(
          df = species_list_with_pct,
          plot_type = "kingdom",
          detailed = "avg_percentage"
        )
        print(percentage_kingdom_plot)
      } else {
        # Create placeholder plot for missing percentage data
        placeholder_pct_plot <- ggplot() + 
          theme_void() + 
          labs(title = "Species Percentages", subtitle = "No percentage data available")
        print(placeholder_pct_plot)
      }
    }, error = function(e) {
      cat("Warning: Error creating percentage plots:", e$message, "\n")
      error_plot <- ggplot() + 
        theme_void() + 
        labs(title = "Species Percentages", subtitle = "Error calculating percentages")
      print(error_plot)
    })
    
    dev.off()
    
    # Add new columns to combined_report_data and export to Excel
    tryCatch({
      # Export to Excel with conditional formatting if openxlsx is available
      if (excel_available) {
        excel_file <- file.path(report_dir, "kraken_stats_report.xlsx")
        
        # Create workbook
        wb <- createWorkbook()
        addWorksheet(wb, "Kraken_Stats")
        
        # Create headers without underscores
        headers <- colnames(combined_report_data)
        clean_headers <- gsub("_", " ", headers)
        clean_headers <- tools::toTitleCase(clean_headers)
        
        # Write headers and data separately
        writeData(wb, "Kraken_Stats", t(clean_headers), startRow = 1, colNames = FALSE)
        writeData(wb, "Kraken_Stats", combined_report_data, startRow = 2, colNames = FALSE)
        
        # Add conditional formatting
        data_rows <- 2:(nrow(combined_report_data) + 1)
        
        # Red-white color scale for unclassified counts
        unclassified_patterns <- c("unclassified_counts", "unclassified.count", "unclassified")
        unclassified_col <- NA
        for (pattern in unclassified_patterns) {
          if (pattern %in% headers) {
            unclassified_col <- which(headers == pattern)
            break
          }
        }
        
        if (!is.na(unclassified_col)) {
          conditionalFormatting(wb, "Kraken_Stats", 
                              cols = unclassified_col, 
                              rows = data_rows,
                              type = "colorScale", 
                              style = c("#FFFFFF", "#FF0000"))
          cat("Applied red-white color scale to unclassified column\n")
        }
        
        # Data bars for species_count
        species_patterns <- c("species_count", "species.count", "species")
        species_col <- NA
        for (pattern in species_patterns) {
          if (pattern %in% headers) {
            species_col <- which(headers == pattern)
            break
          }
        }
        
        if (!is.na(species_col)) {
          conditionalFormatting(wb, "Kraken_Stats", 
                              cols = species_col, 
                              rows = data_rows,
                              type = "dataBar",
                              style = "#0070C0")
          cat("Applied data bars to species_count column\n")
        }
        
        # Auto-size columns
        setColWidths(wb, "Kraken_Stats", cols = 1:length(headers), widths = "auto")
        
        # Add some basic formatting to headers
        header_style <- createStyle(
          textDecoration = "bold",
          fgFill = "#E6E6FA",
          border = "TopBottomLeftRight",
          borderStyle = "thin"
        )
        addStyle(wb, "Kraken_Stats", header_style, rows = 1, cols = 1:length(headers))
        border_style <- createStyle(border = "TopBottomLeftRight", borderStyle = "thin")
        # Add borders to all data cells
        for (col in 1:length(headers)) {
          addStyle(wb, "Kraken_Stats", border_style, rows = data_rows, cols = col)
        }
        
        # Save workbook
        saveWorkbook(wb, excel_file, overwrite = TRUE)
        cat("Kraken statistics exported to Excel:", excel_file, "\n")
      } else {
        cat("Excel export skipped - openxlsx not available\n")
      }
    }, error = function(e) {
      cat("Warning: Failed to modify data or export to Excel:", e$message, "\n")
      print(e)
    })
    
  } else {
    # Per-sample mode - generate plots for individual sample
    sample_merged_long <- unaligned_results  # This is the per-sample data
    sample_name <- sample_merged_long$sample[1]
    kraken_plots <- list()
    if (config$ADD_PARAMS) {
      has_nonhuman <- sample_merged_long$human_reads[1]
    } else {
      has_nonhuman <- nonhuman_results
    }
    
    # Summary plot
    title_text <- paste("Part 3: Kraken2 Statistics -", sample_name)
    
    # Create reference text for individual sample
    ref_text <- paste0(
      "The Kraken2 results for sample ", sample_name, ":\n"
    )
    
    if (has_nonhuman) {
      ref_text <- paste0(ref_text, "- Using Non-human sample results\n")
    } else {
      ref_text <- paste0(ref_text, "- Using Unaligned sample results only (no nonhuman data available)\n")
    }

    ref_text <- paste0(ref_text,
      "Key findings:\n",
      "- Number of species found before secondary filtering: ", combined_report_data$species_count[1], "\n",
      "- Number of species found after secondary filtering: ", combined_report_data$filtered_species_count[1], "\n"
    )

    if (config$INCLUDE_SUBSPECIES) {
      ref_text <- paste0(ref_text, 
        "- Subspecies before secondary filtering: ", combined_report_data$subspecies_count[1], "\n",
        "- Subspecies after secondary filtering: ", combined_report_data$filtered_subspecies_count[1], "\n"
      )
    }

    top_clades <- sample_merged_long %>% 
      filter(!is.na(cladeReads)) %>%
      arrange(desc(cladeReads)) %>%
      slice_head(n = min(config$TOP_N_FREQ, nrow(.)))

    ref_text <- paste0(ref_text, 
      "\nTop species in this sample:\n",
      paste(top_clades$name, top_clades$cladeReads, sep = ": ", collapse = "\n")
    )
    
    summary_plot <- ggplot() + 
      theme_void() + 
      labs(title = title_text, subtitle = ref_text)
    kraken_plots[[paste0(sample_name, "_summary")]] <- summary_plot
    
    # Filter out TaxID 9606 (human) and calculate percentages for kingdom plot
    sample_filtered <- sample_merged_long %>%
      filter(taxID != 9606 & !is.na(cladeReads))
    
    if (nrow(sample_filtered) > 0) {
      # Use calculate_species_percentages for consistent percentage calculation
      tryCatch({
        species_percentages <- calculate_species_percentages(sample_filtered)
        
        # Add percentage data back to sample_filtered
        if (nrow(species_percentages) > 0) {
          sample_filtered <- sample_filtered %>%
            left_join(
              species_percentages %>% select(name, percentage),
              by = "name"
            ) %>%
            mutate(percentage = ifelse(is.na(percentage), 0, percentage))
        } else {
          # Fallback to manual calculation if function fails
          read_col <- if ("bracken_reads" %in% colnames(sample_filtered) && 
                         any(!is.na(sample_filtered$bracken_reads))) {
            "bracken_reads"
          } else {
            "cladeReads"
          }
          total_reads <- sum(sample_filtered[[read_col]], na.rm = TRUE)
          sample_filtered$percentage <- (sample_filtered[[read_col]] / total_reads) * 100
        }
      }, error = function(e) {
        cat("Warning: Error using calculate_species_percentages, falling back to manual calculation:", e$message, "\n")
        # Fallback calculation
        read_col <- if ("bracken_reads" %in% colnames(sample_filtered) && 
                       any(!is.na(sample_filtered$bracken_reads))) {
          "bracken_reads"
        } else {
          "cladeReads"
        }
        total_reads <- sum(sample_filtered[[read_col]], na.rm = TRUE)
        sample_filtered$percentage <- (sample_filtered[[read_col]] / total_reads) * 100
      })
      
      # Source the function-based plotting instead of run_as_script (disable main execution)
      assign("..main_execution_disabled..", TRUE, envir = .GlobalEnv)
      source(file.path(config$RANALYSIS_DIR, "scripts", "Stage1_Data_acquisition", "annotate_species_plots.R"))
      
      # Kingdom basic plot using function interface
      kingdom_basic <- generate_annotation_plot(
        df = sample_filtered,
        plot_type = "kingdom",
        sample_name = sample_name
      )
      kraken_plots[[paste0(sample_name, "_kingdom_basic")]] <- kingdom_basic
      
      # Kingdom detailed plot using read counts (log scale)
      read_col <- if ("bracken_reads" %in% colnames(sample_filtered) && 
                     any(!is.na(sample_filtered$bracken_reads))) {
        "bracken_reads"
      } else {
        "cladeReads"
      }
      
      detailed_kingdom_plot <- generate_annotation_plot(
        df = sample_filtered,
        plot_type = "kingdom",
        detailed = read_col,
        sample_name = sample_name,
        top_n = config$TOP_N_PLOT,
        log = TRUE
      )
      kraken_plots[[paste0(sample_name, "_kingdom_detailed")]] <- detailed_kingdom_plot
      
      # Kingdom percentage plot
      percentage_kingdom_plot <- generate_annotation_plot(
        df = sample_filtered,
        plot_type = "kingdom",
        detailed = "percentage",
        sample_name = sample_name,
        top_n = config$TOP_N_PLOT
      )
      kraken_plots[[paste0(sample_name, "_kingdom_percentage")]] <- percentage_kingdom_plot
      
      # Add bracken comparison plot if bracken data is available
      if ("bracken_reads" %in% colnames(sample_filtered) && 
          any(!is.na(sample_filtered$bracken_reads))) {
        tryCatch({
          # Create comparison data for bracken vs cladeReads
          comparison_data <- sample_filtered %>%
            select(name, cladeReads, bracken_reads) %>%
            filter(!is.na(cladeReads) & !is.na(bracken_reads)) %>%
            arrange(desc(cladeReads)) %>%
            slice_head(n = min(config$TOP_N_PLOT, nrow(.))) %>%
            mutate(species_rank = row_number()) %>%
            pivot_longer(cols = c("cladeReads", "bracken_reads"), 
                         names_to = "read_type", values_to = "reads") %>%
            mutate(read_type = case_when(
              read_type == "cladeReads" ~ "Kraken2 Clade",
              read_type == "bracken_reads" ~ "Bracken"
            )) %>%
            mutate(read_type = factor(read_type, levels = c("Kraken2 Clade", "Bracken")))
          
          if (nrow(comparison_data) > 0) {
            # Create side-by-side bar plot with log scale
            bracken_comparison_plot <- ggplot(comparison_data, 
                                             aes(x = reorder(name, species_rank), y = log10(reads + 1), fill = read_type)) +
              geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
              scale_fill_manual(values = c("Kraken2 Clade" = "steelblue", "Bracken" = "orange")) +
              scale_y_continuous(trans = "log10", 
                                 labels = function(x) format(10^x, scientific = FALSE)) +
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
                    legend.position = "bottom") +
              labs(
                title = paste0("Kraken2 vs Bracken Reads - ", sample_name),
                subtitle = paste("Top species (log10 scale)"),
                x = "Species (ranked by clade reads)",
                y = "Log10(Reads + 1)",
                fill = "Read Type"
              )
            
            kraken_plots[[paste0(sample_name, "_bracken_comparison")]] <- bracken_comparison_plot
          }
        }, error = function(e) {
          cat("Warning: Error creating bracken comparison plot for", sample_name, ":", e$message, "\n")
        })
      }
      
    } else {
      # No data available - create placeholder plots
      empty_plot <- ggplot() + theme_void() + labs(title = paste("No data available for", sample_name))
      
      kraken_plots[[paste0(sample_name, "_kingdom_basic")]] <- empty_plot
      kraken_plots[[paste0(sample_name, "_kingdom_detailed")]] <- empty_plot
    }
    
    return(kraken_plots)
  }
}

species_annotation_report <- function(species_annotation_data, report_dir, batch = TRUE, config) {
  if (batch) {
    pdf(file.path(report_dir, "species_annotation_report.pdf"), width = 8, height = 6)
    
    title_text <- "Part 4: Species Annotation"
    
    # Create reference text instead of displaying the table
    ref_text <- paste0(
      "Species annotation results have been compiled and saved as CSV files.\n\n",
      "The species annotation data includes:\n",
      "- Sample names\n",
      "- Annotated species names\n",
      "- TaxIDs and clade information\n",
      "- Additional metadata columns (if available)\n\n",
      "Data source: Parsed from annotated_species object\n",
      "Output location: Check outputs/ directory for:\n",
      "  - annotated_species.csv (annotated species list)\n\n",
      "Key findings:\n",
      "- Total annotated species: ", nrow(species_annotation_data), "\n"
    )
    
    # Calculate Risk Group summary
    risk_group_summary <- species_annotation_data %>%
      count(RiskGroup, name = "Number_of_species") %>%
      mutate(
        RiskGroup = ifelse(is.na(RiskGroup) | RiskGroup == "", "NotAnnotated", RiskGroup),
        Percentage = round(100 * Number_of_species / sum(Number_of_species), 1)
      ) %>%
      arrange(RiskGroup)
    
    # Add Risk Group summary to reference text
    ref_text <- paste0(ref_text, "\nRisk Group Summary:\n")
    for (i in 1:nrow(risk_group_summary)) {
      ref_text <- paste0(ref_text, 
                        risk_group_summary$RiskGroup[i], ": ", 
                        risk_group_summary$Number_of_species[i], " species (", 
                        risk_group_summary$Percentage[i], "%)\n")
    }
    
    # Create summary plot using ggplot
    summary_plot <- ggplot() + 
      theme_void() + 
      labs(title = title_text, subtitle = ref_text)
    print(summary_plot)
    
    # Generate Summary Plots
    tryCatch({
      # Source the function-based plotting (disable main execution)
      assign("..main_execution_disabled..", TRUE, envir = .GlobalEnv)
      source(file.path(config$RANALYSIS_DIR, "scripts", "Stage1_Data_acquisition", "annotate_species_plots.R"))
      
      # HOMD basic plot
      homd_basic <- generate_annotation_plot(
        df = species_annotation_data,
        plot_type = "homd"
      )
      print(homd_basic)
      
      # Risk group basic plot
      risk_group_basic <- generate_annotation_plot(
        df = species_annotation_data,
        plot_type = "risk_group"
      )
      print(risk_group_basic)
      
      # Summary basic plot
      summary_basic <- generate_annotation_plot(
        df = species_annotation_data,
        plot_type = "summary"
      )
      print(summary_basic)
      
    }, error = function(e) {
      cat("Warning: Error generating annotation plots:", e$message, "\n")
      plot.new()
      text(0.5, 0.5, "Annotation plots could not be generated", cex = 1.2, adj = 0.5)
    })
    
    # Generate Detailed Plots
    tryCatch({
      # Detailed HOMD plot
      detailed_homd_plot <- generate_annotation_plot(
        df = species_annotation_data,
        plot_type = "homd",
        detailed = "Freq"
      )
      # Add frequency subtitle to HOMD plot
      detailed_homd_plot <- detailed_homd_plot + 
        labs(subtitle = paste0(detailed_homd_plot$labels$subtitle, 
                              "\nFrequency in top ", config$TOP_N_FREQ, " species of samples"))
      print(detailed_homd_plot)
      
      # Detailed risk group plot
      detailed_risk_group_plot <- generate_annotation_plot(
        df = species_annotation_data,
        plot_type = "risk_group",
        detailed = "Freq"
      )
      # Add frequency subtitle to Risk Group plot
      detailed_risk_group_plot <- detailed_risk_group_plot + 
        labs(subtitle = paste0(detailed_risk_group_plot$labels$subtitle, 
                              "\nFrequency in top ", config$TOP_N_FREQ, " species of samples"))
      print(detailed_risk_group_plot)
      
    }, error = function(e) {
      cat("Warning: Error generating detailed annotation plots:", e$message, "\n")
      error_plot <- ggplot() + 
        theme_void() + 
        labs(title = "Detailed Annotation Plots", subtitle = "Detailed annotation plots could not be generated")
      print(error_plot)
    })
    
    dev.off()
    
    # Export annotated species to Excel with formatting
    tryCatch({
      if (excel_available) {
        excel_file <- file.path(report_dir, "species_annotation_report.xlsx")
        
        # Create workbook
        wb <- createWorkbook()
        addWorksheet(wb, "Species_Annotation")
        
        # Create headers without underscores and periods
        headers <- colnames(species_annotation_data)
        clean_headers <- gsub("[._]", " ", headers)
        clean_headers <- tools::toTitleCase(clean_headers)
        
        # Prepare data for Excel export - ensure proper data types
        excel_data <- species_annotation_data
        
        # Convert numeric columns to ensure they're properly formatted
        numeric_cols <- c("Freq", "cladeReads_mean", "cladeReads_max", "bracken_reads_mean", "bracken_reads_max")
        for (col in numeric_cols) {
          if (col %in% colnames(excel_data)) {
            excel_data[[col]] <- as.numeric(excel_data[[col]])
          }
        }
        
        # Ensure character columns are properly formatted
        character_cols <- c("name", "taxID", "taxRank", "taxLineage", "RiskGroup", "HOMD", "HOMD.Category")
        for (col in character_cols) {
          if (col %in% colnames(excel_data)) {
            excel_data[[col]] <- as.character(excel_data[[col]])
          }
        }
        
        # Write headers and data separately
        writeData(wb, "Species_Annotation", t(clean_headers), startRow = 1, colNames = FALSE)
        writeData(wb, "Species_Annotation", excel_data, startRow = 2, colNames = FALSE)
        
        # Auto-size columns
        setColWidths(wb, "Species_Annotation", cols = 1:length(headers), widths = "auto")
        
        # Add header formatting
        header_style <- createStyle(
          textDecoration = "bold",
          fgFill = "#E6E6FA",
          border = "TopBottomLeftRight",
          borderStyle = "thin"
        )
        addStyle(wb, "Species_Annotation", header_style, rows = 1, cols = 1:length(headers))
        
        # Add borders to data cells (apply to each row individually)
        data_rows <- 2:(nrow(excel_data) + 1)
        data_style <- createStyle(
          border = "TopBottomLeftRight",
          borderStyle = "thin"
        )
        for (row in data_rows) {
          addStyle(wb, "Species_Annotation", data_style, rows = row, cols = 1:length(headers))
        }
        
        # Save workbook
        saveWorkbook(wb, excel_file, overwrite = TRUE)
        cat("Species annotation exported to Excel:", excel_file, "\n")
      } else {
        cat("Excel export skipped - openxlsx not available\n")
      }
    }, error = function(e) {
      cat("Warning: Failed to export species annotation to Excel:", e$message, "\n")
      print(e)
    })
    
  } else {
    # Per-sample mode - generate plots for individual sample
    sample_data <- species_annotation_data  # This is the per-sample data
    sample_name <- sample_data$sample[1]
    
    # Save sample data to CSV if report_dir is provided (per-sample mode with SPLIT_REPORT)
    if (!is.null(report_dir)) {
      sample_annotated_csv <- file.path(report_dir, paste0(sample_name, "_annotated_species.csv"))
      write.csv(sample_data, sample_annotated_csv, row.names = FALSE)
      cat("Sample annotated species saved to:", sample_annotated_csv, "\n")
    }
    
    # Summary plot
    title_text <- paste("Part 4: Species Annotation -", sample_name)
    
    ref_text <- paste0(
      "Species annotation results for sample ", sample_name, ":\n\n",
      "Key findings:\n",
      "- Total species in sample: ", nrow(sample_data), "\n",
      "- Species with risk annotations: ", sum(sample_data$RiskGroup != "NotAnnotated", na.rm = TRUE), "\n",
      "- Species with HOMD annotations: ", sum(sample_data$HOMD.Category != "NotAnnotated", na.rm = TRUE), "\n",
      "- High-risk species (RG3/4): ", sum(sample_data$RiskGroup %in% c("RG3", "RG4"), na.rm = TRUE), "\n",
      "- Pathogenic species (HOMD): ", sum(sample_data$HOMD.Category %in% c("Pathogen", "Opportunist"), na.rm = TRUE), "\n"
    )
    
    summary_plot <- ggplot() + 
      theme_void() + 
      labs(title = title_text, subtitle = ref_text)
    annotation_plots <- list()
    annotation_plots[[paste0(sample_name, "_summary")]] <- summary_plot
    
    # Source the function-based plotting (disable main execution)
    assign("..main_execution_disabled..", TRUE, envir = .GlobalEnv)
    source(file.path(config$RANALYSIS_DIR, "scripts", "Stage1_Data_acquisition", "annotate_species_plots.R"))
    
    sample_annotated <- sample_data
    # Generate annotation plots if we have annotated data
    if (nrow(sample_annotated) > 0) {
      tryCatch({
        # HOMD basic plot
        homd_basic <- generate_annotation_plot(
          df = sample_annotated,
          plot_type = "homd",
          sample_name = sample_name
        )
        annotation_plots[[paste0(sample_name, "_homd_basic")]] <- homd_basic
        
        # Risk group basic plot
        risk_group_basic <- generate_annotation_plot(
          df = sample_annotated,
          plot_type = "risk_group", 
          sample_name = sample_name
        )
        annotation_plots[[paste0(sample_name, "_risk_group_basic")]] <- risk_group_basic
        
        # Summary basic plot
        summary_basic <- generate_annotation_plot(
          df = sample_annotated,
          plot_type = "summary",
          sample_name = sample_name
        )
        annotation_plots[[paste0(sample_name, "_summary_basic")]] <- summary_basic
        
        # Check for read data for detailed plots
        read_col <- if ("bracken_reads" %in% colnames(sample_annotated) && 
                        any(!is.na(sample_annotated$bracken_reads))) {
          "bracken_reads"
        } else if ("cladeReads" %in% colnames(sample_annotated) && 
                   any(!is.na(sample_annotated$cladeReads))) {
          "cladeReads"
        } else {
          NULL
        }
        
        if (!is.null(read_col)) {
          # Detailed HOMD plot
          detailed_homd_plot <- generate_annotation_plot(
            df = sample_annotated,
            plot_type = "homd",
            detailed = read_col,
            sample_name = sample_name,
            top_n = config$TOP_N_PLOT,
            log = TRUE
          )
          annotation_plots[[paste0(sample_name, "_homd_detailed")]] <- detailed_homd_plot
          
          # Detailed risk group plot
          detailed_risk_group_plot <- generate_annotation_plot(
            df = sample_annotated,
            plot_type = "risk_group",
            detailed = read_col,
            sample_name = sample_name,
            top_n = config$TOP_N_PLOT,
            log = TRUE
          )
          annotation_plots[[paste0(sample_name, "_risk_group_detailed")]] <- detailed_risk_group_plot
          
          # Add percentage plots using calculate_species_percentages
          species_percentages <- calculate_species_percentages(sample_annotated)
          
          if (nrow(species_percentages) > 0) {
            # Add percentage data back to sample_annotated for plotting
            sample_annotated_with_pct <- sample_annotated %>%
              left_join(
                species_percentages %>% select(name, percentage),
                by = "name"
              ) %>%
              mutate(percentage = ifelse(is.na(percentage), 0, percentage))
            
            # Percentage HOMD plot
            percentage_homd_plot <- generate_annotation_plot(
              df = sample_annotated_with_pct,
              plot_type = "homd",
              detailed = "percentage",
              sample_name = sample_name,
              top_n = config$TOP_N_PLOT
            )
            annotation_plots[[paste0(sample_name, "_homd_percentage")]] <- percentage_homd_plot
            
            # Percentage risk group plot
            percentage_risk_group_plot <- generate_annotation_plot(
              df = sample_annotated_with_pct,
              plot_type = "risk_group",
              detailed = "percentage",
              sample_name = sample_name,
              top_n = config$TOP_N_PLOT
            )
            annotation_plots[[paste0(sample_name, "_risk_group_percentage")]] <- percentage_risk_group_plot
          } else {
            cat("Warning: No percentage data available for", sample_name, "\n")
          }
        }
      }, error = function(e) {
        cat("Warning: Error generating annotation plots for", sample_name, ":", e$message, "\n")
        empty_plot <- ggplot() + 
          theme_void() + 
          labs(title = paste("Annotation plots not available for", sample_name))
        annotation_plots[[paste0(sample_name, "_error")]] <- empty_plot
      })
    } else {
      # No annotated data available
      empty_plot <- ggplot() + 
        theme_void() + 
        labs(title = paste("No annotation data available for", sample_name))
      annotation_plots[[paste0(sample_name, "_no_data")]] <- empty_plot
    }
    
    return(annotation_plots)
  }
}

pathogen_detection_report <- function(pathogen_annotation_data, merged_long_data, report_dir, batch, config) {
  if (batch) {
    pdf(file.path(report_dir, "pathogen_detection_report.pdf"), width = 8, height = 6)
    
    title_text <- "Part 5: Pathogen Detection (RG3/4 Species)"
    
    # Filter for pathogenic species (RG3/4 and non-microbiome)
    pathogen_species <- pathogen_annotation_data %>%
      filter(
        RiskGroup %in% c("RG3", "RG4") | 
        (is.na(RiskGroup) & !HOMD.Category %in% c("Microbiome", "microbiome", "MICROBIOME")) |
        (!RiskGroup %in% c("RG1", "RG2") & !HOMD.Category %in% c("Microbiome", "microbiome", "MICROBIOME"))
      ) %>%
      filter(!is.na(name))
    
    # Create reference text
    ref_text <- paste0(
      "Pathogen detection analysis focuses on RG3/4 species and non-microbiome organisms.\n\n",
      "Filtering criteria:\n",
      "- Include: RG3, RG4, risk groups and NotAnnotated species\n",
      "- Exclude: RG1, RG2, risk groups\n",
      "- Exclude: Species annotated as 'Microbiome' in HOMD.Category\n\n",
      "Key findings:\n",
      "- Total species after annotation: ", nrow(pathogen_annotation_data), "\n",
      "- Potential pathogenic species identified: ", nrow(pathogen_species), "\n"
    )
    
    if (nrow(pathogen_species) > 0) {
      ref_text <- paste0(ref_text,
        "- Most abundant pathogen: ", pathogen_species$name[which.max(pathogen_species$cladeReads_mean)], "\n",
        "- Average pathogen clade reads: ", round(mean(pathogen_species$cladeReads_mean, na.rm = TRUE), 1), "\n"
      )
    }
    
    ref_text <- paste0(ref_text,
      "\nData source: Filtered from annotated_species object\n",
      "Output location: Check outputs/ directory for:\n",
      "  - pathogen_detection_report.xlsx (pathogen summary by sample)\n"
    )
    
    # Create summary plot using ggplot
    summary_plot <- ggplot() + 
      theme_void() + 
      labs(title = title_text, subtitle = ref_text)
    print(summary_plot)
    
    # Generate the same annotation plots but filtered for pathogens only
    if (nrow(pathogen_species) > 0) {
      tryCatch({
        # Source the function-based plotting (disable main execution)
        assign("..main_execution_disabled..", TRUE, envir = .GlobalEnv)
        source(file.path(config$RANALYSIS_DIR, "scripts", "Stage1_Data_acquisition", "annotate_species_plots.R"))
        
        # Risk group plot for pathogens
        pathogen_risk_plot <- generate_annotation_plot(
          df = pathogen_species,
          plot_type = "risk_group"
        )
        print(pathogen_risk_plot)
        
        # Detailed risk group plot
        detailed_pathogen_risk_plot <- generate_annotation_plot(
          df = pathogen_species,
          plot_type = "risk_group",
          detailed = "cladeReads_mean",
          log = TRUE
        )
        print(detailed_pathogen_risk_plot)
        
      }, error = function(e) {
        cat("Warning: Error generating pathogen plots:", e$message, "\n")
        plot.new()
        text(0.5, 0.5, "Pathogen plots could not be generated", cex = 1.2, adj = 0.5)
      })
    } else {
      plot.new()
      text(0.5, 0.5, "No pathogenic species identified", cex = 1.5, adj = 0.5, col = "darkgreen")
    }
    
    dev.off()
    
    # Generate Excel report with sample-level pathogen percentages
    tryCatch({
      if (excel_available && nrow(pathogen_species) > 0) {
        # Calculate sample-level percentages for pathogenic species only
        pathogen_sample_data <- calculate_pathogen_percentages(merged_long_data, pathogen_species)
        
        # Add RiskGroup annotation to pathogen sample data
        pathogen_sample_annotation_data <- pathogen_sample_data %>%
          left_join(pathogen_species %>% select(name, RiskGroup), 
                   by = c("RG3/4 and NotAnnotated Species" = "name"))
        
        if (nrow(pathogen_sample_annotation_data) > 0) {
          excel_file <- file.path(report_dir, "pathogen_detection_report.xlsx")
          
          # Create workbook
          wb <- createWorkbook()
          addWorksheet(wb, "Pathogen_Detection")
          
          # Create clean headers
          headers <- c("Sample Name", "RG3/4 and NotAnnotated Species", "Percentage of Reads in Sample", "Risk Group")
          
          # Write headers and data
          writeData(wb, "Pathogen_Detection", t(headers), startRow = 1, colNames = FALSE)
          writeData(wb, "Pathogen_Detection", pathogen_sample_annotation_data, startRow = 2, colNames = FALSE)
          
          # Auto-size columns
          setColWidths(wb, "Pathogen_Detection", cols = 1:length(headers), widths = "auto")
          
          # Add header formatting
          header_style <- createStyle(
            textDecoration = "bold",
            fgFill = "#FFE6E6",  # Light red for pathogen report
            border = "TopBottomLeftRight",
            borderStyle = "thin"
          )
          addStyle(wb, "Pathogen_Detection", header_style, rows = 1, cols = 1:length(headers))
          
          # Add borders to data cells
          data_rows <- 2:(nrow(pathogen_sample_annotation_data) + 1)
          data_style <- createStyle(
            border = "TopBottomLeftRight",
            borderStyle = "thin"
          )
          for (row in data_rows) {
            addStyle(wb, "Pathogen_Detection", data_style, rows = row, cols = 1:length(headers))
          }
          
          # Add conditional formatting for high percentages
          if (nrow(pathogen_sample_annotation_data) > 1) {
            conditionalFormatting(wb, "Pathogen_Detection", 
                                cols = 3, 
                                rows = data_rows,
                                type = "colorScale", 
                                style = c("#FFFFFF", "#FF6B6B"))
          }
          
          addWorksheet(wb, "Pathogen_Annotations")
          writeData(wb, "Pathogen_Annotations", pathogen_species, startRow = 1, colNames = TRUE)
          
          # Save workbook
          saveWorkbook(wb, excel_file, overwrite = TRUE)
          cat("Pathogen detection report exported to Excel:", excel_file, "\n")
        } else {
          cat("No pathogen detections found in samples\n")
        }
      } else if (!excel_available) {
        cat("Excel export skipped - openxlsx not available\n")
      } else {
        cat("No pathogenic species to report\n")
      }
    }, error = function(e) {
      cat("Warning: Failed to export pathogen detection to Excel:", e$message, "\n")
      print(e)
    })
    return(pathogen_species)
  } else {
    # Per-sample mode - generate plots for individual sample
    sample_merged_long <- merged_long_data  # This is the per-sample data
    sample_name <- sample_merged_long$sample[1]
    pathogen_plots <- list()
    
    # Use the pathogen_annotation_data parameter which already has annotations merged
    sample_annotated <- pathogen_annotation_data %>%
      # Fill missing annotations
      mutate(
        RiskGroup = ifelse(is.na(RiskGroup), "NotAnnotated", RiskGroup),
        HOMD.Category = ifelse(is.na(HOMD.Category), "NotAnnotated", HOMD.Category)
      )
    
    # Filter for pathogenic species (same criteria as batch mode)
    sample_pathogen_species <- sample_annotated %>%
      filter(
        RiskGroup %in% c("RG3", "RG4") | 
        (is.na(RiskGroup) & !HOMD.Category %in% c("Microbiome", "microbiome", "MICROBIOME")) |
        (!RiskGroup %in% c("RG1", "RG2") & !HOMD.Category %in% c("Microbiome", "microbiome", "MICROBIOME"))
      ) %>%
      filter(!is.na(name))
    
    # Summary plot
    title_text <- paste("Part 5: Pathogen Detection (RG3/4 Species) -", sample_name)
    
    ref_text <- paste0(
      "Pathogen detection analysis for sample ", sample_name, ":\n\n",
      "Filtering criteria:\n",
      "- Include: RG3, RG4 risk groups\n",
      "- Exclude: RG1, RG2 risk groups\n",
      "- Exclude: Species annotated as 'Microbiome' in HOMD.Category\n\n",
      "Key findings:\n",
      "- Total species in sample: ", nrow(sample_annotated), "\n",
      "- Potential pathogenic species identified: ", nrow(sample_pathogen_species), "\n"
    )
    
    if (nrow(sample_pathogen_species) > 0) {
      # Find the most abundant pathogen based on available read data
      read_col <- if ("bracken_reads" %in% colnames(sample_pathogen_species)) {
        "bracken_reads"
      } else if ("cladeReads" %in% colnames(sample_pathogen_species)) {
        "cladeReads"
      } else {
        NULL
      }
      
      if (!is.null(read_col)) {
        most_abundant_idx <- which.max(sample_pathogen_species[[read_col]])
        ref_text <- paste0(ref_text,
          "- Most abundant pathogen: ", sample_pathogen_species$name[most_abundant_idx], "\n",
          "- Average pathogen reads: ", round(mean(sample_pathogen_species[[read_col]], na.rm = TRUE), 1), "\n"
        )
      }
    }
    
    # Create summary plot using ggplot
    summary_plot <- ggplot() + 
      theme_void() +
      annotate("text", x = 0.5, y = 0.95, 
              label = paste("Part 5: Pathogen Detection (RG3/4 Species) -", sample_name),
              size = 5, hjust = 0.5, vjust = 1, fontface = "bold") +
      annotate("text", x = 0.05, y = 0.85, 
              label = ref_text,
              size = 3, hjust = 0, vjust = 1, family = "mono") +
      xlim(0, 1) + ylim(0, 1)
    
    pathogen_plots[[paste0(sample_name, "_summary")]] <- summary_plot
    
    # Generate pathogen plots if we have pathogenic species
    if (nrow(sample_pathogen_species) > 0) {
      tryCatch({
        # Source the function-based plotting (disable main execution)
        assign("..main_execution_disabled..", TRUE, envir = .GlobalEnv)
        source(file.path(config$RANALYSIS_DIR, "scripts", "Stage1_Data_acquisition", "annotate_species_plots.R"))
        
        # Risk group plot for pathogens
        pathogen_risk_plot <- generate_annotation_plot(
          df = sample_pathogen_species,
          plot_type = "risk_group",
          sample_name = paste(sample_name, "- Pathogens")
        )
        pathogen_plots[[paste0(sample_name, "_pathogen_risk_basic")]] <- pathogen_risk_plot
        
        # Check for read data for detailed plots
        read_col <- if ("bracken_reads" %in% colnames(sample_pathogen_species)) {
          "bracken_reads"
        } else if ("cladeReads" %in% colnames(sample_pathogen_species)) {
          "cladeReads"
        } else {
          NULL
        }
        
        if (!is.null(read_col)) {
          # Detailed risk group plot for pathogens
          detailed_pathogen_risk_plot <- generate_annotation_plot(
            df = sample_pathogen_species,
            plot_type = "risk_group",
            detailed = read_col,
            sample_name = paste(sample_name, "- Pathogens"),
            top_n = config$TOP_N_PLOT,
            log = TRUE
          )
          pathogen_plots[[paste0(sample_name, "_pathogen_risk_detailed")]] <- detailed_pathogen_risk_plot
          
          # Add percentage plots using calculate_species_percentages
          tryCatch({
            # Calculate species percentages for this sample's data
            sample_percentages <- calculate_species_percentages(sample_merged_long)
            
            if (nrow(sample_percentages) > 0) {
              # Filter for pathogenic species only and add percentage data
              pathogen_percentages <- sample_percentages %>%
                filter(taxID %in% sample_pathogen_species$taxID)
              
              if (nrow(pathogen_percentages) > 0) {
                # Add percentage data back to sample_pathogen_species for plotting
                sample_pathogen_with_pct <- sample_pathogen_species %>%
                  left_join(
                    pathogen_percentages %>% select(taxID, percentage),
                    by = "taxID"
                  ) %>%
                  mutate(percentage = ifelse(is.na(percentage), 0, percentage))
                
                # Percentage risk group plot for pathogens
                percentage_pathogen_risk_plot <- generate_annotation_plot(
                  df = sample_pathogen_with_pct,
                  plot_type = "risk_group",
                  detailed = "percentage",
                  sample_name = paste(sample_name, "- Pathogens"),
                  top_n = config$TOP_N_PLOT
                )
                pathogen_plots[[paste0(sample_name, "_pathogen_risk_percentage")]] <- percentage_pathogen_risk_plot
              }
            }
          }, error = function(e) {
            cat("Warning: Error creating pathogen percentage plots for", sample_name, ":", e$message, "\n")
          })
        }
        
      }, error = function(e) {
        cat("Warning: Error generating pathogen plots for", sample_name, ":", e$message, "\n")
        empty_plot <- ggplot() + 
          theme_void() + 
          labs(title = paste("Pathogen plots not available for", sample_name))
        pathogen_plots[[paste0(sample_name, "_error")]] <- empty_plot
      })
    } else {
      # No pathogenic species identified
      no_pathogen_plot <- ggplot() + 
        theme_void() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = paste("No pathogenic species identified in", sample_name), 
                size = 6, color = "darkgreen", hjust = 0.5, vjust = 0.5)
      pathogen_plots[[paste0(sample_name, "_no_pathogens")]] <- no_pathogen_plot
    }
    
    # Add comparison plot: Sample cladeReads vs Batch cladeReads_mean
    tryCatch({
      if (nrow(sample_pathogen_species) > 0 && 
          "cladeReads" %in% colnames(sample_pathogen_species) && 
          "cladeReads_mean" %in% colnames(sample_pathogen_species) &&
          any(!is.na(sample_pathogen_species$cladeReads)) &&
          any(!is.na(sample_pathogen_species$cladeReads_mean))) {
        
        # Prepare data for comparison plot
        comparison_data <- sample_pathogen_species %>%
          filter(!is.na(cladeReads) & !is.na(cladeReads_mean) & 
                 cladeReads > 0 & cladeReads_mean > 0) %>%
          select(name, cladeReads, cladeReads_mean) %>%
          slice_max(cladeReads, n = 15) %>%  # Top 15 by sample reads
          arrange(desc(cladeReads)) %>%
          mutate(name = factor(name, levels = rev(name))) %>%
          pivot_longer(cols = c(cladeReads, cladeReads_mean), 
                      names_to = "read_type", values_to = "reads") %>%
          mutate(read_type = case_when(
            read_type == "cladeReads" ~ paste("Sample", sample_name),
            read_type == "cladeReads_mean" ~ "Batch Mean",
            TRUE ~ read_type
          ))
        
        # Create named vector for scale_fill_manual
        fill_colors <- c("#1f77b4", "#ff7f0e")
        names(fill_colors) <- c(paste("Sample", sample_name), "Batch Mean")
        
        comparison_plot <- ggplot(comparison_data, aes(x = name, y = reads, fill = read_type)) +
          geom_col(position = "dodge", alpha = 0.8) +
          scale_fill_manual(values = fill_colors) +
          scale_y_log10(labels = scales::comma_format()) +
          coord_flip() +
          labs(
            title = paste("Pathogen Detection: Sample vs Batch Comparison"),
            subtitle = paste("Sample:", sample_name, "\n Top 15 pathogenic species by sample cladeReads (log scale)"),
            x = "Species",
            y = "CladeReads (log scale)",
            fill = "Data Source"
          ) +
          theme_minimal() +
          theme(
            axis.text.y = element_text(size = 8),
            legend.position = "top",
            plot.title = element_text(hjust = 0),  # Left-align title
            plot.subtitle = element_text(hjust = 0)  # Left-align subtitle
          )
        
        pathogen_plots[[paste0(sample_name, "_comparison")]] <- comparison_plot
      }
    }, error = function(e) {
      cat("Warning: Could not generate pathogen comparison plot for", sample_name, ":", e$message, "\n")
    })
    
    return(pathogen_plots)
  }
}

# Helper function to calculate species percentages by sample
calculate_species_percentages <- function(merged_long_data) {
  tryCatch({
    # Get list of species taxIDs from the input species data
    species_taxids <- unique(merged_long_data$taxID)
    
    # Determine which read column to use from merged_long_data
    if ("bracken_reads" %in% colnames(merged_long_data) && 
        any(!is.na(merged_long_data$bracken_reads))) {
      read_col <- "bracken_reads"
    } else {
      read_col <- "cladeReads"
    }
    
    # Filter merged_long_data for species and exclude human (taxID 9606) 
    species_reads <- merged_long_data %>%
      filter(taxID %in% species_taxids, taxID != "9606") %>%
      select(sample, name, taxID, all_of(read_col)) %>%
      filter(!is.na(.data[[read_col]]), .data[[read_col]] > 0)
    
    # Calculate total non-human reads per sample
    total_reads_per_sample <- merged_long_data %>%
      filter(taxID != "9606", !is.na(.data[[read_col]])) %>%
      group_by(sample) %>%
      summarise(total_reads = sum(.data[[read_col]], na.rm = TRUE), .groups = 'drop')
    
    # Calculate species percentages
    species_summary <- species_reads %>%
      left_join(total_reads_per_sample, by = "sample") %>%
      mutate(
        percentage = round((.data[[read_col]] / total_reads) * 100, 2)
      ) %>%
      filter(percentage > 0) %>%
      select(sample, name, reads = all_of(read_col), percentage) %>%
      arrange(sample, desc(percentage))
    
    return(species_summary)
    
  }, error = function(e) {
    cat("Error calculating species percentages:", e$message, "\n")
    return(data.frame())
  })
}

# Helper function to calculate pathogen percentages by sample
calculate_pathogen_percentages <- function(merged_long_data, pathogen_species_list) {
  tryCatch({
    trimmed_long <- merged_long_data %>% filter(!is.na(cladeReads))
    
    # Split by sample
    trimmed_by_sample <- split(trimmed_long, trimmed_long$sample)
    trimmed_by_sample <- lapply(trimmed_by_sample, as.data.frame)
    
    # Use calculate_species_percentages as helper function to get percentages 
    # per sample
    all_species_percentages <- do.call(rbind, lapply(trimmed_by_sample, calculate_species_percentages))
    
    if (nrow(all_species_percentages) == 0) {
      return(data.frame())
    }
    
    # If pathogen_species_list provided, filter for only those species
    if (!is.null(pathogen_species_list)) {
      pathogen_taxids <- unique(pathogen_species_list$taxID)
      pathogen_percentages <- all_species_percentages %>%
        # Get taxID from merged_long_data to filter by
        left_join(
          merged_long_data %>% select(name, taxID) %>% distinct(),
          by = "name"
        ) %>%
        filter(taxID %in% pathogen_taxids) %>%
        select(-taxID)  # Remove taxID after filtering
    } else {
      pathogen_percentages <- all_species_percentages
    }
    
    # Format for pathogen detection report Excel output
    if (nrow(pathogen_percentages) > 0) {
      pathogen_summary_formatted <- pathogen_percentages %>%
        mutate(percentage_formatted = paste0(percentage, "%")) %>%
        select(`Sample Name` = sample, 
               `RG3/4 and NotAnnotated Species` = name, 
               `Percentage of Sample Reads` = percentage_formatted) %>%
        arrange(`Sample Name`, desc(as.numeric(gsub("%", "", `Percentage of Sample Reads`))))
      
      return(pathogen_summary_formatted)
    } else {
      return(data.frame())
    }
    
  }, error = function(e) {
    cat("Error calculating pathogen percentages:", e$message, "\n")
    return(data.frame())
  })
}

#### MAIN EXECUTION #### 

# Create configuration from command line arguments
config <- create_config(args)

# Setup parallel processing
if (config$CORES > 1) {
  cl <- setup_parallel(config$CORES)
} else {
  cat("Running in single-core mode\n")
}

# Build arguments for output_processing.R
output_args <- build_output_processing_args(config)

cat("Running output_processing.R with arguments:\n")

# Run output processing
run_script_with_args(file.path(config$RANALYSIS_DIR, "scripts", "Stage1_Data_acquisition", "output_processing.R"), output_args)

if (exists("unaligned_results", envir = .GlobalEnv)) {
  unaligned_results <- get("unaligned_results", envir = .GlobalEnv)
  cat("Unaligned results loaded successfully\n")
} else {
  # Attempt to load in from files
  stop("ERROR: unaligned_results not found after running output_processing.R")
}

if (exists("nonhuman_results", envir = .GlobalEnv)) {
  nonhuman_results <- get("nonhuman_results", envir = .GlobalEnv)
  cat("Non-human results loaded successfully\n")
} else {
  cat("WARNING: nonhuman_results not found after running output_processing.R\n")
  nonhuman_results <- NULL
}

# Run correlation analysis if requested
if (config$PERFORM_CORRELATION) {
  cat("Running correlation analysis...\n")
  
  # Ensure reports directory exists
  if (!dir.exists(config$reports_dir)) {
    dir.create(config$reports_dir, recursive = TRUE)
  }
  
  if (!is.null(nonhuman_results) && exists("nonhuman_results", envir = .GlobalEnv)) {
    correlation_results <- create_correlation_analysis(
      data1 = unaligned_results$merged,
      data2 = nonhuman_results$merged,
      data1_name = "Unaligned",
      data2_name = "Non-Human",
      data1_col_suffix = "cladeReads",
      data2_col_suffix = "cladeReads",
      output_dir = config$reports_dir
    )
    
    # Check if bracken data is available and perform clade vs bracken correlation analysis
    bracken_cols_unaligned <- grep("_bracken_reads$", colnames(unaligned_results$merged), value = TRUE)
    bracken_cols_nonhuman <- if(!is.null(nonhuman_results)) grep("_bracken_reads$", colnames(nonhuman_results$merged), value = TRUE) else character(0)
    
    if (length(bracken_cols_unaligned) > 0) {
      cat("Running clade reads vs bracken reads correlation analysis for unaligned data...\n")
      bracken_correlation_unaligned <- create_correlation_analysis(
        data1 = unaligned_results$merged,
        data2 = unaligned_results$merged,
        data1_name = "Unaligned_clade_reads",
        data2_name = "Unaligned_bracken_reads",
        data1_col_suffix = "cladeReads",
        data2_col_suffix = "bracken_reads",
        output_dir = config$reports_dir
      )
    }
    
    if (length(bracken_cols_nonhuman) > 0) {
      cat("Running clade reads vs bracken reads correlation analysis for non-human data...\n")
      bracken_correlation_nonhuman <- create_correlation_analysis(
        data1 = nonhuman_results$merged,
        data2 = nonhuman_results$merged,
        data1_name = "Non-Human_clade_reads",
        data2_name = "Non-Human_bracken_reads",
        data1_col_suffix = "cladeReads",
        data2_col_suffix = "bracken_reads",
        output_dir = config$reports_dir
      )
    }
    
    # Analysis of host contamination
    cat("Analyzing host contamination...\n")
    host_contamination <- setdiff(unaligned_results$merged$name, nonhuman_results$merged$name)
    cat("Species found only in unaligned data (potential host contamination):", length(host_contamination), "\n")
    
    if (length(host_contamination) > 0) {
      host_contaminated_results <- unaligned_results$merged %>%
        filter(name %in% host_contamination) %>%
        select(name, contains("cladeReads")) %>%
        arrange(desc(rowSums(select(., contains("cladeReads")), na.rm = TRUE)))
      
      cat("Host-contaminated species:\n")
      print(host_contaminated_results$name)
      
      # Save host contamination results
      write.csv(host_contaminated_results, 
                file.path(config$reports_dir, "host_contamination_species.csv"), 
                row.names = FALSE)
    }
  } else {
    cat("Skipping correlation analysis - nonhuman results not available\n")
  }
}

#### Batch Report Exporting ####

# Run read statistics report
if (config$RR_STATS) {
  tryCatch({
    if (exists("rrstats_data", envir = .GlobalEnv)) {
      rrstats_data <- get("rrstats_data", envir = .GlobalEnv)
      cat("Read statistics data loaded successfully\n")
    } else {
      # Load read statistics data from the latest file
      cat("Loading read statistics data from the latest generated file...\n")
      rrstats_data <- read.csv(get_latest_timestamped_file(pattern = "read_statistics", path = config$output_dir), 
                               stringsAsFactors = FALSE)
    }
    runread_stats_report(rrstats_data, config$reports_dir)
  }, error = function(e) {
    cat("Warning: Read statistics report skipped - no read statistics data available:", e$message, "\n")
  })
}

# Run kraken statistics report
if (config$KR_STATS) {
  tryCatch({
    if (exists("combined_report_data", envir = .GlobalEnv)) {
      combined_report_data <- get("combined_report_data", envir = .GlobalEnv)
      cat("Combined report data loaded successfully\n")
    } else {
      combined_report_data <- read.csv(get_latest_timestamped_file(pattern = "sample_report_data", path = config$output_dir), 
                                       stringsAsFactors = FALSE)
    }
    
    kraken_stats_report(unaligned_results, nonhuman_results, combined_report_data, config, config$reports_dir, TRUE)
  }, error = function(e) {
    cat("Warning: Kraken statistics report skipped - no kraken data available:", e$message, "\n")
  })
}

species_df <- unaligned_results$species_list
assign("species_df", species_df, envir = .GlobalEnv)
# Run annotation with custom databases directory
run_script_with_args(file.path(config$RANALYSIS_DIR, "scripts", "Stage1_Data_acquisition", "annotate_species.R"), 
                     c("--df", "species_df", "--output", "annotated_species", 
                       "--databases-dir", config$DATABASES_DIR))

# Run species annotation if annotated species_list is available
if (config$SA_STATS) {
  cat("Running species annotation...\n")

  if (exists("annotated_species", envir = .GlobalEnv)) {
    annotated_species <- get("annotated_species", envir = .GlobalEnv)
    merged_long_data <- if (!is.null(nonhuman_results) && !is.null(nonhuman_results$merged_long)) {
      nonhuman_results$merged_long
    } else {
      unaligned_results$merged_long
    }
    
    # Merge species annotation columns to the merged_long_data joining by taxID
    full_results_with_annotations <- merged_long_data %>%
      left_join(annotated_species %>% select(taxID, RiskGroup, HOMD, HOMD.Category), by = "taxID") %>%
      filter(!is.na(cladeReads))  # Trim rows where cladeReads are NA
    
    # Save the full results with annotations using output_csv_file
    output_csv_file(full_results_with_annotations, "full_results_with_annotated_species", config$reports_dir, "rra")
    cat("Full results with species annotations saved to output directory\n")
    
    species_annotation_report(annotated_species, config$reports_dir, TRUE, config)
    cat("Species annotation completed successfully\n")

  } else {
    stop("ERROR: annotated_species not found after running annotate_species.R")
  }
}

# Run pathogen detection if requested & annotated species_list is available
if (config$PD_STATS) {
  if (exists("annotated_species", envir = .GlobalEnv)) {
    annotated_species <- get("annotated_species", envir = .GlobalEnv)
    cat("Running pathogen detection analysis...\n")
    # Get merged_long data for sample-level calculations
    merged_long_data <- if (!is.null(nonhuman_results) && !is.null(nonhuman_results$merged_long)) {
      nonhuman_results$merged_long
    } else {
      unaligned_results$merged_long
    }
    # pathogen_detection_report(annotated_species, merged_long_data, config$reports_dir, TRUE, config)
    pathogen_species <- pathogen_detection_report(annotated_species, merged_long_data, config$reports_dir, TRUE, config)
  } else {
    stop("ERROR: annotated_species not found after running annotate_species.R")
  }
}

# Clean up parallel cluster
if (exists("cl")) {
  stopCluster(cl)
  rm(cl)
  cat("Parallel cluster stopped\n")
}

cat("Analysis completed successfully!\n")
cat("Output directory:", config$output_dir, "\n")