#!/usr/bin/env Rscript

#### Header ####
# Output Processing Script
# Author: Denis Rivard, modified by GitHub Copilot
# Description: Improved and simplified kraken/bracken output processing script
# Compatible with various kreport file naming patterns

library(here)
library(purrr)
library(tidyr)
library(dplyr)
library(stringr)

# Utility functions from utils.R
output_csv_file <- function(x, file_name, dir, script_name) {
  timestamp <- format(Sys.time(), "%y%m%d")
  # Construct the full file path
  file_path <- file.path(here(dir), paste0(file_name, timestamp, script_name, ".csv"))
  write.csv(x = x, file = file_path, row.names = F)
}

# The following 4 functions were adapted from
# https://github.com/fbreitwieser/pavian/blob/a910755648d6f3d7256fac73315bc837b75e8663/R/datainput-read_report.R
# Florian P Breitwieser, Steven L Salzberg, Pavian: interactive analysis of metagenomics data for microbiome studies and pathogen identification,
# Bioinformatics, Volume 36, Issue 4, February 2020, Pages 1303–1304, https://doi.org/10.1093/bioinformatics/btz715
## recursive function that transforms the kraken dataframe into a cascading list
build_kraken_tree <- function(report) {
  if (nrow(report) == 0 || nrow(report) == 1) {
    # this should only happen if the original input to the function has a size <= 1
    return(list(report))
  }
  
  ## select the current depth as the one of the topmost data.frame row
  sel_depth <- report[,'depth'] == report[1,'depth']
  
  ## partition the data.frame into parts with that depth
  depth_partitions <- cumsum(sel_depth)
  
  ## for each depth partition
  res <- lapply(unique(depth_partitions),
                function(my_depth_partition) {
                  sel <- depth_partitions == my_depth_partition
                  
                  ## return the data.frame row if it is only one row (leaf node, ends recursion)
                  if (sum(sel) == 1)
                    return(report[sel,,drop=F])
                  
                  ## otherwise: take first row as partition descriptor ..
                  first_row <- which(sel)[1]
                  ##  and recurse deeper into the depths with the remaining rows
                  dres <- build_kraken_tree(report[which(sel)[-1],,drop=F])
                  
                  attr(dres,"row") <- report[first_row,,drop=F]
                  dres
                })
  names(res) <- report$name[sel_depth]
  res
}

## Collapse taxonomic taxRanks to only those mentioned in keep_taxRanks
collapse.taxRanks <- function(krakenlist,keep_taxRanks=LETTERS,filter_taxon=NULL) {
  ## input: a list, in which each element is either a
  ##            a list or a data.frame (for the leafs)
  ##   the input has an attribute row that gives details on the current taxRank
  
  ## columns whose values are added to the next taxRank when
  ##  a taxRank is deleted
  cols <- c("taxonReads","n_unique_kmers","n_kmers")
  if (length(krakenlist) == 0 || is.data.frame(krakenlist)) {
    return(krakenlist)
  }
  
  parent_row <- attr(krakenlist,"row")
  all.child_rows <- c()
  
  if (is.null(parent_row)) {
    return(do.call(rbind,lapply(krakenlist,collapse.taxRanks,keep_taxRanks=keep_taxRanks,filter_taxon=filter_taxon)))
  }
  
  ## rm.cladeReads captures the number of cladeReads that are deleted.
  ##  this has to be propagated to the higher taxRank
  rm.cladeReads <- 0
  
  for (kl in krakenlist) {
    if (is.data.frame(kl)) {  ## is a leaf node?
      child_rows <- kl
    } else {                 ## recurse deeper into tree
      child_rows <- collapse.taxRanks(kl,keep_taxRanks,filter_taxon=filter_taxon)
      if ('rm.cladeReads' %in% names(attributes(child_rows))) {
        rm.cladeReads <- rm.cladeReads + attr(child_rows,'rm.cladeReads')
      }
    }
    
    ## check if this taxRank and the taxRanks below should be removed
    delete.taxon <- child_rows[1,'name'] %in% filter_taxon
    if (delete.taxon) {
      rm.cladeReads <- rm.cladeReads + child_rows[1,'cladeReads']
      # dmessage(sprintf("removed %7s cladeReads, including %s childs, for %s",child_rows[1,'"cladeReads"'],nrow(child_rows)-1,child_rows[1,'name']))
      
      ## remove all children
      child_rows <- NULL
      
    } else {
      
      ## check if the last (top-most) row should be kept
      keep_last.child <- child_rows[1,'taxRank'] %in% keep_taxRanks
      
      if (!keep_last.child) {
        cols <- cols[cols %in% colnames(parent_row)]
        
        ## save the specified colum information to the parent
        parent_row[,cols] <- parent_row[,cols] + child_rows[1,cols]
        
        ## remove row
        child_rows <- child_rows[-1,,drop=FALSE]
        
        ## decrease depths of rows below child row
        if (nrow(child_rows) > 0)
          child_rows[,'depth'] <- child_rows[,'depth'] - 1
        
      }
    }
    all.child_rows <- rbind(all.child_rows,child_rows)
  }
  
  ## subtract deleted read count from parent row
  parent_row[,'cladeReads'] <- parent_row[,'cladeReads'] - rm.cladeReads
  res <- rbind(parent_row,all.child_rows)
  
  if (parent_row[,'cladeReads'] < 0)
    stop("mistake made in removing cladeReads")
  #if (parent_row[,'"cladeReads"'] == 0)
  #  res <- c()
  
  if (rm.cladeReads > 0)
    attr(res,'rm.cladeReads') <- rm.cladeReads
  return(res)
}

delete_taxRanks_below <- function(report,taxRank="S") {
  del_taxRank <- 0
  do_del <- FALSE
  del_row <- 0
  
  cols <- c("taxonReads","n_unique_kmers","n_kmers")
  sub.sums <- c(0,0,0)
  
  rows_to_delete <- c()
  for (i in seq_len(nrow(report))) {
    if (report[i,'taxRank'] %in% taxRank) {
      del_depth <- report[i,'depth']
      do_del <- TRUE
      del_row <- i
      sub.sums <- c(0,0,0)
    } else {
      if (do_del) {
        if (report[i,'depth'] > del_taxRank) {
          rows_to_delete <- c(rows_to_delete,i)
          sub.sums <- sub.sums + report[i,cols]
        } else {
          report[del_row,cols] <- report[del_row,cols]+sub.sums
          sub.sums <- c(0,0,0)
          do_del <- FALSE
        }
      }
    }
  }
  report[-rows_to_delete,]
}

#' Read kraken or centrifuge-style report
#'
#' @param myfile kraken report file
#' @param collapse  should the results be collapsed to only those taxRanks specified in keep_taxRanks?
#' @param keep_taxRanks taxRanks to keep when collapse is TRUE
#' @param min.depth minimum depth
#' @param filter_taxon filter certain taxon names
#' @param has_header if the kraken report has a header or not
#' @param add_taxRank_columns if TRUE, for each taxRank columns are added
#'
#' @return report data.frame
#' @export
#'
read_report2 <- function(myfile,collapse=TRUE,keep_taxRanks=c("D","K","P","C","O","F","G","S"),min.depth=0,filter_taxon=NULL,
                         has_header=NULL,add_taxRank_columns=FALSE) {
  first.line <- readLines(myfile,n=1)
  isASCII <-  function(txt) all(charToRaw(txt) <= as.raw(127))
  if (!isASCII(first.line)) {
    # dmessage(myfile," is no valid report - not all characters are ASCII")
    return(NULL)
  }
  if (is.null(has_header)) {
    has_header <- grepl("^[a-zA-Z]",first.line)
  }
  
  if (has_header) {
    report <- utils::read.table(myfile,sep="\t",header = T,
                                quote = "",stringsAsFactors=FALSE, comment.char="#")
    #colnames(report) <- c("percentage","cladeReads","taxonReads","taxRank","taxID","n_unique_kmers","n_kmers","perc_uniq_kmers","name")
    
    ## harmonize column names. TODO: Harmonize them in the scripts!
    colnames(report)[colnames(report)=="clade_perc"] <- "percentage"
    colnames(report)[colnames(report)=="perc"] <- "percentage"
    
    colnames(report)[colnames(report)=="n_reads_clade"] <- "cladeReads"
    colnames(report)[colnames(report)=="n.clade"] <- "cladeReads"
    
    colnames(report)[colnames(report)=="n_reads_taxo"] <- "taxonReads"
    colnames(report)[colnames(report)=="n.stay"] <- "taxonReads"
    
    colnames(report)[colnames(report)=="rank"] <- "taxRank"
    colnames(report)[colnames(report)=="tax_rank"] <- "taxRank"
    
    colnames(report)[colnames(report)=="taxonid"] <- "taxID"
    colnames(report)[colnames(report)=="tax"] <- "taxID"
    
  } else {
    report <- utils::read.table(myfile,sep="\t",header = F,
                                col.names = c("percentage","cladeReads","taxonReads","minimizers", "distinct_minimizers","taxRank","taxID","name"),
                                quote = "",stringsAsFactors=FALSE, comment.char="#")
  }
  
  report$depth <- nchar(gsub("\\S.*","",report$name))/2
  report$name <- gsub("^ *","",report$name)
  report$name <- paste(tolower(report$taxRank),report$name,sep="_")
  
  ## Only stop at certain taxRanks
  ## filter taxon and further up the tree if 'filter_taxon' is defined
  kraken.tree <- build_kraken_tree(report)
  report <- collapse.taxRanks(kraken.tree,keep_taxRanks=keep_taxRanks,filter_taxon=filter_taxon)
  
  ## Add a metaphlan-style taxon string
  if (add_taxRank_columns) {
    report[,keep_taxRanks] <- NA
  }
  report$taxLineage = report$name
  rows_to_consider <- rep(FALSE,nrow(report))
  
  for (i in seq_len(nrow(report))) {
    ## depth > 2 correspond to taxRanks below 'D'
    if (i > 1 && report[i,"depth"] > min.depth) {
      ## find the maximal index of a row below the current depth
      idx <- report$depth < report[i,"depth"] & rows_to_consider
      if (!any(idx)) { next() }
      
      current.taxRank <- report[i,'taxRank']
      my_row <- max(which(idx))
      report[i,'taxLineage'] <- paste(report[my_row,'taxLineage'],report[i,'taxLineage'],sep="|")
      
      if (add_taxRank_columns) {
        if (report[my_row,'taxRank'] %in% keep_taxRanks) {
          taxRanks.cp <- keep_taxRanks[seq(from=1,to=which(keep_taxRanks == report[my_row,'taxRank']))]
          report[i,taxRanks.cp] <- report[my_row,taxRanks.cp]
        }
        
        report[i,report[i,'taxRank']] <- report[i,'name']
      }
    }
    rows_to_consider[i] <- TRUE
  }
  
  report <- report[report$depth >= min.depth,]
  
  report$percentage <- round(report$cladeReads/sum(report$taxonReads),6) * 100
  
  for (column in c("taxonReads", "cladeReads")) 
    if (all(floor(report[[column]]) == report[[column]])) 
      report[[column]] <- as.integer(report[[column]])
  
  if ('n_unique_kmers'  %in% colnames(report))
    report$kmerpercentage <- round(report$n_unique_kmers/sum(report$n_unique_kmers,na.rm=T),6) * 100
  #report$taxRankperc <- 100/taxRank(report$cladeReads)
  
  rownames(report) <- NULL
  
  report
}

# Load parallel libraries conditionally
if ("--parallel" %in% commandArgs(trailingOnly = TRUE)) {
  library(foreach)
  library(doParallel)
}

# Function to display help
show_help <- function() {
  cat("Usage: Rscript output_processing.R [OPTIONS]\n\n")
  cat("Inputs:\n")
  cat("  --unaligned-kreports <DIR> Path to unaligned kreports folder (default: NULL)\n")
  cat("  --unaligned-bracken <DIR>  Path to unaligned bracken folder (default: NULL)\n")
  cat("  --nonhuman-kreports <DIR>  Path to nonhuman kreports folder (default: NULL)\n")
  cat("  --nonhuman-bracken <DIR>   Path to nonhuman bracken folder (default: NULL)\n")
  cat("  --rpm-dir <DIR>           Path to RPM count matrix folder (default: data/METASEQ_pipeline/RPM)\n")
  cat("  --counts-output-dir <DIR> Path to output filtered counts (default: data/processed)\n")
  cat("  --rrna-file <FILE>        Path to rRNA gene list CSV (default: Ranalysis/databases/hgnc_rRNAgenelist.csv)\n")
  cat("  --runtime-dir <DIR>       Path to runtime folder for parsing runtime files (default: NULL)\n")
  cat("  --rrstats-dir <DIR>       Path to read statistics folder for parsing rrstats files (default: NULL)\n")
  cat("  --metadata-dir <DIR>      Path to directory containing sample metadata CSV files (default: NULL)\n")
  cat("Options:\n")
  cat("  (bracken processing is enabled only when a bracken directory is provided)\n")
  cat("  --subspecies              Include subspecies (S1, S2, S3) in addition to species (default: FALSE)\n")
  cat("  --params                  Add confidence_levels, minimum_hit_groups, and human_reads columns to long data (default: FALSE)\n")
  cat("  --parallel                Use parallel processing with cluster from global environment 'cl' (default: FALSE)\n")
  cat("  --skip-counts             Skip RPM count matrix merging and gene filtering (default: FALSE)\n")
  cat("  --output                  Save results to CSV files (default: FALSE)\n")
  cat("  --output-dir <DIR>         Output directory (default: outputs/full_run)\n")
  cat("Filtering:\n")
  cat("  --exclude-taxid <ID>       Exclude species with this taxonomy ID (e.g., 9606 for Homo sapiens)\n")
  cat("  --min-reads <N>           Minimum clade reads threshold (default: 0)\n")
  cat("  --top-n <N>               Number of top species to analyze (default: 3)\n")
  cat("  --minimizer-ratio <R>     Filter by minimum ratio of distinct_minimizers/cladeReads (default: NULL)\n")
  cat("  --minimizer-threshold <N> Filter by minimum distinct_minimizers threshold (default: NULL)\n")
  cat("  --no-protein-coding       Skip protein-coding-only gene filtering (default: FALSE)\n")
  cat("  --no-remove-mt            Skip mitochondrial gene removal (default: FALSE)\n")
  cat("  --no-go-annotation        Skip GO annotation gene filtering (default: FALSE)\n")
  cat("Examples:\n")
  cat("  Rscript output_processing_v5.R                          # Default settings\n")
  cat("  Rscript output_processing_v5.R --exclude-taxid 9606     # Exclude Homo sapiens\n")
  cat("  Rscript output_processing_v5.R --unaligned-kreports /path/to/kreports # Custom kreports path\n")
  cat("  Rscript output_processing_v5.R --output-dir /path/to/output # Custom output directory\n")
  cat("  Rscript output_processing_v5.R --output                 # Save results to CSV files\n")
  cat("  Rscript output_processing_v5.R --params                 # Add parameter columns to long data\n")
  cat("  Rscript output_processing_v5.Rs --subspecies             # Include subspecies in analysis\n")
  cat("  Rscript output_processing_v5.R --spec-list              # Output complete species list\n")
  cat("  Rscript output_processing_v5.R --minimizer-ratio 0.1    # Filter by minimizer ratio >= 0.1\n")
  cat("  Rscript output_processing_v5.R --minimizer-threshold 5  # Filter by distinct_minimizers >= 5\n")
  cat("  Rscript output_processing_v5.R --parallel               # Use parallel processing\n")
  cat("  Rscript output_processing_v5.R --runtime-dir /path/to/runtimes # Parse runtime files\n")
  cat("  Rscript output_processing_v5.R --rrstats-dir /path/to/rrstats # Parse read statistics files\n")
  cat("  Rscript output_processing_v5.R --metadata-dir /path/to/metadata # Parse sample metadata files\n\n")
}
# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check for help request
if ("--help" %in% args || "-h" %in% args) {
  show_help()
  stop("Help requested. Execution stopped.")
}

options(scipen = 999)

# Configuration variables with default values
MIN_CLADE_READS <- 0
EXCLUDE_TAXID <- NULL  # Set to NULL by default (no exclusion)
TOP_SPECIES_N <- 10
PROCESS_BRACKEN <- TRUE   # Set to TRUE to enable bracken processing by default
SAVE_OUTPUT <- FALSE  # Set to FALSE to skip saving CSV files
ADD_PARAMS <- FALSE  # Set to FALSE to skip adding parameter columns
INCLUDE_SUBSPECIES <- FALSE  # Set to FALSE to exclude subspecies
MINIMIZER_RATIO <- NULL  # Set to NULL to skip minimizer ratio filtering
MINIMIZER_THRESHOLD <- NULL  # Set to NULL to skip minimizer threshold filtering
USE_PARALLEL <- FALSE  # Set to FALSE to skip parallel processing
PROCESS_COUNTS <- TRUE  # Set to TRUE to merge RPM count matrices and filter genes
PROTEIN_CODING_FILTER <- TRUE  # Set to TRUE to keep only protein-coding genes
REMOVE_MT_GENES <- TRUE  # Set to TRUE to remove MT- genes
GO_ANNOTATION_FILTER <- TRUE  # Set to TRUE to keep genes with any GO annotation

# Default input directory paths
UNALIGNED_KREPORTS_DIR <- NULL
UNALIGNED_BRACKEN_DIR <- NULL
NONHUMAN_KREPORTS_DIR <- NULL
NONHUMAN_BRACKEN_DIR <- NULL
RUNTIME_DIR <- NULL
RRSTATS_DIR <- NULL
METADATA_DIR <- NULL
OUTPUT_DIR <- here("data", "processed")
RPM_DIR <- here("data", "METASEQ_pipeline", "RPM")
COUNTS_OUTPUT_DIR <- here("data", "processed")
RRNA_FILE <- here("Ranalysis", "databases", "hgnc_rRNAgenelist.csv")

# Parse command line arguments for configuration
if (length(args) > 0) {
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--exclude-taxid" && i < length(args)) {
      EXCLUDE_TAXID <- as.numeric(args[i + 1])
      cat("Excluding taxID:", EXCLUDE_TAXID, "\n")
      i <- i + 2
    } else if (args[i] == "--min-reads" && i < length(args)) {
      MIN_CLADE_READS <- as.numeric(args[i + 1])
      cat("Minimum clade reads set to:", MIN_CLADE_READS, "\n")
      i <- i + 2
    } else if (args[i] == "--top-n" && i < length(args)) {
      TOP_SPECIES_N <- as.numeric(args[i + 1])
      cat("Top N species set to:", TOP_SPECIES_N, "\n")
      i <- i + 2
    } else if (args[i] == "--unaligned-kreports" && i < length(args)) {
      UNALIGNED_KREPORTS_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Unaligned kreports directory set to:", UNALIGNED_KREPORTS_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--unaligned-bracken" && i < length(args)) {
      UNALIGNED_BRACKEN_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Unaligned bracken directory set to:", UNALIGNED_BRACKEN_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--nonhuman-kreports" && i < length(args)) {
      NONHUMAN_KREPORTS_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Nonhuman kreports directory set to:", NONHUMAN_KREPORTS_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--nonhuman-bracken" && i < length(args)) {
      NONHUMAN_BRACKEN_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Nonhuman bracken directory set to:", NONHUMAN_BRACKEN_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--runtime-dir" && i < length(args)) {
      RUNTIME_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Runtime directory set to:", RUNTIME_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--rrstats-dir" && i < length(args)) {
      RRSTATS_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Read statistics directory set to:", RRSTATS_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--metadata-dir" && i < length(args)) {
      METADATA_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Sample metadata directory set to:", METADATA_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--rpm-dir" && i < length(args)) {
      RPM_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("RPM directory set to:", RPM_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--counts-output-dir" && i < length(args)) {
      COUNTS_OUTPUT_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Counts output directory set to:", COUNTS_OUTPUT_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--rrna-file" && i < length(args)) {
      RRNA_FILE <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("rRNA gene list set to:", RRNA_FILE, "\n")
      i <- i + 2
    } else if (args[i] == "--output-dir" && i < length(args)) {
      OUTPUT_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Output directory set to:", OUTPUT_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--output") {
      SAVE_OUTPUT <- TRUE
      cat("CSV output enabled\n")
      i <- i + 1
    } else if (args[i] == "--params") {
      ADD_PARAMS <- TRUE
      cat("Parameter columns will be added to long data\n")
      i <- i + 1
    } else if (args[i] == "--subspecies") {
      INCLUDE_SUBSPECIES <- TRUE
      cat("Subspecies will be included in analysis\n")
      i <- i + 1
    } else if (args[i] == "--minimizer-ratio" && i < length(args)) {
      MINIMIZER_RATIO <- as.numeric(args[i + 1])
      cat("Minimizer ratio filtering set to:", MINIMIZER_RATIO, "\n")
      i <- i + 2
    } else if (args[i] == "--minimizer-threshold" && i < length(args)) {
      MINIMIZER_THRESHOLD <- as.numeric(args[i + 1])
      cat("Minimizer threshold filtering set to:", MINIMIZER_THRESHOLD, "\n")
      i <- i + 2
    } else if (args[i] == "--parallel") {
      USE_PARALLEL <- TRUE
      cat("Parallel processing enabled\n")
      i <- i + 1
    } else if (args[i] == "--skip-counts") {
      PROCESS_COUNTS <- FALSE
      cat("RPM count matrix processing disabled\n")
      i <- i + 1
    } else if (args[i] == "--no-protein-coding") {
      PROTEIN_CODING_FILTER <- FALSE
      cat("Protein-coding gene filtering disabled\n")
      i <- i + 1
    } else if (args[i] == "--no-remove-mt") {
      REMOVE_MT_GENES <- FALSE
      cat("Mitochondrial gene removal disabled\n")
      i <- i + 1
    } else if (args[i] == "--no-go-annotation") {
      GO_ANNOTATION_FILTER <- FALSE
      cat("GO annotation gene filtering disabled\n")
      i <- i + 1
    } else {
      cat("Unknown argument:", args[i], "\n")
      show_help()
      stop("Help requested. Execution stopped.")
    }
  }
}
rm(i, show_help)

# Decide whether to process bracken: enable only if at least one bracken directory was provided
if (is.null(UNALIGNED_BRACKEN_DIR) && is.null(NONHUMAN_BRACKEN_DIR)) {
  PROCESS_BRACKEN <- FALSE
  cat("No bracken directories provided; skipping bracken processing\n")
} else {
  PROCESS_BRACKEN <- TRUE
}

# Display final configuration
cat("\n=== Configuration ===\n")
cat("Minimum clade reads:", MIN_CLADE_READS, "\n")
cat("Top N species:", TOP_SPECIES_N, "\n")
cat("Process bracken:", PROCESS_BRACKEN, "\n")
cat("Save CSV output:", SAVE_OUTPUT, "\n")
cat("Add parameter columns:", ADD_PARAMS, "\n")
cat("Include subspecies:", INCLUDE_SUBSPECIES, "\n")
cat("Minimizer ratio threshold:", if(is.null(MINIMIZER_RATIO)) "None" else MINIMIZER_RATIO, "\n")
cat("Minimizer absolute threshold:", if(is.null(MINIMIZER_THRESHOLD)) "None" else MINIMIZER_THRESHOLD, "\n")
cat("Use parallel processing:", USE_PARALLEL, "\n")
if (is.null(EXCLUDE_TAXID)) {
  cat("Exclude taxID: None (including all species)\n")
} else {
  cat("Exclude taxID:", EXCLUDE_TAXID, "\n")
}
cat("\nInput directories:\n")
cat("  Unaligned kreports:", UNALIGNED_KREPORTS_DIR, "\n")
cat("  Unaligned bracken:", UNALIGNED_BRACKEN_DIR, "\n")
cat("  Nonhuman kreports:", NONHUMAN_KREPORTS_DIR, "\n")
cat("  Nonhuman bracken:", NONHUMAN_BRACKEN_DIR, "\n")
cat("  Runtime directory:", RUNTIME_DIR, "\n")
cat("  Read statistics directory:", RRSTATS_DIR, "\n")
cat("  Sample metadata directory:", METADATA_DIR, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n")
cat("  RPM directory:", RPM_DIR, "\n")
cat("  Counts output directory:", COUNTS_OUTPUT_DIR, "\n")
cat("  rRNA gene list:", RRNA_FILE, "\n")
cat("  Process RPM count matrices:", PROCESS_COUNTS, "\n")
cat("  Protein-coding filter:", PROTEIN_CODING_FILTER, "\n")
cat("  Remove MT genes:", REMOVE_MT_GENES, "\n")
cat("  GO annotation filter:", GO_ANNOTATION_FILTER, "\n")
cat("=====================\n\n")

# Setup parallel processing if requested
if (USE_PARALLEL) {
  if (exists("cl", envir = .GlobalEnv)) {
    cat("Setting up parallel processing with cluster from global environment\n")
    registerDoParallel(cl)
    cat("Number of workers:", getDoParWorkers(), "\n")
  } else {
    warning("Parallel processing requested but no cluster 'cl' found in global environment")
    cat("Setting up parallel processing with new cluster \n")
    cl_setup_start_time <- Sys.time()
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    cat("Cluster setup completed in", round(difftime(Sys.time(), cl_setup_start_time, units = "secs"), 2), "seconds\n")
    rm(cl_setup_start_time)
    # cat("To create a cluster, run: cl <- parallel::makeCluster(8) # or desired number of cores\n")
    # USE_PARALLEL <- FALSE
    # cat("Falling back to sequential processing\n")
  }
}

#### Logic ####
# Enhanced file parsing function that handles multiple naming patterns
parse_sample_info <- function(filename) {
  # Remove file extension
  sample_name <- tools::file_path_sans_ext(basename(filename))
  
  # Initialize default values
  result <- list(
    original_name = sample_name,
    trimmed_name = sample_name,
    confidence = NA,
    minimum_hits = NA,
    database_used = "unknown"
  )
  
  # Pattern 1: Standard pattern like "W50504632_S109_c45m30stdk2dbdl_report"
  if (grepl("_S\\d+_c\\d+m\\d+", sample_name)) {
    result$trimmed_name <- sub("(_S\\d+).*", "", sample_name)
    conf_min_db <- regmatches(sample_name, regexpr("c\\d+m\\d+[A-Za-z0-9_]*", sample_name))
    if (length(conf_min_db) > 0) {
      result$confidence <- as.integer(sub("c(\\d+).*", "\\1", conf_min_db))
      result$minimum_hits <- as.integer(sub(".*m(\\d+).*", "\\1", conf_min_db))
      db_part <- sub(".*m\\d+([A-Za-z0-9_]+).*", "\\1", conf_min_db)
      result$database_used <- sub("_report$", "", db_part)  # Remove trailing _report
    }
  }
  # Pattern 2: Pattern like "temp_ZQ3002179_B_c0m45stdk2dbdl_report" or "ZQ1400570_A_Pos_c45m30stdk2dbdl_report"
  else if (grepl("_[A-Z](_[A-Za-z]*)?_c\\d+m\\d+", sample_name)) {
    result$trimmed_name <- sub("_c\\d+.*", "", sample_name)
    conf_min_db <- regmatches(sample_name, regexpr("c\\d+m\\d+[A-Za-z0-9_]*", sample_name))
    if (length(conf_min_db) > 0) {
      result$confidence <- as.integer(sub("c(\\d+).*", "\\1", conf_min_db))
      result$minimum_hits <- as.integer(sub(".*m(\\d+).*", "\\1", conf_min_db))
      db_part <- sub(".*m\\d+([A-Za-z0-9_]+).*", "\\1", conf_min_db)
      result$database_used <- sub("_report$", "", db_part)  # Remove trailing _report
    }
  }
  # Pattern 3: General pattern that captures any filename with c#m# format
  else if (grepl("c\\d+m\\d+", sample_name)) {
    result$trimmed_name <- sub("_c\\d+.*", "", sample_name)
    conf_min_db <- regmatches(sample_name, regexpr("c\\d+m\\d+[A-Za-z0-9_]*", sample_name))
    if (length(conf_min_db) > 0) {
      result$confidence <- as.integer(sub("c(\\d+).*", "\\1", conf_min_db))
      result$minimum_hits <- as.integer(sub(".*m(\\d+).*", "\\1", conf_min_db))
      db_part <- sub(".*m\\d+([A-Za-z0-9_]+).*", "\\1", conf_min_db)
      result$database_used <- sub("_report$", "", db_part)  # Remove trailing _report
    }
  }
  # Pattern 4: Specific pattern for bracken files with S# like "W50504632_S109_bracken_S" or "W50504632_S109_bracken_S1"
  else if (grepl("_S\\d+_bracken_S\\d*$", sample_name)) {
    result$trimmed_name <- sub("_S\\d+_bracken_S\\d*$", "", sample_name)
    result$confidence <- NA
    result$minimum_hits <- NA
    result$database_used <- "unknown"
  }
  # Pattern 5: Specific pattern for bracken files without S# like "W50504632_bracken_S" or "W50504632_bracken_S1"
  else if (grepl("_bracken_S\\d*$", sample_name)) {
    result$trimmed_name <- sub("_bracken_S\\d*$", "", sample_name)
    result$confidence <- NA
    result$minimum_hits <- NA
    result$database_used <- "unknown"
  }
  # Pattern 6: Simple pattern with just confidence (fallback)
  else if (grepl("c\\d+", sample_name)) {
    result$trimmed_name <- sub("_c\\d+.*", "", sample_name)
    conf_match <- regmatches(sample_name, regexpr("c\\d+", sample_name))
    if (length(conf_match) > 0) {
      result$confidence <- as.integer(sub("c(\\d+)", "\\1", conf_match))
    }
  }
  # Pattern 7: No specific pattern - use full name
  else {
    result$trimmed_name <- sample_name
  }
  return(result)
} # Denis Rivard: LGTM 30/7/2025

# Function to filter kreport data based on clade reads and minimizer criteria
filter_kreport_data <- function(kreport_df, sample_name, min_clade_reads = MIN_CLADE_READS, 
                                minimizer_ratio = MINIMIZER_RATIO, minimizer_threshold = MINIMIZER_THRESHOLD,
                                exclude_taxid = EXCLUDE_TAXID) {
  initial_count <- nrow(kreport_df)
  
  # Apply taxID exclusion filter first if specified
  if (!is.null(exclude_taxid)) {
    excluded_count <- kreport_df %>% filter(taxID == exclude_taxid) %>% nrow()
    kreport_df <- kreport_df %>%
      filter(taxID != exclude_taxid)
    if (excluded_count > 0) {
      cat("Sample", sample_name, ": excluded", excluded_count, "entries with taxID", exclude_taxid, "\n")
    }
  }
  
  # Apply minimum clade reads filter
  if (min_clade_reads > 0) {
    pre_reads_count <- nrow(kreport_df)
    kreport_df <- kreport_df %>%
      filter(cladeReads >= min_clade_reads)
    cat("Sample", sample_name, ": filtered from", pre_reads_count, "to", nrow(kreport_df), 
        "species with >=", min_clade_reads, "clade reads\n")
  }
  
  # Apply minimizer filtering if either parameter is specified
  if (!is.null(minimizer_ratio) || !is.null(minimizer_threshold)) {
    pre_minimizer_count <- nrow(kreport_df)
    
    # Apply ratio filter if specified
    if (!is.null(minimizer_ratio)) {
      kreport_df <- kreport_df %>%
        filter(is.na(distinct_minimizers) | distinct_minimizers == 0 | 
                 (cladeReads / pmax(distinct_minimizers, 1)) >= (1/minimizer_ratio))
      cat("Sample", sample_name, ": applied minimizer ratio filter (>=", minimizer_ratio, 
          "): from", pre_minimizer_count, "to", nrow(kreport_df), "species\n")
    }
    
    # Apply threshold filter if specified
    if (!is.null(minimizer_threshold)) {
      threshold_pre_count <- nrow(kreport_df)
      kreport_df <- kreport_df %>%
        filter(is.na(distinct_minimizers) | distinct_minimizers >= minimizer_threshold)
      cat("Sample", sample_name, ": applied minimizer threshold filter (>=", minimizer_threshold, 
          "): from", threshold_pre_count, "to", nrow(kreport_df), "species\n")
    }
  }
  
  return(kreport_df)
}

# Simplified function to process kreport files
process_kreport_folder <- function(kreports_folder, include_subspecies = INCLUDE_SUBSPECIES) {
  if (!dir.exists(kreports_folder)) {
    stop("Directory does not exist: ", kreports_folder)
  }
  
  kreports_files <- list.files(kreports_folder, pattern = "\\.(kreport|k2report)$", full.names = TRUE)
  
  # Filter out files containing "bracken" in the filename
  kreports_files <- kreports_files[!grepl("bracken_", basename(kreports_files))]
  
  if (length(kreports_files) == 0) {
    warning("No .kreport or .k2report files found in: ", kreports_folder)
    return(list(report_data = report_data_df, kreports_combined = tibble()))
  }
  
  cat("Processing", length(kreports_files), "kreport files from", kreports_folder, "\n")
  
  start_time <- Sys.time()
  # Initialize list to store raw kreports
  kreports_raw_list <- list()
  
  # Read all kreport files (sequential or parallel)
  if (USE_PARALLEL && exists("cl", envir = .GlobalEnv)) {
    cat("Reading kreport files in parallel...\n")
    kreports_raw_list <- foreach(i = seq_along(kreports_files), 
                                 .export = c("read_report2", "build_kraken_tree", "collapse.taxRanks", "delete_taxRanks_below", "parse_sample_info")) %dopar% {
                                   file <- kreports_files[i]
                                   
                                   # Read kreport file
                                   kreport_df <- read_report2(myfile = file, has_header = FALSE, add_taxRank_columns = FALSE, keep_taxRanks = c("D", "K", "P", "C", "O", "F", "G", "S", "S1", "S2", "S3"))
                                   
                                   if (is.null(kreport_df) || nrow(kreport_df) == 0) {
                                     return(NULL)
                                   }
                                   
                                   # Store raw kreport with file info
                                   list(
                                     file = file,
                                     data = kreport_df,
                                     sample_info = parse_sample_info(file)
                                   )
                                 }
  } else {
    cat("Reading kreport files sequentially...\n")
    for (i in seq_along(kreports_files)) {
      file <- kreports_files[i]
      cat("Processing kreport:", file, "\n")
      
      # Read kreport file
      kreport_df <- read_report2(myfile = file, has_header = FALSE, add_taxRank_columns = FALSE, keep_taxRanks = c("D", "K", "P", "C", "O", "F", "G", "S", "S1", "S2", "S3"))
      
      if (is.null(kreport_df) || nrow(kreport_df) == 0) {
        warning("Failed to read or empty file: ", file)
        next
      }
      
      # Store raw kreport with file info
      kreports_raw_list[[i]] <- list(
        file = file,
        data = kreport_df,
        sample_info = parse_sample_info(file)
      )
    }
  }
  
  # Remove NULL entries
  kreports_raw_list <- kreports_raw_list[!sapply(kreports_raw_list, is.null)]
  
  if (length(kreports_raw_list) == 0) {
    warning("No valid kreport files processed")
    return(list(report_data = report_data_df, kreports_combined = tibble()))
  }
  end_time <- Sys.time()
  cat("Finished reading kreport files in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
  
  # Now process each kreport for report data and filtering
  # Initialize fresh report_data for this folder
  if (include_subspecies) {
    report_data_df <- tibble(
      sample = character(), 
      unclassified_counts = integer(), 
      classified_counts = integer(), 
      confidence = integer(), 
      minimum_hits = integer(), 
      database_used = character(),
      species_count = integer(),
      subspecies_count = integer(),
      total_species_count = integer(),
      filtered_species_count = integer(),
      filtered_subspecies_count = integer(),
      filtered_total_species_count = integer()
    )
  } else {
    report_data_df <- tibble(
      sample = character(), 
      unclassified_counts = integer(), 
      classified_counts = integer(), 
      confidence = integer(), 
      minimum_hits = integer(), 
      database_used = character(),
      species_count = integer(),
      filtered_species_count = integer()
    )
  }
  kreports_list <- list()
  
  for (i in seq_along(kreports_raw_list)) {
    kreport_info <- kreports_raw_list[[i]]
    file <- kreport_info$file
    kreport_df <- kreport_info$data
    sample_info <- kreport_info$sample_info
    
    # Get unclassified reads count
    unclassified_reads <- kreport_df$cladeReads[kreport_df$taxID == 0]
    if (length(unclassified_reads) == 0) unclassified_reads <- 0
    
    # Add to report data
    if (include_subspecies) {
      report_data_df <- bind_rows(
        report_data_df,
        tibble(
          sample = sample_info$trimmed_name,
          unclassified_counts = unclassified_reads,
          classified_counts = kreport_df$cladeReads[kreport_df$taxID == 1], 
          confidence = sample_info$confidence,
          minimum_hits = sample_info$minimum_hits,
          database_used = sample_info$database_used,
          species_count = NA_integer_,
          subspecies_count = NA_integer_,
          total_species_count = NA_integer_,
          filtered_species_count = NA_integer_,
          filtered_subspecies_count = NA_integer_,
          filtered_total_species_count = NA_integer_
        )
      )
    } else {
      report_data_df <- bind_rows(
        report_data_df,
        tibble(
          sample = sample_info$trimmed_name,
          unclassified_counts = unclassified_reads,
          classified_counts = kreport_df$cladeReads[kreport_df$taxID == 1], 
          confidence = sample_info$confidence,
          minimum_hits = sample_info$minimum_hits,
          database_used = sample_info$database_used,
          species_count = NA_integer_,
          filtered_species_count = NA_integer_
        )
      )
    }
    
    # Process kreport data - always read all levels, then filter based on include_subspecies
    # First, get all species and subspecies data
    all_species_data <- kreport_df %>%
      filter(taxRank %in% c("S", "S1", "S2", "S3")) %>%
      mutate(
        name = str_remove(name, "^s\\d*_"),  # remove leading s_, s1_, s2_, s3_, etc.
        name = str_remove_all(name, "'") # remove apostrophes
      )
    
    # Filter based on include_subspecies parameter
    if (include_subspecies) {
      # Include both species and subspecies levels
      species_level_data <- all_species_data
    } else {
      # Include only species level
      species_level_data <- all_species_data %>%
        filter(taxRank == "S")
    }
    
    # Calculate initial species counts (before filtering)
    if (include_subspecies) {
      species_count <- all_species_data %>% filter(taxRank == "S") %>% nrow()
      subspecies_count <- all_species_data %>% filter(taxRank %in% c("S1", "S2", "S3")) %>% nrow()
      total_species_count <- nrow(all_species_data)
    } else {
      species_count <- all_species_data %>% filter(taxRank == "S") %>% nrow()
    }
    
    # Apply filtering to species-level data
    filtered_data <- filter_kreport_data(species_level_data, sample_info$trimmed_name)
    
    # Calculate filtered species counts
    if (include_subspecies) {
      filtered_species_count <- filtered_data %>% filter(taxRank == "S") %>% nrow()
      filtered_subspecies_count <- filtered_data %>% filter(taxRank %in% c("S1", "S2", "S3")) %>% nrow()
      filtered_total_species_count <- nrow(filtered_data)
      
      # Update the last row of report_data with species counts
      report_data_df[nrow(report_data_df), "species_count"] <- species_count
      report_data_df[nrow(report_data_df), "subspecies_count"] <- subspecies_count
      report_data_df[nrow(report_data_df), "total_species_count"] <- total_species_count
      report_data_df[nrow(report_data_df), "filtered_species_count"] <- filtered_species_count
      report_data_df[nrow(report_data_df), "filtered_subspecies_count"] <- filtered_subspecies_count
      report_data_df[nrow(report_data_df), "filtered_total_species_count"] <- filtered_total_species_count
    } else {
      filtered_species_count <- nrow(filtered_data)
      
      # Update the last row of report_data with species counts
      report_data_df[nrow(report_data_df), "species_count"] <- species_count
      report_data_df[nrow(report_data_df), "filtered_species_count"] <- filtered_species_count
    }
    
    # Prepare final processed data (common for both subspecies and species-only)
    processed_kreport <- filtered_data %>%
      select(cladeReads, minimizers, distinct_minimizers, taxID, name, taxRank, taxLineage) %>%
      rename_with(~ paste0(sample_info$trimmed_name, "_", .), 
                  c(cladeReads, minimizers, distinct_minimizers))
    
    # Check for duplicate taxIDs
    duplicate_taxids <- processed_kreport %>% 
      count(taxID) %>% 
      filter(n > 1)
    
    if (nrow(duplicate_taxids) > 0) {
      warning("Sample ", sample_info$trimmed_name, " has ", nrow(duplicate_taxids), 
              " duplicate taxIDs: ", paste(duplicate_taxids$taxID, collapse = ", "))
    }
    
    # Add processed kreport to list
    kreports_list[[i]] <- processed_kreport
  }
  
  # Remove NULL entries
  kreports_list <- kreports_list[!sapply(kreports_list, is.null)]
  
  if (length(kreports_list) == 0) {
    warning("No valid kreport files processed")
    # Return properly structured empty tibble
    empty_kreports <- tibble(
      taxID = integer(0),
      name = character(0),
      taxRank = character(0),
      taxLineage = character(0)
    )
    return(list(report_data = report_data_df, kreports_combined = empty_kreports))
  }
  
  # Combine all kreports
  kreports_combined <- reduce(kreports_list, full_join, by = c("taxID", "name", "taxRank", "taxLineage")) %>%
    select(taxID, name, taxRank, taxLineage, everything())
  
  return(list(report_data = report_data_df, kreports_combined = kreports_combined))
}

# Simplified function to process bracken files
process_bracken_folder <- function(bracken_folder) {
  if (!dir.exists(bracken_folder)) {
    warning("Directory does not exist: ", bracken_folder)
    return(tibble())
  }
  
  bracken_files <- list.files(bracken_folder, pattern = ".*bracken.*\\.txt$", full.names = TRUE)
  
  if (length(bracken_files) == 0) {
    warning("No bracken .txt files found in: ", bracken_folder)
    return(tibble())
  }
  
  cat("Processing", length(bracken_files), "bracken files from", bracken_folder, "\n")
  
  bracken_list <- map(bracken_files, function(file) {
    tryCatch({
      # Ensure libraries are available for parallel processing
      
      sample_info <- parse_sample_info(file)
      
      bracken_df <- read.table(file, sep = "\t", header = T, comment.char = "", quote = "") %>%
        select(name, taxonomy_id, new_est_reads) %>%
        rename(taxID = taxonomy_id) %>%
        rename_with(~ paste0(sample_info$trimmed_name, "_bracken_reads"), 
                    .cols = "new_est_reads")
      
      return(bracken_df)
    }, error = function(e) {
      cat("ERROR processing bracken file:", basename(file), "-", e$message, "\n")
      return(NULL)
    })
  })
  
  # Combine all bracken reports
  # Remove NULL entries first
  bracken_list <- bracken_list[!sapply(bracken_list, is.null)]
  
  if (length(bracken_list) > 0) {
    bracken_combined <- reduce(bracken_list, full_join, by = c("name", "taxID"))
  } else {
    bracken_combined <- tibble()
  }
  
  bracken_combined <- bracken_combined %>% mutate(name = str_remove_all(name, "'"))  # remove apostrophes
  
  return(bracken_combined)
}

# Function to merge kreport and bracken data
merge_data <- function(kreports_data, bracken_data = NULL, add_params = ADD_PARAMS, report_data = NULL, dataset_name = NULL) {
  # Merge with bracken data if available
  if (!is.null(bracken_data) && nrow(bracken_data) > 0) {
    # Only join bracken entries whose taxID is present in kreports_data to avoid re-adding filtered species
    bracken_data_filtered <- bracken_data %>%
      filter(taxID %in% kreports_data$taxID)
    merged_data <- kreports_data %>%
      full_join(bracken_data_filtered, by = c("name", "taxID"))
  } else {
    merged_data <- kreports_data
  }
  
  # Check for duplicate taxIDs
  duplicate_taxids <- merged_data %>% 
      count(taxID) %>% 
      filter(n > 1)
  if (nrow(duplicate_taxids) > 0) {
    warning("Duplicate taxIDs found: ", paste(duplicate_taxids$taxID, collapse = ", "))
  }
  
  # Reorder columns
  merged_data <- merged_data %>%
    select(taxID, name, taxRank, taxLineage, everything())
  
  # Create long format
  merged_long <- merged_data %>%
    pivot_longer(
      cols = -c(taxID, name, taxRank, taxLineage),
      names_to = c("sample", "measurement"),
      names_pattern = "^(.*?)_(cladeReads|minimizers|distinct_minimizers|bracken_reads)$",
      values_to = "value"
    ) %>%
    pivot_wider(
      names_from = measurement,
      values_from = value
    ) 
  
  # Filter out rows where there is a bracken_read but cladeReads is NA
  # This is to address the issue that we are filtering kreports by distinct minimizer 
  # but did not apply the same filter to bracken reports
 if ("bracken_reads" %in% names(merged_long)) {
  merged_long <- merged_long |>
    dplyr::filter(!(is.na(cladeReads) & !is.na(bracken_reads)))
} else if ("brackenReads" %in% names(merged_long)) {
  merged_long <- merged_long |>
    dplyr::filter(!(is.na(cladeReads) & !is.na(brackenReads)))
} else {
  merged_long <- merged_long |>
    dplyr::filter(!is.na(cladeReads))
}
  
  # Add parameter columns if requested
  if (add_params) {
    # Join with report_data to get parameter information
    if (!is.null(report_data) && nrow(report_data) > 0) {
      # Create a lookup table from report_data for parameter information
      param_lookup <- report_data %>%
        select(sample, confidence, minimum_hits, database_used) %>%
        distinct() %>%
        rename(
          confidence_levels = confidence,
          minimum_hit_groups = minimum_hits
        ) %>%
        mutate(
          # Determine human_reads based on dataset name
          human_reads = if (!is.null(dataset_name)) {
            grepl("nonhuman", dataset_name, ignore.case = TRUE)
          } else {
            # Fallback to original logic if dataset_name not provided
            grepl("human|nonhuman", database_used, ignore.case = TRUE)
          }
        )
      
      # Join the parameter information with merged_long
      merged_long <- merged_long %>%
        left_join(param_lookup, by = "sample")
      
      cat("Added parameter columns: confidence_levels, minimum_hit_groups, human_reads\n")
    } else {
      warning("report_data is NULL or empty - cannot add parameter columns")
    }
  }
  
  return(list(merged = merged_data, merged_long = merged_long))
}

# Function to add summary statistics
add_summary_stats <- function(merged_data) {
  cladeReads_cols <- grep("_cladeReads$", colnames(merged_data), value = TRUE)
  bracken_reads_cols <- grep("_bracken_reads$", colnames(merged_data), value = TRUE)
  
  result <- merged_data
  
  if (length(cladeReads_cols) > 0) {
    # Create numeric matrix of cladeReads columns
    clade_mat <- as.matrix(select(result, all_of(cladeReads_cols)))
    # Compute row means and max while ignoring NA; if all NA, return NA
    clade_means <- rowMeans(clade_mat, na.rm = TRUE)
    clade_maxs <- apply(clade_mat, 1, function(x) if (all(is.na(x))) 0 else max(x, na.rm = TRUE))
    result <- result %>%
      mutate(
        cladeReads_mean = clade_means,
        cladeReads_max  = clade_maxs
      )
  }
  
  if (length(bracken_reads_cols) > 0) {
    # Extract matrix of bracken read columns
    bracken_mat <- as.matrix(select(result, all_of(bracken_reads_cols)))
    # Compute row means, replace NaN/NA (all-NA rows) with 0
    bracken_means <- rowMeans(bracken_mat, na.rm = TRUE)
    bracken_means[is.nan(bracken_means) | is.na(bracken_means)] <- 0
    # Compute row max, return 0 when all values are NA
    bracken_max <- apply(bracken_mat, 1, function(x) {
      if (all(is.na(x))) 0 else max(x, na.rm = TRUE)
    })
    result <- result %>%
      mutate(
        bracken_reads_mean = bracken_means,
        bracken_reads_max  = bracken_max
      )
    cat("Added bracken summary statistics: bracken_reads_mean, bracken_reads_max\n")
  }
  
  return(result)
}

# Function to get complete species list with topN frequency included
get_species_list_with_frequency <- function(merged_data, count_column_pattern = "_cladeReads$", 
                                            top_n) {
  count_cols <- grep(count_column_pattern, colnames(merged_data), value = TRUE)
  
  if (length(count_cols) == 0) {
    warning("No columns matching pattern found: ", count_column_pattern)
    # Return properly structured empty tibble with required columns
    return(tibble(
      name = character(0),
      taxID = integer(0),
      taxRank = character(0),
      taxLineage = character(0),
      Freq = integer(0)
    ))
  }
  
  # Get top species from each sample
  top_species_list <- map_dfr(count_cols, function(sample_col) {
    temp_df <- merged_data %>%
      select(name, taxID, taxRank, taxLineage, count = all_of(sample_col)) %>%
      arrange(desc(count)) %>%
      slice_head(n = top_n)
    
    return(temp_df)
  })
  
  # Create frequency table for top species
  top_species_freq <- top_species_list %>%
    count(taxID, name = "topN_freq") %>%
    filter(topN_freq >= 1)
  
  # Create complete species list with same format but include topN frequency
  complete_list_df <- merged_data %>%
    select(name, taxID, taxRank, taxLineage) %>%
    distinct() %>%
    left_join(top_species_freq, by = "taxID") %>%
    mutate(Freq = ifelse(is.na(topN_freq), 0, topN_freq)) %>%
    select(-topN_freq) %>%
    arrange(desc(Freq), name)
  
  # Add cladeReads_mean and cladeReads_max if they exist in the data
  if ("cladeReads_mean" %in% colnames(merged_data) && "cladeReads_max" %in% colnames(merged_data)) {
    # Get unique mean and max values for each species
    summary_stats <- merged_data %>% select(taxID, cladeReads_mean, cladeReads_max) %>%
      distinct()
    # Join with the species list
    complete_list_df <- complete_list_df %>%
      left_join(summary_stats, by = "taxID")
  }
  
  # Add bracken_reads_mean and bracken_reads_max if they exist in the data
  if ("bracken_reads_mean" %in% colnames(merged_data) && "bracken_reads_max" %in% colnames(merged_data)) {
    # Get unique bracken mean and max values for each species
    bracken_stats <- merged_data %>% select(taxID, bracken_reads_mean, bracken_reads_max) %>%
      distinct()

    # Join with the species list
    complete_list_df  <- complete_list_df %>%
      left_join(bracken_stats, by = "taxID")
    
    cat("Added bracken statistics to species list: bracken_reads_mean, bracken_reads_max\n")
  }

  return(complete_list_df)
}

# Function to parse sample metadata files
parse_sample_metadata <- function(input_dir) {
  # Look for sample metadata CSV files
  metadata_files <- list.files(input_dir, pattern = "sample_metadata.*\\.csv$", full.names = TRUE, recursive = FALSE)
  
  if (length(metadata_files) == 0) {
    cat("No sample metadata files found in input directory\n")
    return(NULL)
  }
  
  cat("Found", length(metadata_files), "sample metadata files\n")
  
  all_metadata <- list()
  
  for (metadata_file in metadata_files) {
    cat("Processing metadata file:", basename(metadata_file), "\n")
    
    tryCatch({
      metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
      
      # Check for required columns (case-insensitive)
      required_cols <- c("sample", "condition")
      col_names_lower <- tolower(colnames(metadata))
      
      # Find matching columns
      sample_col <- which(col_names_lower %in% c("sample", "samples"))
      condition_col <- which(col_names_lower %in% c("condition", "conditions"))
      
      if (length(sample_col) == 0 || length(condition_col) == 0) {
        cat("Warning: Metadata file missing required columns. Expected 'Sample' and 'Condition'\n")
        cat("Available columns:", paste(colnames(metadata), collapse = ", "), "\n")
        next
      }
      
      # Standardize column names
      colnames(metadata)[sample_col[1]] <- "sample"
      colnames(metadata)[condition_col[1]] <- "condition"
      
      # Remove rows with missing sample names
      metadata <- metadata[!is.na(metadata$sample) & metadata$sample != "", ]
      
      # Apply consistent sample name trimming using parse_sample_info
      if (nrow(metadata) > 0) {
        metadata$sample <- sapply(metadata$sample, function(x) {
          sample_info <- parse_sample_info(x)
          # Remove _S and everything after to match runtime data
          return(sub("(_S\\d+).*", "", sample_info$trimmed_name))
        })
      }
      
      if (nrow(metadata) > 0) {
        all_metadata[[basename(metadata_file)]] <- metadata
      }
      
    }, error = function(e) {
      cat("Error reading", basename(metadata_file), ":", e$message, "\n")
    })
  }
  
  # Combine all metadata
  if (length(all_metadata) > 0) {
    combined_metadata <- do.call(rbind, all_metadata)
    
    # Remove duplicates, keeping latest values
    final_metadata <- combined_metadata %>%
      group_by(sample) %>%
      slice_tail(n = 1) %>%
      ungroup()
    
    return(final_metadata)
  } else {
    return(NULL)
  }
}

# Function to parse runtime files
parse_runtime_files <- function(runtime_folder) {
  if (!dir.exists(runtime_folder)) {
    cat("Runtime folder not found:", runtime_folder, "\n")
    return(NULL)
  }
  
  runtime_files <- list.files(runtime_folder, pattern = "^runtime", full.names = TRUE)
  
  if (length(runtime_files) == 0) {
    cat("No runtime files found in:", runtime_folder, "\n")
    return(NULL)
  }
  
  cat("Found", length(runtime_files), "runtime files\n")
  
  all_runtime_data <- list()
  
  for (runtime_file in runtime_files) {
    cat("Processing runtime file:", basename(runtime_file), "\n")
    
    # Read the file
    lines <- readLines(runtime_file, warn = FALSE)  # Suppress incomplete final line warnings
    
    # Initialize data structures for this file
    runtime_data <- data.frame(
      sample = character(),
      cutadapt_star_samtools_rt = numeric(),
      kraken2_unaligned_rt = numeric(),
      kraken2_nonhuman_rt = numeric(),
      is_first_k2_unaligned = logical(),
      is_first_k2_nonhuman = logical(),
      stringsAsFactors = FALSE
    )
    
    # Track first samples for database loading
    first_k2_unaligned <- TRUE
    first_k2_nonhuman <- TRUE
    
    # Parse each line
    for (line in lines) {
      # Parse Cutadapt/STAR/Samtools runtime
      if (grepl("^Cutadapt/STAR/Samtools \\(.*\\); RT:", line)) {
        sample_match <- regmatches(line, regexpr("\\(([^)]+)\\)", line))
        sample_name_raw <- gsub("[()]", "", sample_match)
        # Use parse_sample_info to get consistent trimmed name, then remove _S and everything after
        sample_info <- parse_sample_info(sample_name_raw)
        sample_name <- sub("(_S\\d+).*", "", sample_info$trimmed_name)
        rt_match <- regmatches(line, regexpr("RT: (\\d+) seconds", line))
        rt_text <- gsub("RT: | seconds", "", rt_match)
        rt_seconds <- if (length(rt_text) > 0 && rt_text != "") as.numeric(rt_text) else NA_real_        # Find existing row or create new one
        existing_row <- which(runtime_data$sample == sample_name)
        if (length(existing_row) > 0) {
          # Update existing row (take latest runtime)
          runtime_data$cutadapt_star_samtools_rt[existing_row] <- rt_seconds
        } else {
          # Add new row
          new_row <- data.frame(
            sample = sample_name,
            cutadapt_star_samtools_rt = rt_seconds,
            kraken2_unaligned_rt = NA,
            kraken2_nonhuman_rt = NA,
            is_first_k2_unaligned = FALSE,
            is_first_k2_nonhuman = FALSE,
            stringsAsFactors = FALSE
          )
          runtime_data <- rbind(runtime_data, new_row)
        }
      }
      
      # Parse Kraken2 unaligned runtime
      if (grepl("^Kraken2 unaligned \\(.*\\); RT:", line)) {
        sample_match <- regmatches(line, regexpr("\\(([^)]+)\\)", line))
        sample_name_raw <- gsub("[()]", "", sample_match)
        # Use parse_sample_info to get consistent trimmed name, then remove _S and everything after
        sample_info <- parse_sample_info(sample_name_raw)
        sample_name <- sub("(_S\\d+).*", "", sample_info$trimmed_name)
        rt_match <- regmatches(line, regexpr("RT: (\\d+) seconds", line))
        rt_text <- gsub("RT: | seconds", "", rt_match)
        rt_seconds <- if (length(rt_text) > 0 && rt_text != "") as.numeric(rt_text) else NA_real_
        
        # Find existing row or create new one
        existing_row <- which(runtime_data$sample == sample_name)
        if (length(existing_row) > 0) {
          # Update existing row (take latest runtime)
          runtime_data$kraken2_unaligned_rt[existing_row] <- rt_seconds
          runtime_data$is_first_k2_unaligned[existing_row] <- first_k2_unaligned
        } else {
          # Add new row
          new_row <- data.frame(
        sample = sample_name,
        cutadapt_star_samtools_rt = NA,
        kraken2_unaligned_rt = rt_seconds,
        kraken2_nonhuman_rt = NA,
        is_first_k2_unaligned = first_k2_unaligned,
        is_first_k2_nonhuman = FALSE,
        stringsAsFactors = FALSE
          )
          runtime_data <- rbind(runtime_data, new_row)
        }
        
        # Mark that we've seen the first K2 unaligned sample
        if (first_k2_unaligned) {
          cat("  First K2 unaligned sample in", basename(runtime_file), ":", sample_name, "\n")
          first_k2_unaligned <- FALSE
        }
      }
      
      # Parse Kraken2 non-human runtime
      if (grepl("^Kraken2 non-human \\(.*\\); RT:", line)) {
        sample_match <- regmatches(line, regexpr("\\(([^)]+)\\)", line))
        sample_name_raw <- gsub("[()]", "", sample_match)
        # Use parse_sample_info to get consistent trimmed name, then remove _S and everything after
        sample_info <- parse_sample_info(sample_name_raw)
        sample_name <- sub("(_S\\d+).*", "", sample_info$trimmed_name)
        rt_match <- regmatches(line, regexpr("RT: (\\d+) seconds", line))
        rt_text <- gsub("RT: | seconds", "", rt_match)
        rt_seconds <- if (length(rt_text) > 0 && rt_text != "") as.numeric(rt_text) else NA_real_
        
        # Find existing row or create new one
        existing_row <- which(runtime_data$sample == sample_name)
        if (length(existing_row) > 0) {
          # Update existing row (take latest runtime)
          runtime_data$kraken2_nonhuman_rt[existing_row] <- rt_seconds
          runtime_data$is_first_k2_nonhuman[existing_row] <- first_k2_nonhuman
        } else {
          # Add new row
          new_row <- data.frame(
            sample = sample_name,
            cutadapt_star_samtools_rt = NA,
            kraken2_unaligned_rt = NA,
            kraken2_nonhuman_rt = rt_seconds,
            is_first_k2_unaligned = FALSE,
            is_first_k2_nonhuman = first_k2_nonhuman,
            stringsAsFactors = FALSE
          )
          runtime_data <- rbind(runtime_data, new_row)
        }
        
        # Mark that we've seen the first K2 non-human sample
        if (first_k2_nonhuman) {
          cat("  First K2 non-human sample in", basename(runtime_file), ":", sample_name, "\n")
          first_k2_nonhuman <- FALSE
        }
      }
    }
    
    all_runtime_data[[basename(runtime_file)]] <- runtime_data
  }
  
  # Combine all runtime data from different files
  if (length(all_runtime_data) > 0) {
    combined_runtime <- do.call(rbind, all_runtime_data)
    
    # Remove any rows with empty sample names
    combined_runtime <- combined_runtime[!is.na(combined_runtime$sample) & combined_runtime$sample != "", ]
    
    # If there are duplicates (same sample from multiple files), keep the latest values
    # Group by sample and take the non-NA values, preferring later files
    # For boolean flags, use OR logic (TRUE if any file marks it as first)
    final_runtime <- combined_runtime %>%
      group_by(sample) %>%
      summarise(
        cutadapt_star_samtools_rt = ifelse(length(cutadapt_star_samtools_rt[!is.na(cutadapt_star_samtools_rt)]) > 0,
                                           last(cutadapt_star_samtools_rt[!is.na(cutadapt_star_samtools_rt)]), NA_real_),
        kraken2_unaligned_rt = ifelse(length(kraken2_unaligned_rt[!is.na(kraken2_unaligned_rt)]) > 0,
                                      last(kraken2_unaligned_rt[!is.na(kraken2_unaligned_rt)]), NA_real_),
        kraken2_nonhuman_rt = ifelse(length(kraken2_nonhuman_rt[!is.na(kraken2_nonhuman_rt)]) > 0,
                                     last(kraken2_nonhuman_rt[!is.na(kraken2_nonhuman_rt)]), NA_real_),
        is_first_k2_unaligned = any(is_first_k2_unaligned, na.rm = TRUE),
        is_first_k2_nonhuman = any(is_first_k2_nonhuman, na.rm = TRUE),
        .groups = 'drop'
      )
    
    return(final_runtime)
  } else {
    return(NULL)
  }
}

# Function to parse read statistics files
parse_rrstats_files <- function(rrstats_folder) {
  if (!dir.exists(rrstats_folder)) {
    cat("Read stats folder not found:", rrstats_folder, "\n")
    return(NULL)
  }
  
  rrstats_files <- list.files(rrstats_folder, pattern = "run_read_stats.*\\.txt$", full.names = TRUE)
  
  if (length(rrstats_files) == 0) {
    cat("No read stats files found in:", rrstats_folder, "\n")
    return(NULL)
  }
  
  cat("Found", length(rrstats_files), "read stats files\n")
  
  all_rrstats_data <- list()
  
  for (rrstats_file in rrstats_files) {
    cat("Processing read stats file:", basename(rrstats_file), "\n")
    
    # Read the file as tab-separated
    tryCatch({
      rrstats_data <- read.table(rrstats_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                                 skip = 0, fill = TRUE, comment.char = "")
      
      # Remove any completely empty rows
      rrstats_data <- rrstats_data[!apply(rrstats_data == "" | is.na(rrstats_data), 1, all), ]
      
      # Skip if no data after cleaning
      if (nrow(rrstats_data) == 0) {
        cat("No valid data found in", basename(rrstats_file), "\n")
        next
      }
      
      # Set column names based on the header structure
      colnames(rrstats_data) <- c("sample", "num_input_reads", "avg_input_read_length", 
                                  "uniquely_mapped_reads", "avg_mapped_length", 
                                  "reads_assigned_to_genes", "percent_mapped_reads", 
                                  "uniquely_mapped_percent", "less_than_10M_unique")
      
      # Convert numeric columns with better error handling
      numeric_cols <- c("num_input_reads", "avg_input_read_length", "uniquely_mapped_reads", 
                        "avg_mapped_length", "reads_assigned_to_genes", "percent_mapped_reads", 
                        "uniquely_mapped_percent")
      
      for (col in numeric_cols) {
        if (col %in% colnames(rrstats_data)) {
          # Convert to numeric with proper error handling
          numeric_vals <- as.numeric(as.character(rrstats_data[[col]]))
          if (any(is.na(numeric_vals) & !is.na(rrstats_data[[col]]))) {
            cat("Note: Some non-numeric values found in column", col, "and converted to NA\n")
          }
          rrstats_data[[col]] <- numeric_vals
        }
      }
      
      # Convert boolean column with better handling
      if ("less_than_10M_unique" %in% colnames(rrstats_data)) {
        rrstats_data$less_than_10M_unique <- as.logical(rrstats_data$less_than_10M_unique)
      }
      
      # Remove rows where sample name is NA or empty
      rrstats_data <- rrstats_data[!is.na(rrstats_data$sample) & rrstats_data$sample != "", ]
      
      # Apply consistent sample name trimming using parse_sample_info
      if (nrow(rrstats_data) > 0) {
        rrstats_data$sample <- sapply(rrstats_data$sample, function(x) {
          sample_info <- parse_sample_info(x)
          # Remove _S and everything after to match runtime data
          return(sub("(_S\\d+).*", "", sample_info$trimmed_name))
        })
      }
      
      if (nrow(rrstats_data) > 0) {
        all_rrstats_data[[basename(rrstats_file)]] <- rrstats_data
      }
      
    }, error = function(e) {
      cat("Error reading", basename(rrstats_file), ":", e$message, "\n")
    })
  }
  
  # Combine all read stats data from different files
  if (length(all_rrstats_data) > 0) {
    combined_rrstats <- do.call(rbind, all_rrstats_data)
    
    # Remove any remaining NA sample names
    combined_rrstats <- combined_rrstats[!is.na(combined_rrstats$sample) & combined_rrstats$sample != "", ]
    
    # If there are duplicates (same sample from multiple files), keep the latest values
    final_rrstats <- combined_rrstats %>%
      group_by(sample) %>%
      summarise(
        num_input_reads = ifelse(length(num_input_reads[!is.na(num_input_reads)]) > 0,
                                 last(num_input_reads[!is.na(num_input_reads)]), NA_real_),
        avg_input_read_length = ifelse(length(avg_input_read_length[!is.na(avg_input_read_length)]) > 0,
                                       last(avg_input_read_length[!is.na(avg_input_read_length)]), NA_real_),
        uniquely_mapped_reads = ifelse(length(uniquely_mapped_reads[!is.na(uniquely_mapped_reads)]) > 0,
                                       last(uniquely_mapped_reads[!is.na(uniquely_mapped_reads)]), NA_real_),
        avg_mapped_length = ifelse(length(avg_mapped_length[!is.na(avg_mapped_length)]) > 0,
                                   last(avg_mapped_length[!is.na(avg_mapped_length)]), NA_real_),
        reads_assigned_to_genes = ifelse(length(reads_assigned_to_genes[!is.na(reads_assigned_to_genes)]) > 0,
                                         last(reads_assigned_to_genes[!is.na(reads_assigned_to_genes)]), NA_real_),
        percent_mapped_reads = ifelse(length(percent_mapped_reads[!is.na(percent_mapped_reads)]) > 0,
                                      last(percent_mapped_reads[!is.na(percent_mapped_reads)]), NA_real_),
        uniquely_mapped_percent = ifelse(length(uniquely_mapped_percent[!is.na(uniquely_mapped_percent)]) > 0,
                                         last(uniquely_mapped_percent[!is.na(uniquely_mapped_percent)]), NA_real_),
        less_than_10M_unique = ifelse(length(less_than_10M_unique[!is.na(less_than_10M_unique)]) > 0,
                                      last(less_than_10M_unique[!is.na(less_than_10M_unique)]), NA),
        .groups = 'drop'
      )
    
    return(final_rrstats)
  } else {
    return(NULL)
  }
}

# RPM count matrix helpers
read_rpm_count_matrix <- function(file_path) {
  if (!file.exists(file_path)) {
    warning("Count matrix file not found: ", file_path)
    return(NULL)
  }
  count_df <- tryCatch({
    utils::read.delim(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }, error = function(e) {
    warning("Failed to read count matrix: ", basename(file_path), " - ", e$message)
    return(NULL)
  })
  if (is.null(count_df) || nrow(count_df) == 0) {
    return(NULL)
  }
  if (!"Geneid" %in% colnames(count_df)) {
    colnames(count_df)[1] <- "Geneid"
  }
  return(count_df)
}

merge_count_matrices <- function(rpm_dir) {
  if (is.null(rpm_dir) || !dir.exists(rpm_dir)) {
    warning("RPM directory not found: ", rpm_dir)
    return(NULL)
  }
  rpm_files <- list.files(rpm_dir, pattern = "countMat.*\\.(txt|csv)$", full.names = TRUE)
  if (length(rpm_files) == 0) {
    warning("No countMat files found in: ", rpm_dir)
    return(NULL)
  }
  cat("Merging", length(rpm_files), "count matrices from", rpm_dir, "\n")
  count_list <- lapply(rpm_files, read_rpm_count_matrix)
  count_list <- count_list[!sapply(count_list, is.null)]
  if (length(count_list) == 0) {
    return(NULL)
  }
  combined <- reduce(count_list, full_join, by = "Geneid")
  combined <- combined %>%
    mutate(Geneid = as.character(Geneid))
  colnames(combined) <- gsub("_count$", "", colnames(combined))
  if (any(duplicated(colnames(combined)))) {
    warning("Duplicate column names detected after merge; making names unique")
    colnames(combined) <- make.unique(colnames(combined))
  }
  return(combined)
}

filter_gene_counts <- function(count_mat, rrna_file,
                               protein_coding_filter = PROTEIN_CODING_FILTER,
                               remove_mt_genes = REMOVE_MT_GENES,
                               go_annotation_filter = GO_ANNOTATION_FILTER) {
  if (is.null(count_mat) || nrow(count_mat) == 0) {
    return(NULL)
  }
  cat("Initial gene count:", nrow(count_mat), "\n")
  # Step 1: Remove rRNA genes
  if (!is.null(rrna_file) && file.exists(rrna_file)) {
    cat("Removing rRNA genes...\n")
    rrnalist <- read.csv(rrna_file)
    if ("Approved.symbol" %in% colnames(rrnalist)) {
      count_mat <- subset(count_mat, !(Geneid %in% rrnalist$Approved.symbol))
    }
    cat("  After rRNA removal:", nrow(count_mat), "genes\n")
  } else {
    cat("  No rRNA gene list found, skipping rRNA filtering\n")
  }
  # Step 2: Drop rows which are all 0
  cat("Removing zero-count genes...\n")
  numeric_cols <- count_mat %>% select(where(is.numeric))
  keep_rows <- rowSums(numeric_cols, na.rm = TRUE) > 0
  count_mat <- count_mat[keep_rows, , drop = FALSE]
  cat("  After zero removal:", nrow(count_mat), "genes\n")
  # Step 3: Merge alternative splice genes
  cat("Merging alternative splice variants (-AS#)...\n")
  count_mat$Geneid_trim <- gsub("-AS[0-9]", "", count_mat$Geneid)
  count_mat <- count_mat %>%
    select(-Geneid) %>%
    group_by(Geneid_trim) %>%
    summarize(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
  colnames(count_mat)[1] <- "Geneid"
  cat("  After AS merging:", nrow(count_mat), "genes\n")
  # Step 4: Merge duplicate gene variants
  cat("Merging duplicate gene variants (_#)...\n")
  count_mat$Geneid_trim2 <- gsub("_[0-9]", "", count_mat$Geneid)
  count_mat <- count_mat %>%
    select(-Geneid) %>%
    group_by(Geneid_trim2) %>%
    summarize(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
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

  # Step 7: Filter to protein-coding genes (optional)
  if (protein_coding_filter) {
    cat("Filtering for protein-coding genes (biomaRt)...\n")
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
      cat("  biomaRt not installed; attempting installation...\n")
      tryCatch({
        install.packages("biomaRt", repos = "https://cloud.r-project.org")
      }, error = function(e) {
        warning("Failed to install biomaRt: ", e$message)
      })
    }

    if (requireNamespace("biomaRt", quietly = TRUE)) {
      gene_symbols <- as.character(count_mat$Geneid)
      gene_info <- tryCatch({
        ensembl_live <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        biomaRt::getBM(
          attributes = c("hgnc_symbol", "gene_biotype"),
          filters = "hgnc_symbol",
          values = gene_symbols,
          mart = ensembl_live
        )
      }, error = function(e) {
        cat("  Live Ensembl unavailable; trying Archive Release 114...\n")
        archive_host <- "https://may2025.archive.ensembl.org"
        ensembl_archive <- biomaRt::useMart(
          biomart = "ensembl",
          dataset = "hsapiens_gene_ensembl",
          host = archive_host
        )
        biomaRt::getBM(
          attributes = c("hgnc_symbol", "gene_biotype"),
          filters = "hgnc_symbol",
          values = gene_symbols,
          mart = ensembl_archive
        )
      })

      if (!is.null(gene_info) && nrow(gene_info) > 0) {
        protein_coding_genes <- unique(gene_info$hgnc_symbol[gene_info$gene_biotype == "protein_coding"])
        protein_coding_genes <- protein_coding_genes[!is.na(protein_coding_genes) & nzchar(protein_coding_genes)]
        if (length(protein_coding_genes) > 0) {
          count_mat <- count_mat[count_mat$Geneid %in% protein_coding_genes, , drop = FALSE]
          cat("  After protein-coding filter:", nrow(count_mat), "genes\n")
        } else {
          warning("Protein-coding filter returned no matching genes; skipping this step")
        }
      } else {
        warning("No gene biotype annotations returned; skipping protein-coding filter")
      }
    } else {
      warning("biomaRt unavailable; skipping protein-coding filter")
    }
  } else {
    cat("Protein-coding filter skipped (--no-protein-coding)\n")
  }

  # Step 8: Remove mitochondrial genes (optional)
  if (remove_mt_genes) {
    cat("Removing mitochondrial genes (MT-)...\n")
    count_mat <- count_mat %>% filter(!grepl("^MT-", Geneid))
    cat("  After MT removal:", nrow(count_mat), "genes\n")
  } else {
    cat("Mitochondrial gene removal skipped (--no-remove-mt)\n")
  }

  # Step 9: Keep only genes with any GO annotation (optional)
  if (go_annotation_filter) {
    cat("Filtering genes with GO annotations (org.Hs.eg.db)...\n")
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      cat("  org.Hs.eg.db not installed; attempting installation...\n")
      tryCatch({
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
        BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)
      }, error = function(e) {
        warning("Failed to install org.Hs.eg.db: ", e$message)
      })
    }

    if (requireNamespace("org.Hs.eg.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
      go_mapped <- tryCatch({
        AnnotationDbi::select(
          org.Hs.eg.db::org.Hs.eg.db,
          keys = as.character(count_mat$Geneid),
          columns = "GO",
          keytype = "SYMBOL"
        )
      }, error = function(e) {
        warning("GO annotation query failed: ", e$message)
        NULL
      })

      if (!is.null(go_mapped) && nrow(go_mapped) > 0) {
        annotated_genes <- unique(go_mapped$SYMBOL[!is.na(go_mapped$GO)])
        annotated_genes <- annotated_genes[!is.na(annotated_genes) & nzchar(annotated_genes)]
        if (length(annotated_genes) > 0) {
          count_mat <- count_mat[count_mat$Geneid %in% annotated_genes, , drop = FALSE]
          cat("  After GO annotation filter:", nrow(count_mat), "genes\n")
        } else {
          warning("GO annotation filter returned no matching genes; skipping this step")
        }
      } else {
        warning("No GO annotation data returned; skipping GO filter")
      }
    } else {
      warning("org.Hs.eg.db or AnnotationDbi unavailable; skipping GO annotation filter")
    }
  } else {
    cat("GO annotation filter skipped (--no-go-annotation)\n")
  }

  return(as.data.frame(count_mat))
}

process_gene_counts <- function(rpm_dir, counts_output_dir, rrna_file,
                                protein_coding_filter = PROTEIN_CODING_FILTER,
                                remove_mt_genes = REMOVE_MT_GENES,
                                go_annotation_filter = GO_ANNOTATION_FILTER) {
  combined_counts <- merge_count_matrices(rpm_dir)
  if (is.null(combined_counts)) {
    warning("No count matrices merged; skipping gene filtering")
    return(NULL)
  }
  if (!dir.exists(counts_output_dir)) {
    dir.create(counts_output_dir, recursive = TRUE)
  }
  combined_path <- file.path(counts_output_dir, "combined_star_gene_counts.csv")
  write.csv(combined_counts, combined_path, row.names = FALSE)
  cat("Combined gene count matrix saved to:", combined_path, "\n")
  filtered_counts <- filter_gene_counts(
    combined_counts,
    rrna_file,
    protein_coding_filter = isTRUE(protein_coding_filter),
    remove_mt_genes = isTRUE(remove_mt_genes),
    go_annotation_filter = isTRUE(go_annotation_filter)
  )
  if (!is.null(filtered_counts)) {
    filtered_path <- file.path(counts_output_dir, "filtered_gene_counts.csv")
    write.csv(filtered_counts, filtered_path, row.names = FALSE)
    cat("Filtered gene count matrix saved to:", filtered_path, "\n")
  }
  return(list(combined = combined_counts, filtered = filtered_counts))
}

# Main processing function
process_dataset <- function(dataset_name, kreports_folder, bracken_folder = NULL, 
                            process_bracken = PROCESS_BRACKEN, add_params = ADD_PARAMS, 
                            include_subspecies = INCLUDE_SUBSPECIES) {
  cat("\n=== Processing", dataset_name, "dataset ===\n")
  
  # Process kreports
  kreport_results <- process_kreport_folder(kreports_folder, include_subspecies = include_subspecies)
  
  # Process bracken if folder provided and processing is enabled
  bracken_data <- NULL
  if (process_bracken && !is.null(bracken_folder) && dir.exists(bracken_folder)) {
    bracken_data <- process_bracken_folder(bracken_folder)
    cat("Bracken processing enabled for", dataset_name, "\n")
  } else {
    cat("Bracken processing skipped for", dataset_name, "\n")
  }
  
  # Merge data
  merge_results <- merge_data(kreport_results$kreports_combined, bracken_data, 
                              add_params = add_params, report_data = kreport_results$report_data,
                              dataset_name = dataset_name)
  
  # Add summary statistics
  merged_with_stats <- add_summary_stats(merge_results$merged)
  
  # Get species list with topN frequency
  species_list <- get_species_list_with_frequency(merged_with_stats, top_n = TOP_SPECIES_N)
  cat("Generated species list with", nrow(species_list), "species\n")
  
  result <- list(
    report_data = kreport_results$report_data,
    merged = merged_with_stats,
    merged_long = merge_results$merged_long,
    species_list = species_list,
    bracken_data = bracken_data
  )
  return(result)
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

cat("Starting kraken/bracken output processing...\n")

# Check if unaligned kreports directory is specified
if (is.null(UNALIGNED_KREPORTS_DIR)) {
  cat("ERROR: No unaligned kreports directory specified. Use --unaligned-kreports to specify the path.\n")
  stop("Required unaligned kreports directory not specified.")
}

# Process unaligned dataset
unaligned_bracken_folder <- if (PROCESS_BRACKEN) UNALIGNED_BRACKEN_DIR else NULL
unaligned_results <- process_dataset(
  dataset_name = "Unaligned",
  kreports_folder = UNALIGNED_KREPORTS_DIR,
  bracken_folder = unaligned_bracken_folder,
  process_bracken = PROCESS_BRACKEN,
  add_params = ADD_PARAMS,
  include_subspecies = INCLUDE_SUBSPECIES
)

# Skip processing nonhuman dataset if NONHUMAN_KREPORTS_DIR is NULL
if (!is.null(NONHUMAN_KREPORTS_DIR)) {
  # Process nonhuman dataset  
  nonhuman_bracken_folder <- if (PROCESS_BRACKEN) NONHUMAN_BRACKEN_DIR else NULL
  nonhuman_results <- process_dataset(
    dataset_name = "Nonhuman",
    kreports_folder = NONHUMAN_KREPORTS_DIR, 
    bracken_folder = nonhuman_bracken_folder,
    process_bracken = PROCESS_BRACKEN,
    add_params = ADD_PARAMS,
    include_subspecies = INCLUDE_SUBSPECIES
  )
} else {
  cat("\nSkipping nonhuman dataset processing as NONHUMAN_KREPORTS_DIR is NULL\n")
}

# Only save if SAVE_OUTPUT is TRUE
if (SAVE_OUTPUT) {
  cat("\n=== Saving outputs ===\n")
  output_dir <- OUTPUT_DIR
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # Save unaligned results
  output_csv_file(unaligned_results$merged, "unaligned_merged", output_dir, "op")
  output_csv_file(unaligned_results$merged_long, "unaligned_merged_long", output_dir, "op")
  output_csv_file(unaligned_results$species_list, "species_list_unaligned", output_dir, "op")
  processed_output_dir <- here("data", "processed")
  if (!dir.exists(processed_output_dir)) {
    dir.create(processed_output_dir, recursive = TRUE)
  }
  write.csv(unaligned_results$merged, file.path(processed_output_dir, "unaligned_merged.csv"), row.names = FALSE)
  write.csv(unaligned_results$merged_long, file.path(processed_output_dir, "unaligned_merged_long.csv"), row.names = FALSE)
  write.csv(unaligned_results$species_list, file.path(processed_output_dir, "species_list_unaligned.csv"), row.names = FALSE)
  
  # Save unaligned bracken data if it exists
  if (!is.null(unaligned_results$bracken_data) && nrow(unaligned_results$bracken_data) > 0) {
    output_csv_file(unaligned_results$bracken_data, "unaligned_bracken", output_dir, "op")
  }
  
  # Save nonhuman results only if they exist
  if (!is.null(NONHUMAN_KREPORTS_DIR) && exists("nonhuman_results")) {
    output_csv_file(nonhuman_results$merged, "nonhuman_merged", output_dir, "op")
    output_csv_file(nonhuman_results$merged_long, "nonhuman_merged_long", output_dir, "op")
    output_csv_file(nonhuman_results$species_list, "species_list_nonhuman", output_dir, "op")
    write.csv(nonhuman_results$merged, file.path(processed_output_dir, "nonhuman_merged.csv"), row.names = FALSE)
    write.csv(nonhuman_results$merged_long, file.path(processed_output_dir, "nonhuman_merged_long.csv"), row.names = FALSE)
    write.csv(nonhuman_results$species_list, file.path(processed_output_dir, "species_list_nonhuman.csv"), row.names = FALSE)
    
    # Save nonhuman bracken data if it exists
    if (!is.null(nonhuman_results$bracken_data) && nrow(nonhuman_results$bracken_data) > 0) {
      output_csv_file(nonhuman_results$bracken_data, "nonhuman_bracken", output_dir, "op")
    }
  }
  # Save combined report data
  if (!is.null(NONHUMAN_KREPORTS_DIR) && exists("nonhuman_results")) {
    combined_report_data <- bind_rows(
      unaligned_results$report_data %>% mutate(dataset = "unaligned"),
      nonhuman_results$report_data %>% mutate(dataset = "nonhuman")
    )
  } else {
    combined_report_data <- unaligned_results$report_data %>% mutate(dataset = "unaligned")
  }
  output_csv_file(combined_report_data, "sample_report_data", output_dir, "op")
  
  # Parse and save metadata files if requested
  cat("\n=== Parsing additional metadata ===\n")
  
  # Parse runtime files if directory is provided
  if (!is.null(RUNTIME_DIR)) {
    cat("Parsing runtime files...\n")
    runtime_data <- parse_runtime_files(RUNTIME_DIR)
    if (!is.null(runtime_data)) {
      output_csv_file(runtime_data, "runtime_info", output_dir, "op")
      cat("Runtime data saved successfully\n")
    }
  }
  
  # Parse read statistics files if directory is provided
  if (!is.null(RRSTATS_DIR)) {
    cat("Parsing read statistics files...\n")
    rrstats_data <- parse_rrstats_files(RRSTATS_DIR)
    if (!is.null(rrstats_data)) {
      output_csv_file(rrstats_data, "read_statistics", output_dir, "op")
      cat("Read statistics data saved successfully\n")
    }
  }
  
  # Parse sample metadata files if directory is provided
  if (!is.null(METADATA_DIR)) {
    cat("Parsing sample metadata files...\n")
    metadata_data <- parse_sample_metadata(METADATA_DIR)
    if (!is.null(metadata_data)) {
      output_csv_file(metadata_data, "sample_metadata", output_dir, "op")
      cat("Sample metadata saved successfully\n")
    }
  }
  
  # If multiple metadata types were parsed, create a combined metadata file
  parsed_metadata <- list()
  if (!is.null(RUNTIME_DIR) && exists("runtime_data") && !is.null(runtime_data)) {
    parsed_metadata$runtime <- runtime_data
  }
  if (!is.null(RRSTATS_DIR) && exists("rrstats_data") && !is.null(rrstats_data)) {
    parsed_metadata$rrstats <- rrstats_data
  }
  if (!is.null(METADATA_DIR) && exists("metadata_data") && !is.null(metadata_data)) {
    parsed_metadata$metadata <- metadata_data
  }
  
  if (length(parsed_metadata) > 1) {
    cat("Creating combined metadata file...\n")
    # Start with the first dataset as base
    combined_metadata <- parsed_metadata[[1]]
    
    # Merge with other datasets
    for (i in 2:length(parsed_metadata)) {
      combined_metadata <- merge(combined_metadata, parsed_metadata[[i]], by = "sample", all = TRUE)
    }
    
    output_csv_file(combined_metadata, "combined_metadata", output_dir, "op")
    cat("Combined metadata saved successfully\n")
  }
  
  cat("Processing completed successfully!\n")
  cat("Output files saved to:", output_dir, "\n")
} else {
  cat("CSV output skipped (--output not specified)\n")
  cat("Processing completed successfully!\n")
}

# Process RPM count matrices and filter genes
if (PROCESS_COUNTS) {
  cat("\n=== Processing RPM count matrices ===\n")
  count_results <- process_gene_counts(
    RPM_DIR,
    COUNTS_OUTPUT_DIR,
    RRNA_FILE,
    protein_coding_filter = PROTEIN_CODING_FILTER,
    remove_mt_genes = REMOVE_MT_GENES,
    go_annotation_filter = GO_ANNOTATION_FILTER
  )
  if (is.null(count_results)) {
    cat("No RPM count matrices processed\n")
  }
}

# Clean up configuration variables
rm(args, PROCESS_BRACKEN, SAVE_OUTPUT, ADD_PARAMS, INCLUDE_SUBSPECIES, USE_PARALLEL, OUTPUT_DIR, MIN_CLADE_READS, RUNTIME_DIR, RRSTATS_DIR, METADATA_DIR, TOP_SPECIES_N, parsed_metadata, RPM_DIR, COUNTS_OUTPUT_DIR, RRNA_FILE, PROCESS_COUNTS, PROTEIN_CODING_FILTER, REMOVE_MT_GENES, GO_ANNOTATION_FILTER)
