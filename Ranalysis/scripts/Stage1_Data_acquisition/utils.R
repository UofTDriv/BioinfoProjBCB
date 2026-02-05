library(here)
library(stringr)
library(ggplot2)

output_csv_file <- function(x, file_name, dir, script_name) {
  timestamp <- format(Sys.time(), "%y%m%d")
  # Construct the full file path
  file_path <- file.path(here(dir), paste0(file_name, timestamp, script_name, ".csv"))
  write.csv(x = x, file = file_path, row.names = F)
}

output_sample_set <- function(x, script_name, file_name = "Uploaded_sample_set-matrix-all-", dir = "inputs") {
  timestamp <- format(Sys.time(), "%y%m%d")
  # Construct the full file path
  file_path <- file.path(here(dir), paste0(file_name, timestamp, script_name, ".tsv"))
  write.table(x = x, file = file_path, sep = "\t", row.names = F)
}

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

## GPTo3 written function
build_lineage <- function(df, taxRankColName, lvl_regex = "^S\\d*$") {
  ## how many leading spaces precede each name?
  leading  <- nchar(df$name) - nchar(str_trim(df$name, side = "left"))
  ## indentation width (usually 2 spaces); fall back to 1 if none found
  indent   <- min(leading[leading > 0], default = 1)
  depth    <- leading / indent           # 0 = root, 1 = kingdom, …
  
  lin_vec  <- character(nrow(df))
  stack    <- character(max(depth) + 1)  # running “path” buffer
  
  for (i in seq_len(nrow(df))) {
    d <- depth[i] + 1                    # R is 1-based
    stack[d] <- str_trim(df$name[i])     # replace/extend current level
    if (length(stack) > d) stack <- stack[seq_len(d)]      # discard deeper levels
    if (grepl(lvl_regex, df[[taxRankColName]][i])) {
      lin_vec[i] <- paste(stack[seq_len(d)], collapse = ">")
    }
  }
  # TODO: Fix the function to not need to remove the root and NA
  ret <- gsub("NA>", "", lin_vec, fixed = TRUE)
  ret <- gsub("root>", "", ret, fixed = TRUE)
  ret
}

run_as_script <- function(script_path, ...) {
  # Simulate command-line arguments
  args <- c(...)
  assign("commandArgs", function(trailingOnly = TRUE) args, envir = .GlobalEnv)
  # Run the script
  source(script_path)
}

my_ggplot_theme <- theme_bw() + 
  theme(legend.position = "none", 
        plot.background = element_blank(), 
        axis.text.x = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black", 
                                    fill=NA, 
                                    linewidth=1),
        axis.text.y = element_text(size = 12))


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
      dmessage(sprintf("removed %7s cladeReads, including %s childs, for %s",child_rows[1,'"cladeReads"'],nrow(child_rows)-1,child_rows[1,'name']))
      
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
    dmessage(myfile," is no valid report - not all characters are ASCII")
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

