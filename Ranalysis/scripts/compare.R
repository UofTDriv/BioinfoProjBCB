# Define the directories
dir1 <- "data/covid_rerun_STAR"
dir2 <- "data/covid2021_STAR"

# Get the list of files
files1 <- list.files(dir1, full.names = TRUE)
files2 <- list.files(dir2, full.names = TRUE)

# Ensure the file lists are the same
if (!identical(basename(files1), basename(files2))) {
  stop("File lists are not identical.")
}

# Compare each pair of files
for (i in seq_along(files1)) {
  file1 <- files1[i]
  file2 <- files2[i]

  cat("Comparing", basename(file1), "and", basename(file2), "\n")

  # Read the files as data frames, using the first row as header and the first column as row names
  df1 <- read.table(file1, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
  df2 <- read.table(file2, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

  # Compare the data frames
  if (identical(df1, df2)) {
    cat("  Files are identical.\n")
  } else {
    cat("  Files are different.\n")

    # Compare dimensions
    if (nrow(df1) != nrow(df2)) {
      cat("    - Different number of rows:", nrow(df1), "vs", nrow(df2), "\n")
    }
    if (ncol(df1) != ncol(df2)) {
      cat("    - Different number of columns:", ncol(df1), "vs", ncol(df2), "\n")
    }

    # If dimensions are the same, compare further
    if (nrow(df1) == nrow(df2) && ncol(df1) == ncol(df2)) {
      # Compare column names
      if (!identical(colnames(df1), colnames(df2))) {
        cat("    - Column names are different.\n")
        # Find and print differing column names
        diff_cols <- setdiff(colnames(df1), colnames(df2))
        if(length(diff_cols) > 0) cat("      Columns in first file but not in second:", paste(diff_cols, collapse=", "), "\n")
        diff_cols <- setdiff(colnames(df2), colnames(df1))
        if(length(diff_cols) > 0) cat("      Columns in second file but not in first:", paste(diff_cols, collapse=", "), "\n")
      }

      # Compare row names
      if (!identical(rownames(df1), rownames(df2))) {
        cat("    - Row names are different.\n")
        # Find and print differing row names
        diff_rows <- setdiff(rownames(df1), rownames(df2))
        if(length(diff_rows) > 0) cat("      Rows in first file but not in second:", paste(diff_rows, collapse=", "), "\n")
        diff_rows <- setdiff(rownames(df2), rownames(df1))
        if(length(diff_rows) > 0) cat("      Rows in second file but not in first:", paste(diff_rows, collapse=", "), "\n")
      }
      
      # Use all.equal for a detailed comparison of the data
      # Capture the output of all.equal to prevent it from printing to the console directly
      comparison <- all.equal(df1, df2, check.attributes = FALSE)
      if (!isTRUE(comparison)) {
        cat("    - Data differences found:\n")
        
        # To provide a more specific summary of where the differences are
        # we can check which cells are different
        diff_matrix <- df1 != df2
        
        # This can be memory intensive for large dataframes, so we will only summarize
        num_diffs <- sum(diff_matrix, na.rm = TRUE)
        cat("      - Found", num_diffs, "differing data points.\n")
        
        # Optionally, show the first few differences
        if (num_diffs > 0) {
          cat("      - Some examples of differences:\n")
          
          # Find indices of first few differences
          diff_indices <- which(diff_matrix, arr.ind = TRUE)
          
          # Print info about the first 5 differences
          for(j in 1:min(5, nrow(diff_indices))) {
            row_idx <- diff_indices[j, 1]
            col_idx <- diff_indices[j, 2]
            
            gene <- rownames(df1)[row_idx]
            sample <- colnames(df1)[col_idx]
            
            val1 <- df1[row_idx, col_idx]
            val2 <- df2[row_idx, col_idx]
            
            cat("        - Gene '", gene, "', Sample '", sample, "': ", val1, " vs ", val2, "\n", sep="")
          }
        }
      }
    }
  }
}