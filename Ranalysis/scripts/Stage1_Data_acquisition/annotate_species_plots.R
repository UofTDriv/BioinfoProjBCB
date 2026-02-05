#!/usr/bin/env Rscript

# Species Annotation Plotting Script
# Author: Denis Rivard, modified by GitHub Copilot
# Description: Generates various plots for annotated species data

# Load required libraries
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
})

# Define a standard ggplot theme
my_ggplot_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    strip.text = element_text(size = 12, face = "bold")
  )

# Define color scheme for Risk Groups
category_colors <- c(
  "RG1" = "#4CAF50",      # Green - Low risk
  "RG2" = "#FF9800",      # Orange - Moderate risk  
  "RG3" = "#E53935",      # Red - High risk
  "RG4" = "#B71C1C",      # Dark red - Very high risk
  "NotAnnotated" = "#9E9E9E"   # Gray - Not annotated
)

#' Create Risk Group Bar Chart
#' 
#' @param df Input dataframe with RiskGroup and other annotation columns
#' @param sample_name Custom sample name for title
#' @return ggplot object
create_risk_group_plot <- function(df, sample_name = "All Samples") {
  
  # Count species by risk group
  rg_levels <- c("RG1", "RG2", "RG3", "RG4", "NotAnnotated")
  
  rg_counts <- df %>%
    count(RiskGroup) %>%
    mutate(RiskGroup = factor(RiskGroup, levels = rg_levels)) %>%
    complete(RiskGroup, fill = list(n = 0)) %>%
    arrange(RiskGroup)
  
  # Get example species for each risk group (top 3)
  example_species <- df %>%
    group_by(RiskGroup) %>%
    slice_head(n = 3) %>%
    summarise(examples = paste(name, collapse = "\n"), .groups = 'drop')
  
  # Merge with counts
  rg_counts <- left_join(rg_counts, example_species, by = "RiskGroup")
  
  # Create plot
  rg_plot <- ggplot(rg_counts, aes(x = RiskGroup, y = n, fill = RiskGroup)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = category_colors) +
    labs(
      title = paste("Species per Risk Group -", sample_name),
      subtitle = paste("Total species:", sum(rg_counts$n)),
      x = "Risk Group",
      y = "Number of Species"
    ) +
    my_ggplot_theme +
    theme(
      axis.text.x = element_text(angle = 0),
      legend.position = "none",
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    ) +
    geom_text(data = subset(rg_counts, RiskGroup %in% c("RG1", "RG2", "NotAnnotated") & !is.na(examples)),
              aes(label = examples, y = n/2), color = "black", size = 2.5, fontface = "italic", lineheight = 0.9) +
    geom_text(data = subset(rg_counts, RiskGroup == "RG3" & !is.na(examples)),
              aes(label = examples, y = n + 0.05*max(rg_counts$n)), color = "black", size = 2.5, fontface = "italic", lineheight = 0.9)
  
  return(rg_plot)
}

#' Create HOMD Category Bar Chart
#' 
#' @param df Input dataframe with HOMD.Category and other annotation columns
#' @param sample_name Custom sample name for title
#' @return ggplot object
create_homd_plot <- function(df, sample_name = "All Samples") {
  
  # Count species by HOMD category
  homd_levels <- c("Pathogen", "Opportunist", "Microbiome", "NotAnnotated")
  
  homd_counts <- df %>%
    count(HOMD.Category) %>%
    mutate(HOMD.Category = factor(HOMD.Category, levels = homd_levels)) %>%
    complete(HOMD.Category, fill = list(n = 0)) %>%
    arrange(HOMD.Category)
  
  # Define colors for HOMD categories
  homd_colors <- c(
    "Pathogen" = "#E53935",      # Red
    "Opportunist" = "#FF9800",   # Orange
    "Microbiome" = "#4CAF50",    # Green
    "NotAnnotated" = "#9E9E9E"   # Gray
  )
  
  # Get example species for each category
  example_species <- df %>%
    group_by(HOMD.Category) %>%
    slice_head(n = 3) %>%
    summarise(examples = paste(name, collapse = "\n"), .groups = 'drop')
  
  # Merge with counts
  homd_counts <- left_join(homd_counts, example_species, by = "HOMD.Category")
  
  # Create plot
  homd_plot <- ggplot(homd_counts, aes(x = HOMD.Category, y = n, fill = HOMD.Category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = homd_colors) +
    labs(
      title = paste("Species per HOMD Category -", sample_name),
      subtitle = paste("Total species:", sum(homd_counts$n)),
      x = "HOMD Category",
      y = "Number of Species"
    ) +
    my_ggplot_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    ) +
    geom_text(data = subset(homd_counts, !is.na(examples)),
              aes(label = examples, y = n/2), color = "black", size = 2.5, fontface = "italic", lineheight = 0.9)
  
  return(homd_plot)
}

#' Create Kingdom Distribution Plot
#' 
#' @param df Input dataframe with taxLineage column
#' @param sample_name Custom sample name for title
#' @return ggplot object
create_kingdom_plot <- function(df, sample_name = "All Samples") {
  
  # Extract kingdom from taxLineage - only consider entries that start with "r_root|k_"
  df_with_kingdom <- df %>%
  mutate(
    kingdom = case_when(
    # Extract kingdom from taxLineage only if it starts with r_root|k_
    str_detect(taxLineage, "^r_root\\|k_") ~ str_extract(taxLineage, "(?<=r_root\\|k_)[^|]+"),
    # Group all other entries as "Other non-kingdom roots"
    TRUE ~ "Other non-kingdom roots"
    )
  )
  
  # Count species by kingdom
  kingdom_counts <- df_with_kingdom %>%
  count(kingdom) %>%
  arrange(desc(n))
  
  # Define color palette for kingdoms
  kingdom_color_palette <- c(
  "#2E8B57",  # Sea Green
  "#FF6347",  # Tomato
  "#4169E1",  # Royal Blue
  "#DAA520",  # Goldenrod
  "#8B4513",  # Saddle Brown
  "#9370DB",  # Medium Purple
  "#20B2AA",  # Light Sea Green
  "#CD853F",  # Peru
  "#B22222",  # Fire Brick
  "#4682B4",  # Steel Blue
  "#808080"   # Gray for "Other non-kingdom roots"
  )
  
  # Create color mapping for kingdoms
  unique_kingdoms <- unique(kingdom_counts$kingdom)
  # Put "Other non-kingdom roots" at the end with gray color
  if ("Other non-kingdom roots" %in% unique_kingdoms) {
  other_kingdoms <- unique_kingdoms[unique_kingdoms != "Other non-kingdom roots"]
  ordered_kingdoms <- c(other_kingdoms, "Other non-kingdom roots")
  kingdom_colors <- setNames(kingdom_color_palette[seq_along(ordered_kingdoms)], ordered_kingdoms)
  } else {
  kingdom_colors <- setNames(kingdom_color_palette[seq_along(unique_kingdoms)], unique_kingdoms)
  }
  
  # Get example species for each kingdom
  example_species <- df_with_kingdom %>%
  group_by(kingdom) %>%
  slice_head(n = 3) %>%
  summarise(examples = paste(name, collapse = "\n"), .groups = 'drop')
  
  # Merge with counts
  kingdom_counts <- left_join(kingdom_counts, example_species, by = "kingdom")
  
  # Create plot
  kingdom_plot <- ggplot(kingdom_counts, aes(x = reorder(kingdom, n), y = n, fill = kingdom)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = kingdom_colors, na.value = "gray") +
  labs(
    title = paste("Species per Kingdom -", sample_name),
    subtitle = paste("Total species:", sum(kingdom_counts$n)),
    x = "Kingdom",
    y = "Number of Species"
  ) +
  my_ggplot_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  ) +
  geom_text(data = kingdom_counts,
        aes(label = paste0("n=", n), y = n + 0.02*max(n)), 
        color = "black", size = 3, fontface = "bold") +
  coord_flip()
  
  return(kingdom_plot)
}

#' Create Combined Risk Group and HOMD Summary Plot
#' 
#' @param df Input dataframe with annotation columns
#' @param sample_name Custom sample name for title
#' @return ggplot object
create_summary_plot <- function(df, sample_name = "All Samples") {
  
  # Create summary data
  summary_data <- df %>%
    summarise(
      total_species = n(),
      risk_annotated = sum(RiskGroup != "NotAnnotated"),
      homd_annotated = sum(HOMD.Category != "NotAnnotated"),
      high_risk = sum(RiskGroup %in% c("RG3", "RG4")),
      pathogens = sum(HOMD.Category %in% c("Pathogen", "Opportunist"))
    ) %>%
    pivot_longer(everything(), names_to = "category", values_to = "count") %>%
    mutate(
      category = factor(category, levels = c("total_species", "risk_annotated", "homd_annotated", "high_risk", "pathogens")),
      percentage = round(count / max(count) * 100, 1)
    )
  
  # Create plot
  summary_plot <- ggplot(summary_data, aes(x = category, y = count, fill = category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c(
      "total_species" = "#2196F3",
      "risk_annotated" = "#4CAF50", 
      "homd_annotated" = "#FF9800",
      "high_risk" = "#E53935",
      "pathogens" = "#9C27B0"
    )) +
    labs(
      title = paste("Species Annotation Summary -", sample_name),
      x = "Category",
      y = "Number of Species"
    ) +
    my_ggplot_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    ) +
    geom_text(aes(label = paste0(count, "\n(", percentage, "%)")), 
              vjust = -0.3, color = "black", size = 3)
  
  return(summary_plot)
}

#' Create Detailed Species Plot
#' 
#' @param df Input dataframe with annotation and value columns
#' @param detailed_col Column name containing values (e.g., Frequency, CladeReads)
#' @param plot_type Plot type to determine category column and colors
#' @param sample_name Custom sample name for title
#' @param top_n Number of top species to show (default: 50)
#' @param log Whether to create log-transformed plots (default: FALSE)
#' @return ggplot object
create_detailed_plot <- function(df, detailed_col, plot_type, sample_name = "All Samples", top_n = 50, log = FALSE) {
  
  # Check if detailed column exists
  if (!detailed_col %in% colnames(df)) {
  stop(paste("Column", detailed_col, "not found in dataframe. Available columns:", 
         paste(colnames(df), collapse = ", ")))
  }
  
  # Determine category column and colors based on plot_type
  if (plot_type == "risk_group") {
  category_col <- "RiskGroup"
  plot_colors <- category_colors
  plot_title <- paste("Species by Risk Group -", sample_name)
  } else if (plot_type == "homd") {
  category_col <- "HOMD.Category"
  plot_colors <- c(
    "Pathogen" = "#E53935",      # Red
    "Opportunist" = "#FF9800",   # Orange
    "Microbiome" = "#4CAF50",    # Green
    "NotAnnotated" = "#9E9E9E"   # Gray
  )
  plot_title <- paste("Species by HOMD Category -", sample_name)
  } else if (plot_type == "kingdom") {
  # Extract kingdom for kingdom plots - only consider entries that start with "r_root|k_"
  df <- df %>%
    mutate(kingdom = case_when(
    # Extract kingdom from taxLineage only if it starts with r_root|k_
    str_detect(taxLineage, "^r_root\\|k_") ~ str_extract(taxLineage, "(?<=r_root\\|k_)[^|]+"),
    # Group all other entries as "Other non-kingdom roots"
    TRUE ~ "Other non-kingdom roots"
    ))
  category_col <- "kingdom"
  
  # Create dynamic colors for kingdoms
  unique_kingdoms <- unique(df$kingdom)
  kingdom_color_palette <- c(
    "#2E8B57", "#FF6347", "#4169E1", "#DAA520", "#8B4513",
    "#9370DB", "#20B2AA", "#CD853F", "#B22222", "#4682B4", "#808080"
  )
  # Put "Other non-kingdom roots" at the end with gray color
  if ("Other non-kingdom roots" %in% unique_kingdoms) {
    other_kingdoms <- unique_kingdoms[unique_kingdoms != "Other non-kingdom roots"]
    ordered_kingdoms <- c(other_kingdoms, "Other non-kingdom roots")
    plot_colors <- setNames(kingdom_color_palette[seq_along(ordered_kingdoms)], ordered_kingdoms)
  } else {
    plot_colors <- setNames(kingdom_color_palette[seq_along(unique_kingdoms)], unique_kingdoms)
  }
  plot_title <- paste("Species by Kingdom -", sample_name)
  } else {
  stop("Detailed plots are only available for plot_type: risk_group, homd, or kingdom")
  }
  
  # Check if category column exists
  if (!category_col %in% colnames(df)) {
  stop(paste("Category column", category_col, "not found in dataframe for plot_type", plot_type))
  }
  
  # Prepare data: sort by detailed column and take top_n
  plot_data <- df %>%
  arrange(desc(get(detailed_col))) %>%
  slice_head(n = top_n) %>%
  mutate(
    name_f = factor(name, levels = name),  # Keep descending order (largest to smallest left to right)
    category = factor(get(category_col), levels = names(plot_colors))
  )
  
  # Apply log transformation only if requested
  if (log) {
  plot_data[[detailed_col]] <- log10(plot_data[[detailed_col]] + 1)
  y_label <- paste0("Log10(", detailed_col, " + 1)")
  } else {
  y_label <- detailed_col
  }
  
  detailed_plot <- ggplot(plot_data, aes(x = name_f, y = get(detailed_col), fill = category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Category", values = plot_colors, na.value = "black") +
  labs(
    title = plot_title,
    subtitle = paste("Top", min(top_n, nrow(plot_data)), "species by", y_label),
    x = "Species",
    y = y_label
  ) +
  my_ggplot_theme +
  theme(
    legend.position = "right",
    plot.background = element_blank(),
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1, size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.y = element_text(size = 10, color = "black")
  )
  
  return(detailed_plot)
}

#' Parse command line arguments
#' 
#' @return List of parsed arguments
parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Default configuration
  config <- list(
    df = NULL,
    plot_type = "risk_group",
    detailed = NULL,
    sample_name = "All Samples",
    output = NULL,
    width = 12,
    height = 8,
    help = FALSE,
    top_n = 30,
    log = FALSE
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
      } else if (arg == "--plot-type") {
        if (i + 1 <= length(args)) {
          i <- i + 1
          config$plot_type <- args[i]
        } else {
          stop("--plot-type requires an argument")
        }
      } else if (arg == "--detailed") {
        if (i + 1 <= length(args)) {
          i <- i + 1
          config$detailed <- args[i]
        } else {
          stop("--detailed requires an argument")
        }
      } else if (arg == "--sample-name") {
        if (i + 1 <= length(args)) {
          i <- i + 1
          config$sample_name <- args[i]
        } else {
          stop("--sample-name requires an argument")
        }
      } else if (arg == "--output") {
        if (i + 1 <= length(args)) {
          i <- i + 1
          config$output <- args[i]
        } else {
          stop("--output requires an argument")
        }
      } else if (arg == "--width") {
        if (i + 1 <= length(args)) {
          i <- i + 1
          config$width <- as.numeric(args[i])
        } else {
          stop("--width requires an argument")
        }
      } else if (arg == "--height") {
        if (i + 1 <= length(args)) {
          i <- i + 1
          config$height <- as.numeric(args[i])
        } else {
          stop("--height requires an argument")
        }
      } else if (arg == "--top-n") {
        if (i + 1 <= length(args)) {
          i <- i + 1
          config$top_n <- as.numeric(args[i])
        } else {
          stop("--top-n requires an argument")
        }
      } else if (arg == "--log") {
        config$log <- TRUE
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
  cat("Usage: Rscript annotate_species_plots.R [OPTIONS]\n\n")
  cat("Options:\n")
  cat("  --df <object>             R object name containing the annotated dataframe (required for standard plots)\n")
  cat("                            Must have RiskGroup, HOMD.Category, and other annotation columns\n")
  cat("  --plot-type <type>        Type of plot to generate (default: risk_group)\n")
  cat("                            Options: risk_group, homd, kingdom, summary\n")
  cat("  --detailed <column>       Create detailed species plot using specified column values\n")
  cat("                            (e.g., Frequency, CladeReads) - overrides standard plot\n")
  cat("  --sample-name <name>      Custom sample name for plot title (default: 'All Samples')\n")
  cat("  --output <file>           Output PDF file path (default: display in R environment)\n")
  cat("  --width <number>          Plot width in inches (default: 12)\n")
  cat("  --height <number>         Plot height in inches (default: 8)\n")
  cat("  --top-n <number>          Number of top species to show (default: 30)\n")
  cat("  --log                     Create log-transformed plots (default: FALSE)\n")
  cat("  --help, -h               Show this help message\n\n")
  cat("Plot Types:\n")
  cat("  risk_group:              Bar chart showing species counts by Risk Group (RG1-RG4)\n")
  cat("  homd:                    Bar chart showing species counts by HOMD Category\n")
  cat("  kingdom:                 Bar chart showing species counts by taxonomic kingdom\n")
  cat("  summary:                 Summary plot with annotation statistics\n\n")
  cat("Detailed Plot:\n")
  cat("  When --detailed is used, creates a scatter plot showing individual species\n")
  cat("  ranked by the specified column value and colored by category\n")
  cat("  (similar to reportRanalysis.R style)\n\n")
  cat("Input Requirements:\n")
  cat("  The input dataframe must contain:\n")
  cat("  - 'name': Species name column\n")
  cat("  - 'taxID': Taxonomy ID column\n")
  cat("  - 'RiskGroup': Risk group annotation from annotate_species.R\n")
  cat("  - 'HOMD.Category': HOMD category annotation from annotate_species.R\n")
  cat("  - 'taxLineage': Taxonomic lineage (for kingdom plots)\n\n")
  cat("Examples:\n")
  cat("  run_as_script(here('scripts', 'annotate_species_plots.R'), '--df', 'annotated_species')\n")
  cat("  run_as_script(here('scripts', 'annotate_species_plots.R'), '--df', 'my_data', '--plot-type', 'homd')\n")
  cat("  run_as_script(here('scripts', 'annotate_species_plots.R'), '--df', 'my_data', '--plot-type', 'kingdom')\n")
  cat("  run_as_script(here('scripts', 'annotate_species_plots.R'), '--df', 'my_data', '--plot-type', 'homd', '--detailed', 'Frequency')\n\n")
}

#' Main function to orchestrate plot generation
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
    cat("Error: --df argument is required for plots.\n\n")
    show_help()
    stop("Missing required arguments")
  }
  
  # Validate plot type
  valid_plot_types <- c("risk_group", "homd", "kingdom", "summary")
  if (!args$plot_type %in% valid_plot_types) {
    cat("Error: Invalid plot type. Must be one of:", paste(valid_plot_types, collapse = ", "), "\n\n")
    show_help()
    stop("Invalid plot type")
  }
  
  # Get dataframe from environment
  cat("Loading input dataframe from object:", args$df, "\n")
  if (!exists(args$df, envir = .GlobalEnv)) {
    stop(paste("Object", args$df, "not found in environment"))
  }
  df <- get(args$df, envir = .GlobalEnv)
  
  # Validate required columns based on plot type
  required_base_cols <- c("name", "taxID")
  
  # Check if detailed plot is requested
  if (!is.null(args$detailed)) {
    # For detailed plots, we need the appropriate category column
    if (args$plot_type == "risk_group") {
      required_cols <- c(required_base_cols, "RiskGroup")
    } else if (args$plot_type == "homd") {
      required_cols <- c(required_base_cols, "HOMD.Category")
    } else if (args$plot_type == "kingdom") {
      required_cols <- c(required_base_cols, "taxLineage")
    } else if (args$plot_type == "summary") {
      cat("Warning: Detailed plots are not available for summary plot type. Using standard summary plot.\n")
      required_cols <- c(required_base_cols, "RiskGroup", "HOMD.Category")
    }
  } else {
    # Standard plots
    if (args$plot_type == "risk_group") {
      required_cols <- c(required_base_cols, "RiskGroup")
    } else if (args$plot_type == "homd") {
      required_cols <- c(required_base_cols, "HOMD.Category")
    } else if (args$plot_type == "kingdom") {
      required_cols <- c(required_base_cols, "taxLineage")
    } else if (args$plot_type == "summary") {
      required_cols <- c(required_base_cols, "RiskGroup", "HOMD.Category")
    }
  }
  
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns for", args$plot_type, "plot:", paste(missing_cols, collapse = ", ")))
  }
  
  # Generate plot based on type and detailed flag
  if (!is.null(args$detailed) && args$plot_type != "summary") {
    cat("Generating detailed", args$plot_type, "plot using column:", args$detailed, "\n")
    plot_result <- create_detailed_plot(df, args$detailed, args$plot_type, args$sample_name, log = args$log)
  } else {
    cat("Generating", args$plot_type, "plot...\n")
    
    if (args$plot_type == "risk_group") {
      plot_result <- create_risk_group_plot(df, args$sample_name)
    } else if (args$plot_type == "homd") {
      plot_result <- create_homd_plot(df, args$sample_name)
    } else if (args$plot_type == "kingdom") {
      plot_result <- create_kingdom_plot(df, args$sample_name)
    } else if (args$plot_type == "summary") {
      plot_result <- create_summary_plot(df, args$sample_name)
    }
  }
  
  # Output or display
  if (!is.null(args$output)) {
    cat(paste("Saving plot to:", args$output, "\n"))
    ggsave(
      filename = args$output,
      plot = plot_result,
      width = args$width,
      height = args$height,
      device = "pdf"
    )
    cat("Plot saved successfully!\n")
  } else {
    cat("Storing plot in global environment...\n")
  }
  
  # Assign to global environment for access after run_as_script
  if (!is.null(args$detailed) && args$plot_type != "summary") {
    if(!is.null(args$sample_name) && args$sample_name != "All Samples") {
      plot_var_name <- paste0("detailed_", args$plot_type, "_plot_", args$sample_name)
    } else {
      plot_var_name <- paste0("detailed_", args$plot_type, "_plot")
    }
  } else {
    if(!is.null(args$sample_name) && args$sample_name != "All Samples") {
      plot_var_name <- paste0(args$plot_type, "_plot_", args$sample_name)
    } else {
       plot_var_name <- paste0(args$plot_type, "_plot")
    }
  }
  assign(plot_var_name, plot_result, envir = .GlobalEnv)
  cat("Plot stored as:", plot_var_name, "\n")
  
  return(plot_result)
}

#' Function-based plotting interface (recommended for programmatic use)
#' 
#' @param df Input annotated dataframe
#' @param plot_type Type of plot to generate
#' @param detailed Optional column name for detailed plot
#' @param sample_name Optional sample name for plot titles
#' @param top_n Number of top species to show (default: 30)
#' @param log Whether to create log-transformed plots (default: FALSE)
#' @return ggplot object
generate_annotation_plot <- function(df, plot_type = "risk_group", detailed = NULL, 
                                   sample_name = "All Samples", top_n = 30, log = FALSE) {
  
  if (is.null(df)) {
    stop("df argument is required")
  }
  
  # Validate plot type
  valid_plot_types <- c("risk_group", "homd", "kingdom", "summary")
  if (!plot_type %in% valid_plot_types) {
    stop("Invalid plot type. Must be one of:", paste(valid_plot_types, collapse = ", "))
  }
  
  # Validate required columns based on plot type
  required_base_cols <- c("name", "taxID")
  
  # Check if detailed plot is requested
  if (!is.null(detailed)) {
    # For detailed plots, we need the appropriate category column
    if (plot_type == "risk_group") {
      required_cols <- c(required_base_cols, "RiskGroup")
    } else if (plot_type == "homd") {
      required_cols <- c(required_base_cols, "HOMD.Category")
    } else if (plot_type == "kingdom") {
      required_cols <- c(required_base_cols, "taxLineage")
    } else if (plot_type == "summary") {
      warning("Detailed plots are not available for summary plot type. Using standard summary plot.")
      required_cols <- c(required_base_cols, "RiskGroup", "HOMD.Category")
    }
  } else {
    # Standard plots
    if (plot_type == "risk_group") {
      required_cols <- c(required_base_cols, "RiskGroup")
    } else if (plot_type == "homd") {
      required_cols <- c(required_base_cols, "HOMD.Category")
    } else if (plot_type == "kingdom") {
      required_cols <- c(required_base_cols, "taxLineage")
    } else if (plot_type == "summary") {
      required_cols <- c(required_base_cols, "RiskGroup", "HOMD.Category")
    }
  }
  
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns for", plot_type, "plot:", paste(missing_cols, collapse = ", ")))
  }
  
  # Generate plot based on type and detailed flag
  if (!is.null(detailed) && plot_type != "summary") {
    plot_result <- create_detailed_plot(df, detailed, plot_type, sample_name, top_n, log)
  } else {
    if (plot_type == "risk_group") {
      plot_result <- create_risk_group_plot(df, sample_name)
    } else if (plot_type == "homd") {
      plot_result <- create_homd_plot(df, sample_name)
    } else if (plot_type == "kingdom") {
      plot_result <- create_kingdom_plot(df, sample_name)
    } else if (plot_type == "summary") {
      plot_result <- create_summary_plot(df, sample_name)
    }
  }
  
  return(plot_result)
}

#' Example usage function for backward compatibility with run_as_script
#' 
#' @param df Input annotated dataframe
#' @param plot_type Type of plot to generate
#' @param detailed Optional column name for detailed plot
#' @return ggplot object or list of ggplot objects
example_usage <- function(df = NULL, plot_type = "risk_group", detailed = NULL) {
  
  if (is.null(df)) {
    stop("Standard mode requires df argument")
  }
  
  # Store object in global environment for run_as_script to access
  assign("my_annotated_data", df, envir = .GlobalEnv)
    
  # Build arguments
  args_list <- c(
    "--df", "my_annotated_data",
    "--plot-type", plot_type
  )
  
  if (!is.null(detailed)) {
    args_list <- c(args_list, "--detailed", detailed)
  }
  
  # Use run_as_script to call the plotting script
  run_as_script(
    here("scripts", "annotate_species_plots.R"),
    args_list
  )
  
  # Return the generated plot
  plot_var_name <- paste0(plot_type, "_plot")
  if (!is.null(detailed)) {
    plot_var_name <- paste0("detailed_", plot_type, "_plot")
  }
  return(get(plot_var_name, envir = .GlobalEnv))
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# Only execute main function when script is run directly, not when sourced
if (!exists("..main_execution_disabled..")) {
  # Check for help request first
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args || "-h" %in% args) {
    show_help()
    stop("Help requested. Execution stopped.")
  }
  
  # Run main function
  main()
}