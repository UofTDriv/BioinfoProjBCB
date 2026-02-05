# R analysis package for RNAseq pipeline with Kraken2 metatranscriptomic data.
Denis R | Jeff Wrana Lab

## Project Description

## Pipeline Dependencies
The following (along with all dependencies) need to be installed in the computer node image – contact cluster admin if you receive any “ModuleNotFoundError”.
1.	R (v4.0+)
  a.	R packages includes: `here`, `ggplot2`, `dplyr`, `tidyr`, `stringr`, `purrr`, `tibble`, `ggrepel`, `kableExtra`, `gridExtra`, `openxlsx`, `foreach`, `doParallel`
## Input Folder Structure
`[not_accessed_by_Ranalysis]` 
```
project-name_date/
├── pipeline_log
├── version_time_log
├── runtime
├── QC/
│   └── run_read_stats.txt
├── kraken2/
│   ├── c<confidence>m<min-hits><db-dir>_unaligned/
│   │   ├── <sampleid>_<params>_report.k2report
│   │   ├── <sampleid>_<params>_report_bracken_species.k2report
│   │   ├── output/
│   │   └── unclassified/
│   └── c<confidence>m<min-hits><db-dir>_nonhuman/
│       ├── <sampleid>_<params>_report.k2report
│       ├── <sampleid>_<params>_report_bracken_species.k2report
│       ├── output/
│       └── unclassified/
├── bracken/
│   ├── c<confidence>m<min-hits><db-dir>_unaligned/
│   │   └── <sampleid>_bracken_S.txt
│   └── c<confidence>m<min-hits><db-dir>_nonhuman/
│       └── <sampleid>_bracken_S.txt
├── Ranalysis/
│   ├── install.R
│   ├── scripts/
│   │   ├── utils.R
│   │   ├── output_processing.R
│   │   ├── annotate_species.R
│   │   ├── annotate_species_plots.R
│   │   └── reportRanalysis.R
│   └── databases/
│       ├── epathogen-\<timestamp>-result.csv
│       └── HOMD_taxon_table\<timestamp>.csv
├── [STAR]/
├── [BAM]/
├── [RPM]/
├── [jobs]/
├── [resources]/
└── [raw_fastqs]/
```

## Output Folder Structure
All outputs into Ranalysis/outputs

## R analyis
Outputs of the RNAseq pipeline are processed downstream in R. The analysis and report generation is orchestrated by the reportRanalysis.R script which uses multiple other R scripts as modules to carry out specific functions. These are executed by using source(“scripts”) in R; hence, all variables are shared via the GlobalEnvir in R.  
Nonetheless each script was built to be fully modular and usable independently (except using shared functions stored in utils.R).  
**reportRanalysis.R usage:**  
`Rscript reportRanalysis.R --input-dir <path> --proj-name <name> [options]`

### **a.**     **Install.R**
- Initiates a renv for improved package management in the project.
- Lists package installs for the packages used in the project.

### **b.**     **utils.R**
- This script contains functions that are used across multiple files. It contains the run_as_script function to run R scripts as standalone scripts within R (traditionally requiring using Rscript in the terminal).
- It also contains read_report2 and its helper functions: build_kraken_tree, collapse.taxRanks, delete_taxRanks_below which are used to ingest k2reports and are taken directly from the Pavian metagenomics viewer codebase: [https://github.com/fbreitwieser/pavian/blob/master/R/datainput-read_report.R](https://github.com/fbreitwieser/pavian/blob/master/R/datainput-read_report.R)
- Finally, it contains utility functions for standardized file I/O.

### **c.**     **output_processing.R**
- Processes Kraken2 and Bracken output files (kreport and bracken tables) from both unaligned and nonhuman datasets.
- Parses and filters kreport/bracken files, for samples supporting various naming conventions.
- Merges and summarizes species abundance data across samples.
- Supports filtering by clade reads, minimizer ratio, and exclusion of specific taxIDs.
- Optionally processes runtime and read statistics, and can run in parallel.
- Outputs merged data for downstream annotation and reporting.
- Will generate the following files:
	- `combined_metadata<date>op.csv`
		- Includes “sample” and the associated runtime, readrun stats, and sample metadata into a single csv file where each row is a sample.
	- `unaligned_bracken<date>op.csv & nonhuman_bracken<date>op.csv`
		- Includes (species) “name”, “taxID”, and a column for each sample’s Bracken Reads attained from “new_est_reads” column in Bracken report outputs.
	- `unaligned_merged_long<date>op.csv & nonhuman_merged_long<date>op.csv`
		- The pivoted nonhuman (host-removed) merged data and the included columns are "taxID","name","taxRank","taxLineage","sample","cladeReads","minimizers","distinct_minimizers",\[**"**bracken_reads" (if available)\]**,** and if params are requested it will have \["confidence_levels","minimum_hit_groups","database_used","human_reads”\].
	- `sample_report_data<date>op.csv`
		- Includes summary data for each sample’s kreport with the included columns "unclassified_counts","classified_counts",\["confidence","minimum_hits","database_used"\] if params requested, "species_count","filtered_species_count", \[“subspecies_count”, “filtered subspecies_count”, “total_species_count”, “filtered_total_species_count”\] if subspecies are requested, "dataset" whether the sample was unaligned or nonhuman.
	- `read_statistics<date>op.csv`
		- Identical to the run_read_stats.txt in QC except outputted to csv as a backup.
	- `unaligned_merged<date>op.csv & nonhuman_merged<date>op.csv`
		- For both datasets: each sample’s kreport ingest into a matrix format that was merged by "taxID","name","taxRank","taxLineage" of each species. Each sample will have 3 columns with the suffixes: `“_cladeReads”, “_minimizers”, “_distinct_minimizers”.`
	- `runtime_info<date>op.csv`
		- The ingested data from the pipeline runtime as a csv with the following columns "sample","cutadapt_star_samtools_rt","kraken2_unaligned_rt","kraken2_nonhuman_rt","is_first_k2_unaligned", "is_first_k2_nonhuman" a TRUE or FALSE for whether this sample was the first to be processed by the kraken script and therefore was involved in loading in the database.
	- `species_list_unaligned<date>op.csv & species_list_nonhuman<date>op.csv`
		- Includes the list of all species in a dataset (unaligned or nonhuman) with "name","taxID","taxRank","taxLineage" as well as summary information such as: “Freq” appearance of the sample in the “--top-n-freq” species per sample; as well as "cladeReads_mean","cladeReads_max","bracken_reads_mean","bracken_reads_max" across all samples in the unaligned dataset.
		- species_list_unaligned is used downstream for species annotation and batch analysis reporting.

### **d.**     **annotate_species.R**
- Annotates species abundance tables with risk group (epathogen) and HOMD (Human Oral Microbiome Database) classifications.
- Loads the latest epathogen and HOMD database files.
- Annotates each species by taxID or name with risk group (RG1–RG4, NotAnnotated).
- Adds HOMD body site and category (Pathogen, Opportunist, Microbiome, NotAnnotated).
- Designed for both command-line and programmatic use; works by taking in a species_list dataframe from the global environment and outputs an annotated dataframe.

### **e.**     **annotate_species_plots.R**
- Takes in an annotated dataframe from annotate_species.R and generates annotation summary and rank order plots for annotated species data.
- Creates bar charts and summary plots by risk group, HOMD category, and kingdom.
- Customizable color schemes for risk groups, HOMD categories, and kingdoms.
- Can be run interactively or as a script for batch plot generation.

### **f.**      **reportRanalysis.R**
- Orchestrates the full downstream analysis and reporting workflow for Kraken2/Bracken results.
- Parses command-line arguments to configure analysis options and input/output locations.
- Runs output_processing.R to generate merged abundance tables.
- Annotates species using annotate_species.R.
- Generates species summary statistics and rank order visualizations using annotate_species_plots.R.
- Generates summary statistics, correlation analyses, and visual reports (PDF/plots).
- Supports parallel processing and modular report generation (runtime, read stats, kraken stats, species annotation, pathogen detection).
- Handles both batch and per-sample reporting, with flexible output directory management.

## Troubleshooting
```
{Ranalysis} r_analysis_<jobID>.out:
- The project is out-of-sync -- use `renv::status()` for details.
Error in library(here) : there is no package called ‘here’
Execution halted{samtools}
```
Solution:  
R packages did not install properly on the galen cluster. 
Make sure install.R was ran beforehand.
1.	Switch to devhouse server: cd into the project directory’s Ranalysis folder. Start a new R session. 
```
# Bootstrapping renv 1.1.5 ---------------------------------------------------
- Downloading renv ... OK
- Installing renv  ... OK

- Project '~/RNASeqProject/SPARKEN/RNASeqPipeline/Ranalysis' loaded. [renv 1.1.5]
- The project is out-of-sync -- use `renv::status()` for details.
```
After the R welcome message, you should see the R project be loaded. 
\[Running renv::status() should tell you the “The following package(s) are used in this project, but are not installed:” followed by the R packages listed in dependencies\]
2.	Inside the R session run
renv::snapshot()
You will be told what packages weren’t installed followed by these options:
What do you want to do? 
1: Snapshot, just using the currently installed packages.
2: Install the packages, then snapshot.
3: Cancel, and resolve the situation on your own. 
3.	For Selection type 2 “Install the packages, then snapshot.” It will then tell you 
The following package(s) will be updated in the lockfile:
```
# CRAN -------------------------------------------------------------
<packages>
The version of R recorded in the lockfile will be updated:
- R              [4.4.1 -> 4.4.2]
Do you want to proceed? [Y/n]:
```
The version will be updated because there are different versions between Galen & devhouse. 
4.	Select Y and it will inform you that. 
`- Lockfile written to "<project-dir>/Ranalysis/renv.lock".`
You can then quit the R session without saving the workspace image.
5.	\[Using FileZilla\] inspect “`project-dir/Ranalysis/renv/library/`”, since the install was attempted on both servers, in there should be linux-rocky-9.4/R-4.4 and linux-rocky-9.5/R-4.4 for galen and devhouse respectively.
Move/replace R-4.4 from linux-rocky-9.5 to linux-rocky-9.4. 
6.	\[Start an R session inside /Ranalysis/ on galen and do renv::snapshot() to update the lockfile to galen’s R version\]

## References

HOMD database:
Isabel Fernández Escapa, Tsute Chen, Yanmei Huang, Prasad Gajare, Floyd E Dewhirst, and Katherine P Lemon. New insights into human nostril microbiome from the expanded Human Oral Microbiome Database (eHOMD): a resource for species-level identification of microbiome data from the aerodigestive tract. DOI: 10.1128/mSystems.00187-18.
Online Open Access: https://msystems.asm.org/content/3/6/e00187-18

Riskgroup database: 
Health Canada. ePathogen – Pathogen Risk Assessment Database. Government of Canada.
Online Open Access: https://health.canada.ca/en/epathogen
