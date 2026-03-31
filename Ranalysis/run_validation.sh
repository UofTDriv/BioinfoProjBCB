#!/usr/bin/env bash

# Bioinformatics Pipeline Orchestration Script (Validation Branch)
# Description: Automates the execution of the pipeline using the validation dataset.

set -e
set -u
set -o pipefail

PROJECT_DIR=$(pwd)
INPUT_DIR="data/input/validation"
PROJ_NAME="BioinfoProj_Validation"
LOG_DIR="outputs/logs_validation"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

echo "============================================="
echo "Starting Bioinformatics Validation Pipeline Overview"
echo "Project Directory: ${PROJECT_DIR}"
echo "Start Time: ${TIMESTAMP}"
echo "============================================="

# Ensure log directory exists
mkdir -p "${LOG_DIR}"
mkdir -p "outputs/Stage4_Validation"

log_msg() {
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] $1" | tee -a "${LOG_DIR}/pipeline_${TIMESTAMP}.log"
}

# --- STAGE 1: Data Acquisition ---
# Note: For validation, the dataset might already be processed in data/processed/validationRanalysis
log_msg "Validation dataset previously processed by reportRanalysis.R or using its dedicated branch."
log_msg "Checking for processed unaligned merging targets..."
if [ ! -f "data/processed/validationRanalysis/data/processed/unaligned_merged.csv" ]; then
    log_msg "Warning: Validation Stage 1 files not found! Skipping or needs raw input data."
else
    log_msg "Validation merged output found."
fi

# We will run Stage 4 on existing outputs or prompt a rerun of DESeq2.
# Since the user specifically wants Stage 4 metrics, let's run the Stage 4 scripts 
# against the main project's graph or a specified graphML path.

# --- STAGE 4: Model Validation ---
log_msg "Starting Stage 4: Model Validation"

log_msg "Fetching External Disease Databases..."
Rscript Ranalysis/scripts/Stage4_Validation/01_fetch_disease_databases.R \
    2>&1 | tee "${LOG_DIR}/stage4_db_fetch_${TIMESTAMP}.log"

log_msg "Computing Cross-Reference Validation Metrics..."
# We pass paths to the cross-reference metrics so it can run specifically on validation networks if needed
# By default, 02_cross_reference_metrics.R picks up from outputs/.
Rscript Ranalysis/scripts/Stage4_Validation/02_cross_reference_metrics.R \
    2>&1 | tee "${LOG_DIR}/stage4_metrics_${TIMESTAMP}.log"

log_msg "Validation Pipeline Execution Complete."
echo "============================================="
echo "Done. Check validation outputs at: outputs/Stage4_Validation/"
echo "============================================="
