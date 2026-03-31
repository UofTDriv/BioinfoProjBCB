#!/usr/bin/env bash

# Bioinformatics Pipeline Orchestration Script
# Author: Denis Rivard
# Description: Automates the execution of Stage 1, Stage 2, and Stage 3 R scripts.
# Automatically detects if Bracken output exists and runs the appropriate Stage 2 & Stage 3 scripts.

set -e # Exit immediately if a command exits with a non-zero status
set -u # Treat unset variables as an error
set -o pipefail

# Define variables
PROJECT_DIR=$(pwd)
INPUT_DIR="data/input"
PROJ_NAME="BioinfoProj"
LOG_DIR="outputs/logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

echo "============================================="
echo "Starting Bioinformatics Pipeline Overview"
echo "Project Directory: ${PROJECT_DIR}"
echo "Start Time: ${TIMESTAMP}"
echo "============================================="

# Ensure log directory exists
mkdir -p "${LOG_DIR}"

log_msg() {
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] $1" | tee -a "${LOG_DIR}/pipeline_${TIMESTAMP}.log"
}

# --- STAGE 1: Data Acquisition ---
log_msg "Starting Stage 1: Data Acquisition (reportRanalysis.R)"
Rscript Ranalysis/scripts/Stage1_Data_acquisition/reportRanalysis.R \
    --input-dir "${INPUT_DIR}" \
    --proj-name "${PROJ_NAME}" \
    2>&1 | tee "${LOG_DIR}/stage1_${TIMESTAMP}.log"

# Check if bracken successfully ran
BRACKEN_RAN=false
if ls data/processed/brackenRanalysis/outputs/nonhuman_bracken*.csv 1> /dev/null 2>&1; then
    BRACKEN_RAN=true
    log_msg "Detected Bracken output. Will execute Bracken-specific tracks."
else
    log_msg "No Bracken output detected. Will proceed with unaligned track only."
fi

# --- STAGE 2: Data Processing (DESeq2) ---
log_msg "Starting Stage 2: Data Processing (DESeq2)"

# Unaligned track
log_msg "Running unaligned DESeq2..."
Rscript Ranalysis/scripts/Stage2_Data_Processing/run_deseq2_species_gene.R \
    2>&1 | tee "${LOG_DIR}/stage2_unaligned_${TIMESTAMP}.log"

# Bracken track
if [ "$BRACKEN_RAN" = true ]; then
    log_msg "Running Bracken DESeq2..."
    Rscript Ranalysis/scripts/Stage2_Data_Processing/run_deseq2_species_gene_bracken.R \
        2>&1 | tee "${LOG_DIR}/stage2_bracken_${TIMESTAMP}.log"
fi

# --- STAGE 3: Integration ---
log_msg "Starting Stage 3: Integration (Network Generation & Fusion)"

log_msg "Running Combined Network Fusion (Unaligned)..."
Rscript Ranalysis/scripts/Stage3_Integration/Stage3_Combined_Network_Fusion.R \
    --deseq-dir=outputs/DESeq2_results \
    --out-graphml=outputs/Stage3_Network_for_RISK.graphml \
    2>&1 | tee "${LOG_DIR}/stage3_network_fusion_${TIMESTAMP}.log"

if [ "$BRACKEN_RAN" = true ]; then
    log_msg "Running Combined Network Fusion (Bracken)..."
    Rscript Ranalysis/scripts/Stage3_Integration/Stage3_Combined_Network_Fusion_bracken.R \
        --deseq-dir=outputs/DESeq2_results_bracken \
        --out-graphml=outputs/Stage3_Network_for_RISK_bracken.graphml \
        2>&1 | tee "${LOG_DIR}/stage3_network_fusion_bracken_${TIMESTAMP}.log"
fi

# --- STAGE 4: Model Validation ---
log_msg "Starting Stage 4: Model Validation"

log_msg "Fetching External Disease Databases..."
Rscript Ranalysis/scripts/Stage4_Validation/01_fetch_disease_databases.R \
    2>&1 | tee "${LOG_DIR}/stage4_db_fetch_${TIMESTAMP}.log"

log_msg "Computing Cross-Reference Validation Metrics..."
Rscript Ranalysis/scripts/Stage4_Validation/02_cross_reference_metrics.R \
    2>&1 | tee "${LOG_DIR}/stage4_metrics_${TIMESTAMP}.log"

log_msg "Validation Pipeline Execution Complete."

log_msg "Pipeline execution complete!"
echo "============================================="
echo "Done. Check logs at: ${LOG_DIR}"
echo "============================================="
