#!/bin/bash
#SBATCH --job-name=st3_worker_single
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=Ranalysis/scripts/Stage3_Integration/slurm/logs/st3_worker_single_%j.out

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

module purge
module load gnu15/15.1.0
module load R/4.5.0

# Set environment for system dependencies
export PKG_CONFIG_PATH=/usr/lib64/pkgconfig:/usr/lib/pkgconfig:$PKG_CONFIG_PATH
export INCLUDE_DIR=/usr/include
export LIB_DIR=/usr/lib64
export CFLAGS="-I/usr/include -I/usr/include/curl"
export LDFLAGS="-L/usr/lib64"
export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH

DATA_DIR="${DATA_DIR:-data/processed/stage3_batches_bracken_hpc}"
GENE_CORR_THRESHOLD="${GENE_CORR_THRESHOLD:-0.6}"

# Accept batch ID as first argument instead of array task ID
BATCH_ID="${1:-}"
if [[ -z "${BATCH_ID}" ]]; then
  echo "BATCH_ID not provided. Usage: sbatch stage3_worker_single.sbatch <batch_id>"
  exit 1
fi

Rscript Ranalysis/scripts/Stage3_Integration/stage3_worker.R \
  batch_id=${BATCH_ID} \
  data_dir=${DATA_DIR} \
  gene_corr_threshold=${GENE_CORR_THRESHOLD}
