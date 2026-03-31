#!/bin/bash
#SBATCH --job-name=st3_worker
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=Ranalysis/scripts/Stage3_Integration/slurm/logs/st3_worker_%A_%a.out

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

module purge
module load R

DATA_DIR="${DATA_DIR:-data/processed/stage3_batches_bracken_hpc}"
GENE_CORR_THRESHOLD="${GENE_CORR_THRESHOLD:-0.6}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "SLURM_ARRAY_TASK_ID is not set. Submit this script with --array=1-N."
  exit 1
fi

Rscript Ranalysis/scripts/Stage3_Integration/stage3_worker.R \
  batch_id=${SLURM_ARRAY_TASK_ID} \
  data_dir=${DATA_DIR} \
  gene_corr_threshold=${GENE_CORR_THRESHOLD}
