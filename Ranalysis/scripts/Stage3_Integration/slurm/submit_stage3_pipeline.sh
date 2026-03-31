#!/bin/bash
set -euo pipefail

MODE="${1:-bracken}"
BASE_DIR="/home/wranalab/derivard/BCB330Proj/BioinfoProjBCB/Ranalysis/scripts/Stage3_Integration/slurm"

if [[ "${MODE}" == "bracken" ]]; then
  SCATTER_SCRIPT="${BASE_DIR}/stage3_scatter_bracken.sbatch"
  GATHER_SCRIPT="${BASE_DIR}/stage3_gather_bracken.sbatch"
  DATA_DIR="data/processed/stage3_batches_bracken_hpc"
elif [[ "${MODE}" == "unaligned" ]]; then
  SCATTER_SCRIPT="${BASE_DIR}/stage3_scatter_unaligned.sbatch"
  GATHER_SCRIPT="${BASE_DIR}/stage3_gather_unaligned.sbatch"
  DATA_DIR="data/processed/stage3_batches_unaligned_hpc"
else
  echo "Usage: bash ${BASE_DIR}/submit_stage3_pipeline.sh [bracken|unaligned]"
  exit 1
fi

mkdir -p "${BASE_DIR}/logs"

echo "Submitting scatter job (${MODE})..."
SCATTER_JID=$(sbatch --parsable "${SCATTER_SCRIPT}")
echo "Scatter job id: ${SCATTER_JID}"

echo "Submitting worker array after scatter completion..."
# Keep array size synced with scatter num_batches value in sbatch scripts.
WORKER_JID=$(sbatch --parsable \
  --dependency=afterok:${SCATTER_JID} \
  --array=1-50 \
  --export=ALL,DATA_DIR=${DATA_DIR},GENE_CORR_THRESHOLD=0.6 \
  "${BASE_DIR}/stage3_worker_array.sbatch")

echo "Worker array job id: ${WORKER_JID}"

echo "Submitting gather job after worker array completion..."
GATHER_JID=$(sbatch --parsable \
  --dependency=afterok:${WORKER_JID} \
  "${GATHER_SCRIPT}")

echo "Gather job id: ${GATHER_JID}"
echo "Pipeline submitted: scatter -> worker array -> gather"
