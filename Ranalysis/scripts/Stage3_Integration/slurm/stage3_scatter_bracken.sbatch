#!/bin/bash
#SBATCH --job-name=st3_scatter_b
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=Ranalysis/scripts/Stage3_Integration/slurm/logs/st3_scatter_b_%j.out

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

module purge
module load R

Rscript Ranalysis/scripts/Stage3_Integration/stage3_scatter.R \
  filtered_gene_counts=data/processed/filtered_gene_counts.csv \
  species_abundance=data/processed/brackenRanalysis/outputs/nonhuman_bracken260330op.csv \
  species_list_nonhuman=data/processed/brackenRanalysis/outputs/species_list_nonhuman260330op.csv \
  deseq_dir=outputs/DESeq2_results_bracken \
  padj_threshold=0.05 \
  log2fc_threshold=1.0 \
  num_batches=200 \
  out_dir=data/processed/stage3_batches_bracken_hpc
