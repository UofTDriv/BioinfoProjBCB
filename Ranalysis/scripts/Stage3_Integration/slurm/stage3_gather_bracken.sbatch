#!/bin/bash
#SBATCH --job-name=st3_gather_b
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=Ranalysis/scripts/Stage3_Integration/slurm/logs/st3_gather_b_%j.out

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

Rscript Ranalysis/scripts/Stage3_Integration/stage3_gather.R \
  data_dir=data/processed/stage3_batches_bracken_hpc \
  annotated_species=outputs/brackenRanalysisAllReports/outputs/annotated_species260330.csv \
  species_corr_threshold=0.5 \
  padj_threshold=0.05 \
  log2fc_threshold=1.0 \
  out_graphml=outputs/Stage3_Network_for_RISK_bracken.graphml \
  out_homd_graphml=outputs/Stage3_Network_HOMD_Niche_bracken.graphml \
  out_annotation_csv=outputs/Stage3_RISK_Annotations_bracken.csv \
  out_validation_csv=outputs/Stage3_Jaccard_Validation_bracken.csv \
  out_fgsea_csv=outputs/Stage3_fgsea_go_bp_results_bracken.csv \
  out_species_layer=outputs/Stage3_Species_only_bracken.graphml \
  out_gene_layer=outputs/Stage3_Gene_only_bracken.graphml \
  export_phases=true
