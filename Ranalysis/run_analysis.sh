#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name=r_analysis
#SBATCH --mail-user=denis.rivard@mail.utoronto.ca
#SBATCH --mail-type=FAIL,END
#SBATCH -o /home/wranalab/derivard/2026metaseq_pipeline/new_pipeline/jobs/%x_%j.out

# Run the R analysis script
Rscript /home/wranalab/derivard/2026metaseq_pipeline/new_pipeline/Ranalysis/scripts/Stage3_Integration/main_gene_path_associations.R
