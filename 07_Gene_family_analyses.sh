#!/bin/bash

#################################################
## CAFE_Run
###############################################

#SBATCH --job-name=cafe5
#SBATCH --output=logs/cafe5_%j.log
#SBATCH --error=logs/cafe5_%j.err
#SBATCH --ntasks=16
#SBATCH --mem=244G
#SBATCH --time=48:00:00


module purge
module load bear-apps/2023a
module load CAFE5/5.1.0-gfbf-2023a

set -euo pipefail

CAFE_DIR="/rds/projects/s/shorinoo-teff/genomes/pangenome/orthogroup/run5/cafe_input"
mkdir -p "${CAFE_DIR}/cafe_output"

sed -i 's/\r$//' cafe_gene_counts.tsv

cafe5 \
    -i          "${CAFE_DIR}/cafe_gene_counts.tsv" \
    -t           "${CAFE_DIR}/cafe_species_tree.nwk" \
    -k        4 \
    -c          "${SLURM_CPUS_PER_TASK:-16}" \
    -o  "${CAFE_DIR}/cafe_output/cafe_run"


######################################
## CAFE_Plot
#####################################
#SBATCH --job-name=plot_cafe
#SBATCH --output=logs/plot_cafe_%j.log
#SBATCH --error=logs/plot_cafe_%j.err
#SBATCH --ntasks=4
#SBATCH --time=05:00:00


module purge
module load bear-apps/2023a
module load Biopython/1.86-gfbf-2023a
module load Python/3.11.3-GCCcore-12.3.0
module load matplotlib/3.7.2-gfbf-2023a

cafeplotter -i . -o cafe_plot/png/ --expansion_color 'red' --contraction_color 'blue' --format 'png'
cafeplotter -i . -o cafe_plot/svg/ --expansion_color 'red' --contraction_color 'blue' --format 'svg'
cafeplotter -i . -o cafe_plot/pdf/ --expansion_color 'red' --contraction_color 'blue' --format 'pdf'

