#!/usr/bin/env bash
##############################################################################
# functional_annotation_pipeline.sh
#
# Teff pangenome: proteome extraction -> OrthoFinder -> pull out the
# core / dispensable / unique gene sets -> pick one representative gene
# per orthogroup -> functional annotation (InterProScan + eggNOG-mapper).
#
# Requirements:
#   - conda environments: seqkit, agat, orthofinder, eggnog (eggnog-mapper)
#   - a local InterProScan install
#   - gene lists produced upstream (see check_core_disp_uniqu.R):
#       <LISTS_DIR>/core_genes.txt, disp_genes.txt, uniq_genes.txt
#       <LISTS_DIR>/single_core_in_og.txt, single_disp_in_og.txt, single_uniq_in_og.txt
#
# Usage:
#   ./functional_annotation_pipeline.sh              # run every step
#   ./functional_annotation_pipeline.sh proteomes     # run one step
#   ./functional_annotation_pipeline.sh og_fastas annotate_representatives
#
# Steps: proteomes orthofinder gene_set_fastas og_fastas
#        representative_fastas annotate_representatives
#
# All paths below can be overridden from the environment, e.g.:
#   THREADS=64 ROOT=/data/e.teff ./functional_annotation_pipeline.sh
##############################################################################
set -euo pipefail

############################ CONFIGURATION ###################################
: "${ROOT:=/home/shared/shared_all/e.teff}"

: "${GENOME_DIR:=${ROOT}/genomes/teff/assembly/chr_genome/fasta}"   # input *.fasta
: "${ANNO_DIR:=${ROOT}/genes_anno/helixer_annotation_EDTAfilt}"     # <geno>.filt.EDTA.gff3
: "${PROT_DIR:=${ROOT}/genomes/teff/protein_sequences/protein_fasta_files}"
: "${FUNC_DIR:=${ROOT}/functional_anno}"
: "${LISTS_DIR:=${FUNC_DIR}/v1}"          # gene-ID lists from upstream classification
: "${OG_DIR:=${LISTS_DIR}/OG}"            # per-orthogroup fastas
: "${INTERPRO_OUT:=${LISTS_DIR}/interpro_out}"
: "${EGGNOG_OUT:=${LISTS_DIR}/eggnog_out}"

: "${THREADS:=28}"
: "${OF_THREADS:=96}"
: "${OF_ALIGN_THREADS:=8}"
: "${OF_INFLATION:=1.5}"

: "${CONDA_SH:=${HOME}/miniforge3/etc/profile.d/conda.sh}"
: "${ENV_SEQKIT:=seqkit}"
: "${ENV_AGAT:=agat}"
: "${ENV_ORTHOFINDER:=orthofinder}"
: "${ENV_EGGNOG:=eggnog}"

: "${INTERPROSCAN_SH:=${FUNC_DIR}/v0/interpro/my_interproscan/interproscan-5.77-108.0/interproscan.sh}"
: "${EGGNOG_DIR:=${LISTS_DIR}/eggnog-mapper-2.1.14}"

# gene-set name -> gene-ID list file, used by both gene_set_fastas and
# annotate_representatives (representatives use the "single_..._in_og" lists)
declare -A GENE_SETS=(
  [core_genes]="${LISTS_DIR}/core_genes.txt"
  [disp_genes]="${LISTS_DIR}/disp_genes.txt"
  [uniq_genes]="${LISTS_DIR}/uniq_genes.txt"
)
declare -A REPRESENTATIVE_SETS=(
  [single_core_in_og]="${LISTS_DIR}/single_core_in_og.txt"
  [single_disp_in_og]="${LISTS_DIR}/single_disp_in_og.txt"
  [single_uniq_in_og]="${LISTS_DIR}/single_uniq_in_og.txt"
)

############################### HELPERS ######################################
log() { printf '\n[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

activate() {                       # conda activate that survives 'set -u'
  set +u
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "$1"
  set -u
}

mkdir -p "${PROT_DIR}" "${OG_DIR}/fasta" "${INTERPRO_OUT}" "${EGGNOG_OUT}"

################################ STEPS #######################################

# 1. FASTA -> 60 chars/line (AGAT requirement) -> protein FASTA per genome
step_proteomes() {
  log "Extracting proteomes"
  activate "${ENV_SEQKIT}"
  for file in "${GENOME_DIR}"/*.fasta; do
    geno=$(basename "${file}" | cut -d'_' -f1)
    gff=${ANNO_DIR}/${geno}.filt.EDTA.gff3
    [[ -f ${gff} ]] || { echo "  !! no annotation for ${geno}, skipping"; continue; }

    wrapped=${GENOME_DIR}/${geno}.60w.fasta
    seqkit seq -w 60 "${file}" > "${wrapped}"

    activate "${ENV_AGAT}"
    agat_sp_extract_sequences.pl \
      -g "${gff}" -f "${wrapped}" -p \
      -o "${PROT_DIR}/${geno}.protein.fa"
    activate "${ENV_SEQKIT}"

    echo "  ${geno}: $(grep -c '^>' "${PROT_DIR}/${geno}.protein.fa") proteins"
  done
}

# 2. OrthoFinder: one full run (gene trees) over the whole proteome set
step_orthofinder() {
  log "Running OrthoFinder"
  activate "${ENV_ORTHOFINDER}"
  orthofinder -f "${PROT_DIR}" \
    -M msa -S diamond -T iqtree -A mafft \
    -I "${OF_INFLATION}" -t "${OF_THREADS}" -a "${OF_ALIGN_THREADS}"
}

# 3. Pull the full core / dispensable / unique gene sets out of all proteomes
step_gene_set_fastas() {
  log "Extracting core / dispensable / unique gene sets"
  activate "${ENV_SEQKIT}"

  local all_proteins=${PROT_DIR}/all_proteins.fa
  cat "${PROT_DIR}"/*.protein.fa > "${all_proteins}"

  for set in "${!GENE_SETS[@]}"; do
    list=${GENE_SETS[$set]}
    [[ -f ${list} ]] || { echo "  !! missing ${list}, skipping ${set}"; continue; }
    seqkit grep -f "${list}" "${all_proteins}" > "${PROT_DIR}/${set}.fa"
    sed -i 's/\*//g' "${PROT_DIR}/${set}.fa"   # strip stop-codon characters
    echo "  ${set}: $(grep -c '^>' "${PROT_DIR}/${set}.fa") sequences"
  done
}

# 4. One FASTA per orthogroup (used to pick representatives downstream)
step_og_fastas() {
  log "Building per-orthogroup FASTA files"
  activate "${ENV_SEQKIT}"
  local all_proteins=${PROT_DIR}/all_proteins.fa

  for file in "${LISTS_DIR}"/*in_og.txt; do
    og=$(basename "${file}" .txt)
    seqkit grep -f "${file}" "${all_proteins}" > "${OG_DIR}/fasta/${og}.fa"
  done
}

# 5. Representative sequences (one gene per orthogroup, per category)
step_representative_fastas() {
  log "Extracting representative sequences per orthogroup"
  activate "${ENV_SEQKIT}"
  local all_proteins=${PROT_DIR}/all_proteins.fa

  for set in "${!REPRESENTATIVE_SETS[@]}"; do
    list=${REPRESENTATIVE_SETS[$set]}
    [[ -f ${list} ]] || { echo "  !! missing ${list}, skipping ${set}"; continue; }
    seqkit grep -f "${list}" "${all_proteins}" > "${LISTS_DIR}/${set}.fa"
    sed -i 's/\*//g' "${LISTS_DIR}/${set}.fa"
    echo "  ${set}: $(grep -c '^>' "${LISTS_DIR}/${set}.fa") sequences"
  done
}

# 6. Functional annotation of representatives only (full gene sets are too slow)
step_annotate_representatives() {
  log "Annotating representative sequences (InterProScan + eggNOG-mapper)"
  for set in "${!REPRESENTATIVE_SETS[@]}"; do
    fa=${LISTS_DIR}/${set}.fa
    [[ -f ${fa} ]] || { echo "  !! missing ${fa}, run representative_fastas first"; continue; }

    sh "${INTERPROSCAN_SH}" -i "${fa}" -f tsv -goterms --cpu "${THREADS}" \
       -o "${INTERPRO_OUT}/interpro.${set}.tsv"

    activate "${ENV_EGGNOG}"
    "${EGGNOG_DIR}/emapper.py" --data_dir "${EGGNOG_DIR}/data" \
       -i "${fa}" --output "${set}.eggnog" --output_dir "${EGGNOG_OUT}" \
       --cpu "${THREADS}" --override
  done
}

############################### MAIN #########################################
# One-time eggNOG-mapper environment setup (not part of this pipeline):
#   conda create -n eggnog -c conda-forge -c bioconda eggnog-mapper
#   download_eggnog_data.py --data_dir "${EGGNOG_DIR}/data"

ALL_STEPS=(proteomes orthofinder gene_set_fastas og_fastas
           representative_fastas annotate_representatives)
STEPS=("${@:-all}")
[[ ${STEPS[0]} == all ]] && STEPS=("${ALL_STEPS[@]}")

for s in "${STEPS[@]}"; do
  declare -F "step_${s}" >/dev/null || { echo "Unknown step: ${s}" >&2; exit 1; }
  "step_${s}"
done
log "Done."