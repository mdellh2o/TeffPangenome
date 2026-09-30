
#!/bin/bash
#SBATCH --job-name=orthofinder_cafe
#SBATCH --output=logs/orthofinder_cafe_%j.log
#SBATCH --error=logs/orthofinder_cafe_%j.err
#SBATCH --ntasks=16
#SBATCH --mem=244G
#SBATCH --time=48:00:00

####################################
### Ultrametric_Species_Tree
#######################################

# =============================================================================
# OrthoFinder → Ultrametric Species Tree (ape::chronos) → CAFE5 Pipeline
#
# Usage:
#   sbatch run_orthofinder_cafe.sh sp1.fa sp2.fa sp3.fa sp4.fa sp5.fa sp6.fa
#
# Ultrametric dating:
#   ape::chronos() uses penalised likelihood to fit a molecular clock and
#   produces a fully ultrametric tree (all tips equidistant from the root)
#   with no calibration file required. Branch lengths are in relative time
#   units. Three clock models are available: "correlated" (default, relaxed),
#   "relaxed" (uncorrelated), and "discrete". Lambda controls smoothing
#   (higher = closer to strict clock; 0 = fully relaxed).
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# 0. Load modules (adjust names/versions to match your HPC environment)
# ---------------------------------------------------------------------------

module purge
module load bear-apps/2023a
module load OrthoFinder/2.5.5-foss-2023a
module load MAFFT/7.520-GCC-12.3.0-with-extensions
module load FastTree/2.1.11-GCCcore-12.3.0
module load R/4.4.1-gfbf-2023a # ape::chronos for ultrametric conversion
module load Python/3.11.3-GCCcore-12.3.0
module load DIAMOND/2.1.8-GCC-12.3.0 # OrthoFinder dependency


# ---------------------------------------------------------------------------
# 1. Validate input
# ---------------------------------------------------------------------------
if [[ $# -ne 7 ]]; then # Seyi Note: Changed to 7
    echo "ERROR: Exactly 8 protein FASTA files required."
    echo "Usage: sbatch $0 sp1.fa sp2.fa sp3.fa sp4.fa sp5.fa sp6.fa sp7.fa sp8.fa"
    exit 1
fi

PROTEINS=("$@")
THREADS=16

# ---------------------------------------------------------------------------
# 2. Set up directory structure
# ---------------------------------------------------------------------------
WORKDIR="${SLURM_SUBMIT_DIR}"
FASTA_DIR="${WORKDIR}/proteomes"
OF_OUTDIR="${WORKDIR}/orthofinder_results"
TREE_DIR="${WORKDIR}/species_tree"
CAFE_DIR="${WORKDIR}/cafe_input"
LOG_DIR="${WORKDIR}/logs"

mkdir -p "${FASTA_DIR}" "${TREE_DIR}" "${CAFE_DIR}" "${LOG_DIR}"

# ---------------------------------------------------------------------------
# 3. Copy and validate protein FASTAs
# ---------------------------------------------------------------------------
echo "[$(date)] Copying input FASTA files..."
for FA in "${PROTEINS[@]}"; do
    if [[ ! -f "${FA}" ]]; then
        echo "ERROR: File not found: ${FA}"
        exit 1
    fi
    SPECIES=$(basename "${FA}" | sed 's/\.[^.]*$//')
    DEST="${FASTA_DIR}/${SPECIES}.fa"
    cp "${FA}" "${DEST}"
    echo "  Copied: ${FA} → ${DEST}"
done

# ---------------------------------------------------------------------------
# 4. Run OrthoFinder
# ---------------------------------------------------------------------------
echo "[$(date)] Running OrthoFinder..."
orthofinder \
    -f "${FASTA_DIR}" \
    -t "${THREADS}" \
    -a "${THREADS}" \
    -o "${OF_OUTDIR}" \
    -S diamond \
    -M msa \
    -T fasttree

# Locate the most recent OrthoFinder results folder
OF_RESULTS=$(ls -td "${OF_OUTDIR}"/Results_* 2>/dev/null | head -1)
if [[ -z "${OF_RESULTS}" ]]; then
    echo "ERROR: OrthoFinder output directory not found."
    exit 1
fi
echo "[$(date)] OrthoFinder results: ${OF_RESULTS}"

OF_SPECIES_TREE="${OF_RESULTS}/Species_Tree/SpeciesTree_rooted.txt"
if [[ ! -f "${OF_SPECIES_TREE}" ]]; then
    echo "ERROR: OrthoFinder species tree not found: ${OF_SPECIES_TREE}"
    exit 1
fi
cp "${OF_SPECIES_TREE}" "${TREE_DIR}/species_tree_raw.nwk"

# ---------------------------------------------------------------------------
# 5. Build an ultrametric species tree with ape::chronos()
#
#    chronos() fits a penalised-likelihood clock model with no calibration
#    needed. Key parameters:
#      lambda : smoothing penalty (1 = default; increase toward strict clock,
#               decrease toward fully relaxed; try 0.1–10 if tree looks odd)
#      model  : "correlated"  — rates autocorrelated along branches (default)
#               "relaxed"     — rates uncorrelated (each branch independent)
#               "discrete"    — finite number of discrete rate classes
#    The resulting branch lengths are in relative time units (dimensionless),
#    which is sufficient for CAFE5 to estimate gene family evolution rates.
# ---------------------------------------------------------------------------
echo "[$(date)] Converting species tree to ultrametric with ape::chronos()..."

cat > "${TREE_DIR}/make_ultrametric.R" << 'RSCRIPT'

#install.packages("ape", repos = "https://cran.r-project.org")
library(ape)

args     <- commandArgs(trailingOnly = TRUE)
in_tree  <- args[1]
out_tree <- args[2]

tree <- read.tree(in_tree)

# Ensure the tree is rooted (OrthoFinder output should be, but check anyway)
if (!is.rooted(tree)) {
    stop("Input tree is unrooted. Please provide a rooted Newick file.")
}

# Remove any zero-length branches that can cause chronos() to fail
tree$edge.length[tree$edge.length <= 0] <- 1e-6

# Fit penalised-likelihood clock — no calibration object needed
# (chronos() defaults to placing all tips at equal depth)
ultra <- chronos(
    tree,
    lambda = 1,        # smoothing parameter; adjust if convergence fails
    model  = "correlated",
    control = chronos.control(
        iter.max   = 1e4,
        eval.max   = 1e4,
        tol        = 1e-6
    )
)
class(ultra) <- "phylo"

# Fix any tiny negative branch lengths from numerical rounding
ultra$edge.length[ultra$edge.length < 0] <- 0


# Verify ultrametricity before writing (using ape's built-in check)
if (!is.ultrametric(ultra, tol = 1e-4)) {
    warning("Tree may not be fully ultrametric. Consider adjusting lambda.")
} else {
    cat("OK: Tree is ultrametric\n")
}

# Report tip depth range for diagnostic purposes
tip_depths <- node.depth.edgelength(ultra)[1:length(ultra$tip.label)]
range_d <- max(tip_depths) - min(tip_depths)
cat(sprintf("Tip depth range: %.2e\n", range_d))



write.tree(ultra, file = out_tree)
cat("Ultrametric tree written to:", out_tree, "\n")
RSCRIPT

Rscript "${TREE_DIR}/make_ultrametric.R" \
    "${TREE_DIR}/species_tree_raw.nwk" \
    "${TREE_DIR}/species_tree_ultrametric.nwk" \
    2>&1 | tee "${LOG_DIR}/chronos.log"

if [[ ! -f "${TREE_DIR}/species_tree_ultrametric.nwk" ]]; then
    echo "ERROR: Ultrametric tree not produced. Check ${LOG_DIR}/chronos.log"
    exit 1
fi
echo "[$(date)] Ultrametric tree: ${TREE_DIR}/species_tree_ultrametric.nwk"

# ---------------------------------------------------------------------------
# 6. Prepare CAFE5 input files
#    (a) Gene-family count table — tab-separated:
#        Desc <tab> Family ID <tab> sp1 <tab> … <tab> spN
#    (b) Ultrametric species tree in Newick format
# ---------------------------------------------------------------------------
echo "[$(date)] Preparing CAFE5 input..."

OF_ORTHOGROUPS="${OF_RESULTS}/Orthogroups/Orthogroups.GeneCount.tsv"
if [[ ! -f "${OF_ORTHOGROUPS}" ]]; then
    echo "ERROR: OrthoFinder gene-count table not found: ${OF_ORTHOGROUPS}"
    exit 1
fi

cat > "${CAFE_DIR}/prepare_cafe.py" << 'PYSCRIPT'
#!/usr/bin/env python3
"""
Convert OrthoFinder Orthogroups.GeneCount.tsv → CAFE5 input format.

CAFE5 format (tab-separated):
  Desc        Family ID   Species1  Species2  ...
  (null)      OG0000001   3         1         ...

Rules applied:
  - All-zero families are dropped.
  - Families present in at least one species are kept.
  - CAFE5 handles zero counts in individual species natively.
"""
import sys
import csv

in_file, out_file = sys.argv[1], sys.argv[2]

with open(in_file, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    species_cols = [f for f in reader.fieldnames
                    if f not in ("Orthogroup", "Total")]
    rows_out = [row for row in reader
                if sum(int(row[sp]) for sp in species_cols) > 0]

with open(out_file, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["Desc", "Family ID"] + species_cols)
    for row in rows_out:
        writer.writerow(["(null)", row["Orthogroup"]]
                        + [int(row[sp]) for sp in species_cols])

print(f"Written {len(rows_out)} orthogroups to {out_file}")
PYSCRIPT

python3 "${CAFE_DIR}/prepare_cafe.py" \
    "${OF_ORTHOGROUPS}" \
    "${CAFE_DIR}/cafe_gene_counts.tsv"

cp "${TREE_DIR}/species_tree_ultrametric.nwk" "${CAFE_DIR}/cafe_species_tree.nwk"

# ---------------------------------------------------------------------------
# 7. Generate a ready-to-submit CAFE5 job script
# ---------------------------------------------------------------------------
cat > "${CAFE_DIR}/run_cafe.sh" << 'CAFESCRIPT'
#!/bin/bash
#SBATCH --job-name=cafe5
#SBATCH --output=logs/cafe5_%j.log
#SBATCH --error=logs/cafe5_%j.err
#SBATCH --ntasks=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00


module purge
module load bear-apps/2023a
module load CAFE5/5.1.0-gfbf-2023a

set -euo pipefail

CAFE_DIR="$(dirname "$(realpath "$0")")"
mkdir -p "${CAFE_DIR}/cafe_output"

sed -i 's/\r$//' cafe_gene_counts.tsv

cafe5 \
    -i          "${CAFE_DIR}/cafe_gene_counts.tsv" \
    -t           "${CAFE_DIR}/cafe_species_tree.nwk" \
    -k        4 \
    -c          "${SLURM_CPUS_PER_TASK:-8}" \
    -o  "${CAFE_DIR}/cafe_output/cafe_run"
CAFESCRIPT

chmod +x "${CAFE_DIR}/run_cafe.sh"

# ---------------------------------------------------------------------------
# 8. Summary
# ---------------------------------------------------------------------------
echo ""
echo "============================================================"
echo " Pipeline complete — $(date)"
echo "============================================================"
echo ""
echo " OrthoFinder results  : ${OF_RESULTS}"
echo " Raw species tree     : ${TREE_DIR}/species_tree_raw.nwk"
echo " chronos log          : ${LOG_DIR}/chronos.log"
echo " Ultrametric tree     : ${TREE_DIR}/species_tree_ultrametric.nwk"
echo ""
echo " CAFE5 input files:"
echo "   Gene counts table  : ${CAFE_DIR}/cafe_gene_counts.tsv"
echo "   Species tree       : ${CAFE_DIR}/cafe_species_tree.nwk"
echo "   CAFE5 run script   : ${CAFE_DIR}/run_cafe.sh"
echo ""
echo " To run CAFE5:"
echo "   sbatch ${CAFE_DIR}/run_cafe.sh"
echo "============================================================"
