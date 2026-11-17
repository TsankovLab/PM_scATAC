#!/bin/bash
# project_myeloid_on_pbmc.sh  (R1_Q4)
#
# Reference-maps all myeloid scATAC cells onto the Greenleaf PBMC/BM
# IterativeLSI space to obtain BioClassification labels.
#
# - Query  : myeloid_cells/scatac_ArchR (celltype_main="Myeloid", all cells)
# - Ref    : greenleaf_PBMC_BM (BioClassification column)
# - Method : uwot transform on IterativeLSI-based UMAP (run_projection_default.sh)
#
# All myeloid cells are projected in a SINGLE call using celltype_main=Myeloid,
# so Arrow files are copied only once (vs 7x when looping per celltype).
# Prerequisite: run submit_add_myeloid_compartment.sh first to add that column.
#
# Outputs → R1_Q4/output/pbmc_projection/
#   myeloid_all_knn_labels.csv   KNN labels + confidence for all cells
#   myeloid_all_UMAP_*.pdf       query cells overlaid on reference UMAP
#
# Requires: ~64 GB RAM, 4 cores, ~2 h
# Submit  : bsub < project_myeloid_on_pbmc.sh
#BSUB -J myeloid_pbmc_proj_R1Q4
#BSUB -P <LSF_ACCOUNT>
#BSUB -q premium
#BSUB -n 4
#BSUB -R "rusage[mem=16000]"
#BSUB -W 12:00
#BSUB -o logs/myeloid_pbmc_proj_R1Q4_%J.out
#BSUB -e logs/myeloid_pbmc_proj_R1Q4_%J.err

source "${MESO_CONFIG:-$(for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$PWD"; do while [ ! -f "$d/config.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; [ -f "$d/config.sh" ] && { echo "$d/config.sh"; break; }; done)}"   # data / results / reference roots: see config.sh at the repo root
set -euo pipefail

# ── Rscript from meso_scatac conda env ───────────────────────────────────────
export PATH=${CONDA_ENVS}/meso_scatac/bin:${PATH}
export OPENBLAS_NUM_THREADS=1

# ── Scratch tempdir (uwot::load_uwot untars a large model) ───────────────────
export TMPDIR=${SCRATCH_DIR}/tmp
mkdir -p "${TMPDIR}"

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECTION_WRAPPER=${SCRATCH_DIR}/MesoRNA/ProjectionCode/run_projection_default.sh

QUERY_DIR=${DATA_DIR}/myeloid_cells/scatac_ArchR
REF_DIR=${REF_DIR}/Xmen/greenleaf_PBMC_BM

OUTPUT_DIR=${RESULTS_DIR}/R1_Q4/output/pbmc_projection
mkdir -p "${OUTPUT_DIR}"

# ── Project all myeloid cells in one call ────────────────────────────────────
echo "========================================================"
echo " Projecting all myeloid cells onto Greenleaf PBMC/BM reference"
echo "========================================================"

QUERY_DIR="${QUERY_DIR}" \
REF_DIR="${REF_DIR}" \
SUBSET_COL=celltype_main \
SUBSET_VALUE=Myeloid \
REF_LABEL_COL=BioClassification \
OUTPUT_DIR="${OUTPUT_DIR}" \
OUTPUT_PREFIX=myeloid_all \
CORES=4 \
GENOME=hg19 \
TMPDIR="${TMPDIR}" \
    bash "${PROJECTION_WRAPPER}" --save_rds

echo "========================================================"
echo " Done. Outputs in: ${OUTPUT_DIR}"
echo "========================================================"
