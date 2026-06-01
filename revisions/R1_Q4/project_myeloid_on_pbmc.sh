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
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 4
#BSUB -R "rusage[mem=16000]"
#BSUB -W 12:00
#BSUB -o /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/myeloid_pbmc_proj_R1Q4_%J.out
#BSUB -e /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/myeloid_pbmc_proj_R1Q4_%J.err

set -euo pipefail

# ── Rscript from meso_scatac conda env ───────────────────────────────────────
export PATH=/sc/arion/work/giottb01/conda/envs/meso_scatac/bin:${PATH}
export OPENBLAS_NUM_THREADS=1

# ── Scratch tempdir (uwot::load_uwot untars a large model) ───────────────────
export TMPDIR=/sc/arion/scratch/giottb01/tmp
mkdir -p "${TMPDIR}"

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECTION_WRAPPER=/sc/arion/scratch/leew17/MesoRNA/ProjectionCode/run_projection_default.sh

QUERY_DIR=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR
REF_DIR=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Xmen/greenleaf_PBMC_BM

OUTPUT_DIR=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/pbmc_projection
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
