#!/bin/bash
# modisco_postmerge.sh
# Runs after all modisco merge batch jobs complete.
# Args: $1=chromBPdir  $2=model_head  $3=repodir  $4=celltypes (space-sep, quoted)
#
# Steps:
#   4a. Compile merged h5s into one
#   4b. TOMTOM matching against known motifs
#   5.  Re-run FI-NeMo per cell type with compiled merged motifs

set -euo pipefail

ml anaconda3/latest

CHROMBPDIR=${1}
MODEL_HEAD=${2}
REPODIR=${3}
shift 3
CELLTYPES=("$@")

MERGED_DIR=${CHROMBPDIR}/modisco_merged_${MODEL_HEAD}
CLUSTER_KEY=${MERGED_DIR}/cluster_key_combined.txt
MOTIF_DB=/sc/arion/work/giottb01/conda/envs/chrombpnet/lib/python3.8/site-packages/chrombpnet/data/motifs.meme.txt
HDMA_UTILS=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/HDMA/code/03-chrombpnet
PYTHON=/sc/arion/work/giottb01/conda/envs/chrombpnet/bin/python

echo "=== Post-merge pipeline: model_head=${MODEL_HEAD} | cell types: ${CELLTYPES[*]} ==="

# ---- 4a. Compile all merged h5s into one ------------------------------------
echo "=== Step 4a: compile merged modisco ==="
source activate chrombpnet
mkdir -p ${MERGED_DIR}/compiled

PYTHONPATH=${HDMA_UTILS} ${PYTHON} ${REPODIR}/utils/04a-compile_modisco_obj.py \
    --output_dir ${MERGED_DIR} \
    --cluster_key_path ${CLUSTER_KEY} \
    --compiled_modisco_h5_path ${MERGED_DIR}/compiled/modisco_compiled.h5 \
    --compiled_modisco_tsv_path ${MERGED_DIR}/compiled/modisco_compiled.tsv \
    --modisco_merged_dir ${MERGED_DIR}/merged_motifs

echo "  Compiled: $(ls ${MERGED_DIR}/compiled/)"

# ---- 4b. TOMTOM matching ----------------------------------------------------
echo "=== Step 4b: TOMTOM matching ==="
mkdir -p ${MERGED_DIR}/compiled_tomtom

PYTHONPATH=${HDMA_UTILS} ${PYTHON} ${REPODIR}/utils/04b-get_tomtom_matches.py \
    --modisco-h5 ${MERGED_DIR}/compiled/modisco_compiled.h5 \
    --out-dir ${MERGED_DIR}/compiled_tomtom \
    --meme-db ${MOTIF_DB} \
    --verbose True

echo "  TOMTOM done: $(ls ${MERGED_DIR}/compiled_tomtom/)"

# ---- 5. FI-NeMo per cell type with compiled merged motifs -------------------
echo "=== Step 5: FI-NeMo on merged motifs ==="
source activate finemo

for ct in "${CELLTYPES[@]}"; do
    echo "  finemo: ${ct}"
    bash ${REPODIR}/utils/finemo_motif_calls_on_merged.sh \
        "${CHROMBPDIR}" "${ct}" "${MODEL_HEAD}"
done

echo "=== Post-merge pipeline complete ==="
