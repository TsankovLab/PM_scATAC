#!/bin/bash
# run_footprint_generic.sh  (R1_Q4)
# Generic ChromBPNet footprint script for one modisco pattern across all 6 NKT cell types.
# Args: $1=pattern_id (e.g. pos_patterns.pattern_5)
#       $2=core_seq   (e.g. TGAGTCAT)
#       $3=core_rc    (e.g. ATGACTCA)
#
# Runs chrombpnet footprints for 5 folds × 6 cell types, averages, and plots.

set -euo pipefail

ml cuda/11.7.0
ml cudnn/8.9.5-11
ml proxies
ml tensorrt/8.5.3.1

source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh
set +u; conda activate chrombpnet; set -u
export MPLBACKEND=Agg

PATTERN_ID="${1}"   # e.g. pos_patterns.pattern_5
CORE_SEQ="${2}"
CORE_RC="${3}"

# Sanitise for directory name: replace dots/slashes with underscores
PATTERN_SAFE=$(echo "${PATTERN_ID}" | tr './' '__')

MESO_DIR=/sc/arion/scratch/giottb01/Xmen/meso
GREFDIR=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet
OUT_DIR=${MESO_DIR}/footprints_${PATTERN_SAFE}
PLOT_SCRIPT=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/plot_footprint_generic.py

CELLTYPES=(CD4 CD8 CD8_exhausted FGFBP2_NK KLRC1_NK Tregs)

mkdir -p "${OUT_DIR}"

# Motif file
MOTIF_FILE=${OUT_DIR}/motif_file.txt
printf '%s\t%s\n%s\t%s\n' "${PATTERN_SAFE}" "${CORE_SEQ}" "${PATTERN_SAFE}" "${CORE_RC}" \
    > "${MOTIF_FILE}"

echo "=== Footprints for ${PATTERN_ID} (${CORE_SEQ}) ==="

for ct in "${CELLTYPES[@]}"; do
    echo "--- ${ct} ---"
    REGIONS=${MESO_DIR}/${ct}/${ct}_peakset_all_no_blacklist.bed
    for fold in 0 1 2 3 4; do
        OUT_PREFIX=${OUT_DIR}/${ct}_fold${fold}
        [ -f "${OUT_PREFIX}_footprints.h5" ] && echo "  fold ${fold}: skip" && continue
        chrombpnet footprints \
            -m  "${MESO_DIR}/${ct}/no_bias_model/fold_${fold}/models/chrombpnet_nobias.h5" \
            -r  "${REGIONS}" \
            -g  "${GREFDIR}/genome_references/hg38.genome.fa" \
            -fl "${GREFDIR}/folds/fold_${fold}.json" \
            -op "${OUT_PREFIX}" \
            -pwm_f "${MOTIF_FILE}"
        echo "  fold ${fold}: done"
    done
done

export HDF5_PLUGIN_PATH=${CONDA_PREFIX}/lib/hdf5/plugin

echo "=== Averaging and plotting ==="
python3 "${PLOT_SCRIPT}" \
    --out_dir     "${OUT_DIR}" \
    --pattern_id  "${PATTERN_SAFE}" \
    --motif_group "${PATTERN_SAFE}"

echo "=== Done: ${PATTERN_ID} ==="
