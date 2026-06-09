#!/bin/bash
# run_pattern24_footprints.sh  (R1_Q4)
#
# Computes ChromBPNet footprints for pos_patterns.pattern_24 (AP-1 + NR4A2
# composite) across all 6 NKT cell types, using KLRC1_NK peaks as the
# common region set so results are directly comparable.
#
# Follows the chrombpnet_footprint.sh pattern in git_repo/utils/:
#   1. Run chrombpnet footprints for each fold (0-4) for each cell type
#   2. Average across folds per cell type
#   3. Overlay all cell types in one plot
#
# Submit: bsub < run_pattern24_footprints.sh
#BSUB -J pattern24_footprints
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q gpu
#BSUB -n 4
#BSUB -W 12:00
#BSUB -gpu "num=1"
#BSUB -R "a100"
#BSUB -R "rusage[mem=32000]"
#BSUB -R "span[hosts=1]"
#BSUB -o /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/pattern24_footprints_%J.out
#BSUB -e /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/pattern24_footprints_%J.err

set -euo pipefail

# Load system CUDA modules required by TensorFlow 2.8 (no CUDA in conda env)
ml cuda/11.7.0
ml cudnn/8.9.5-11
ml proxies
ml tensorrt/8.5.3.1

source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh
set +u; conda activate chrombpnet; set -u
export MPLBACKEND=Agg

MESO_DIR=/sc/arion/scratch/giottb01/Xmen/meso
GREFDIR=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet
OUT_DIR=/sc/arion/scratch/giottb01/Xmen/meso/footprints_pattern24
MOTIF_FILE=${OUT_DIR}/motif_footprints_pattern24.txt

# Use KLRC1_NK peaks as the common region set for all cell types
REGIONS=${MESO_DIR}/KLRC1_NK/KLRC1_NK_peakset_all_no_blacklist.bed

CELLTYPES=(CD4 CD8 CD8_exhausted FGFBP2_NK KLRC1_NK Tregs)

mkdir -p "${OUT_DIR}"

# ── Run chrombpnet footprints per celltype per fold ──────────────────────────
for ct in "${CELLTYPES[@]}"; do
    echo "=== ${ct} ==="
    for fold in 0 1 2 3 4; do
        MODEL=${MESO_DIR}/${ct}/no_bias_model/fold_${fold}/models/chrombpnet_nobias.h5
        OUT_PREFIX=${OUT_DIR}/${ct}_fold${fold}

        if [ -f "${OUT_PREFIX}_footprints.h5" ]; then
            echo "  fold ${fold}: already done, skipping"
            continue
        fi

        echo "  fold ${fold}: running footprints..."
        chrombpnet footprints \
            -m  "${MODEL}" \
            -r  "${REGIONS}" \
            -g  "${GREFDIR}/genome_references/hg38.genome.fa" \
            -fl "${GREFDIR}/folds/fold_${fold}.json" \
            -op "${OUT_PREFIX}" \
            -pwm_f "${MOTIF_FILE}"

        echo "  fold ${fold}: done"
    done
done

export HDF5_PLUGIN_PATH=${CONDA_PREFIX}/lib/hdf5/plugin

# ── Average folds + plot ──────────────────────────────────────────────────────
PLOT_SCRIPT=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/plot_pattern24_footprints.py

echo "=== Averaging folds and plotting ==="
python3 "${PLOT_SCRIPT}" --out_dir "${OUT_DIR}"

echo "=== Done. ==="
