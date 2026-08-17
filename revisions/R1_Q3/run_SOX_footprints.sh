#!/bin/bash
# run_SOX_footprints.sh  (R1_Q3)
#
# ChromBPNet footprints for SOX/HMG motifs (de novo + canonical SOX9_MA0077.1)
# in both SOX9_high_P23 and SOX9_low_P23 models, 5 folds each.
#
# Submit: bsub < run_SOX_footprints.sh
#BSUB -J SOX_footprints
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q gpu
#BSUB -n 4
#BSUB -W 12:00
#BSUB -gpu "num=1"
#BSUB -R "a100"
#BSUB -R "rusage[mem=32000]"
#BSUB -R "span[hosts=1]"
#BSUB -o /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/SOX_footprints_%J.out
#BSUB -e /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/SOX_footprints_%J.err

set -euo pipefail

ml cuda/11.7.0
ml cudnn/8.9.5-11
ml proxies
ml tensorrt/8.5.3.1

source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh
set +u; conda activate chrombpnet; set -u
export MPLBACKEND=Agg

MESO_DIR=/sc/arion/scratch/giottb01/meso_SOX9_P23_chromBPnet
GREFDIR=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet
MOTIF_FILE=${MESO_DIR}/SOX_motifs_footprint.txt
OUT_DIR=${MESO_DIR}/footprints_SOX_motifs
PLOT_SCRIPT=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3/plot_SOX_footprints.py

CELLTYPES=(SOX9_high_P23 SOX9_low_P23)

mkdir -p "${OUT_DIR}"

for ct in "${CELLTYPES[@]}"; do
    echo "=== ${ct} ==="
    REGIONS=${MESO_DIR}/${ct}/${ct}_peakset_all_no_blacklist.bed
    for fold in 0 1 2 3 4; do
        OUT_PREFIX=${OUT_DIR}/${ct}_fold${fold}
        if [ -f "${OUT_PREFIX}_footprints.h5" ]; then
            echo "  fold ${fold}: already done, skipping"
            continue
        fi
        echo "  fold ${fold}: running..."
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

echo "=== Averaging folds and plotting ==="
python3 "${PLOT_SCRIPT}" --out_dir "${OUT_DIR}"

echo "=== Done ==="
