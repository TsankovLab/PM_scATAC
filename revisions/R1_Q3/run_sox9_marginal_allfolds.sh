#!/bin/bash
# run_sox9_marginal_allfolds.sh
# Extend the "is SOX9 learned?" marginal-footprint test to ALL trained folds of
# the SOX9pos and SOX9neg ChromBPNet models, then aggregate (mean +/- SD across
# folds, paired t-test pos vs neg). Uses only the trained models — no contrib
# scores / modisco needed. Fold 0 footprints already exist and are reused.
#
# Submit:  bsub < run_sox9_marginal_allfolds.sh
#BSUB -J sox9_marginal_allfolds
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q gpu
#BSUB -n 4
#BSUB -W 8:00
#BSUB -gpu "num=1"
#BSUB -R "a100"
#BSUB -R "rusage[mem=32000]"
#BSUB -R "span[hosts=1]"
#BSUB -o /sc/arion/scratch/giottb01/meso_SOX9_module_chromBPnet/sox9_marginal_allfolds_%J.out
#BSUB -e /sc/arion/scratch/giottb01/meso_SOX9_module_chromBPnet/sox9_marginal_allfolds_%J.err

set -euo pipefail
ml cuda/11.7.0 2>/dev/null || true
ml cudnn/8.9.5-11 2>/dev/null || true
ml proxies 2>/dev/null || true
ml tensorrt/8.5.3.1 2>/dev/null || true
source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh
set +u; conda activate chrombpnet; set -u
export MPLBACKEND=Agg

WORKDIR=/sc/arion/scratch/giottb01/meso_SOX9_module_chromBPnet
GREFDIR=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet
MOTIF_FILE=${WORKDIR}/SOX_motifs_footprint.txt
OUT_DIR=${WORKDIR}/sox9_marginal_test
PLOT=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3/plot_sox9_marginal_allfolds.py
mkdir -p "${OUT_DIR}"

# compute marginal footprints for every fold trained in BOTH cell types
DONE_FOLDS=()
for f in 0 1 2 3 4; do
    POS_M="${WORKDIR}/SOX9pos/no_bias_model/fold_${f}/models/chrombpnet_nobias.h5"
    NEG_M="${WORKDIR}/SOX9neg/no_bias_model/fold_${f}/models/chrombpnet_nobias.h5"
    if [ ! -f "${POS_M}" ] || [ ! -f "${NEG_M}" ]; then
        echo "=== skip fold ${f}: model not trained in both cell types ==="
        continue
    fi
    for ct in SOX9pos SOX9neg; do
        OP="${OUT_DIR}/${ct}_fold${f}"
        if [ -f "${OP}_footprints.h5" ]; then
            echo "=== reuse existing footprints: ${ct} fold ${f} ==="
        else
            echo "=== marginal footprints: ${ct} fold ${f} ==="
            chrombpnet footprints \
                -m  "${WORKDIR}/${ct}/no_bias_model/fold_${f}/models/chrombpnet_nobias.h5" \
                -r  "${WORKDIR}/${ct}/${ct}_peakset_all_no_blacklist.bed" \
                -g  "${GREFDIR}/genome_references/hg38.genome.fa" \
                -fl "${GREFDIR}/folds/fold_${f}.json" \
                -op "${OP}" \
                -pwm_f "${MOTIF_FILE}"
        fi
    done
    DONE_FOLDS+=("${f}")
done

echo "=== Folds processed: ${DONE_FOLDS[*]} ==="

export HDF5_PLUGIN_PATH="${CONDA_PREFIX}/lib/hdf5/plugin"
python3 "${PLOT}" \
    --dir     "${OUT_DIR}" \
    --out_pdf "${OUT_DIR}/sox9_marginal_pos_vs_neg_allfolds.pdf" \
    --out_csv "${OUT_DIR}/sox9_marginal_comparison_allfolds.csv"

echo "=== SOX9 all-folds marginal test complete ==="
