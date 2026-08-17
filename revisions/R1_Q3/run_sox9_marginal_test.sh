#!/bin/bash
# run_sox9_marginal_test.sh
# Early "is SOX9 learned?" test: marginal footprints of the SOX9 motif panel on
# the SOX9pos and SOX9neg ChromBPNet models (first fold trained in BOTH), using
# only the trained models — no contribution scores / modisco needed.
#
# Submit:  bsub < run_sox9_marginal_test.sh   (or via the monitor)
#BSUB -J sox9_marginal
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q gpu
#BSUB -n 4
#BSUB -W 4:00
#BSUB -gpu "num=1"
#BSUB -R "a100"
#BSUB -R "rusage[mem=32000]"
#BSUB -R "span[hosts=1]"
#BSUB -o /sc/arion/scratch/giottb01/meso_SOX9_module_chromBPnet/sox9_marginal_%J.out
#BSUB -e /sc/arion/scratch/giottb01/meso_SOX9_module_chromBPnet/sox9_marginal_%J.err

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
PLOT=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3/plot_sox9_marginal.py
mkdir -p "${OUT_DIR}"

# find the lowest fold trained in BOTH cell types
FOLD=""
for f in 0 1 2 3 4; do
    if [ -f "${WORKDIR}/SOX9pos/no_bias_model/fold_${f}/models/chrombpnet_nobias.h5" ] && \
       [ -f "${WORKDIR}/SOX9neg/no_bias_model/fold_${f}/models/chrombpnet_nobias.h5" ]; then
        FOLD="${f}"; break
    fi
done
[ -z "${FOLD}" ] && echo "ERROR: no fold trained in both cell types yet" && exit 1
echo "=== Using fold ${FOLD} ==="

for ct in SOX9pos SOX9neg; do
    OP="${OUT_DIR}/${ct}_fold${FOLD}"
    if [ ! -f "${OP}_footprints.h5" ]; then
        echo "=== marginal footprints: ${ct} fold ${FOLD} ==="
        chrombpnet footprints \
            -m  "${WORKDIR}/${ct}/no_bias_model/fold_${FOLD}/models/chrombpnet_nobias.h5" \
            -r  "${WORKDIR}/${ct}/${ct}_peakset_all_no_blacklist.bed" \
            -g  "${GREFDIR}/genome_references/hg38.genome.fa" \
            -fl "${GREFDIR}/folds/fold_${FOLD}.json" \
            -op "${OP}" \
            -pwm_f "${MOTIF_FILE}"
    fi
done

export HDF5_PLUGIN_PATH="${CONDA_PREFIX}/lib/hdf5/plugin"
python3 "${PLOT}" \
    --pos "${OUT_DIR}/SOX9pos_fold${FOLD}_footprints.h5" \
    --neg "${OUT_DIR}/SOX9neg_fold${FOLD}_footprints.h5" \
    --out_pdf "${OUT_DIR}/sox9_marginal_pos_vs_neg.pdf" \
    --out_csv "${OUT_DIR}/sox9_marginal_comparison.csv" \
    --fold "${FOLD}"
echo "=== SOX9 marginal test complete (fold ${FOLD}) ==="
