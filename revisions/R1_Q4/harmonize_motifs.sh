#!/bin/bash
# harmonize_motifs.sh
#
# Harmonizes ChromBPNet motifs across momac, resident, and PBMC_CD14_Mono
# by clustering and merging modisco patterns following the HDMA pipeline:
#
#   1. Convert per-cell-type modisco patterns to PFMs  (modisco_to_pfm.py)
#   2. Cluster PFMs across cell types                  (gimme cluster)
#   3. Merge modisco seqlets within each cluster       (03-merge_modisco.py, batched via LSF)
#   4. Compile all merged clusters into one h5/tsv     (04a-compile_modisco_obj.py)
#   5. Match to known TF motifs via TOMTOM             (04b-get_tomtom_matches.py)
#   6. Re-run FI-NeMo with compiled merged motifs      (finemo_motif_calls_on_merged.sh)
#
# Run interactively on login node (steps 1-2, 4-5 are fast; step 3 submits LSF jobs).
# Usage: bash harmonize_motifs.sh [profile|counts]

set -euo pipefail

MODEL_HEAD="${1:-profile}"

CHROMBPDIR=/sc/arion/scratch/giottb01/scatac_meso/revisions/R1_Q4
REPODIR=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo
HDMA_UTILS=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/HDMA/code/03-chrombpnet
MOTIF_DB=/sc/arion/work/giottb01/conda/envs/chrombpnet/lib/python3.8/site-packages/chrombpnet/data/motifs.meme.txt

CELLTYPES=(momac resident PBMC_CD14_Mono)
MERGED_DIR=${CHROMBPDIR}/modisco_merged_${MODEL_HEAD}

PYTHON=/sc/arion/work/giottb01/conda/envs/chrombpnet/bin/python

ml anaconda3/2022.10
source activate chrombpnet

echo "=== Harmonizing motifs | model_head=${MODEL_HEAD} | cell types: ${CELLTYPES[*]} ==="
mkdir -p ${MERGED_DIR}

# ---- 1. Build modisco h5 paths TSV ------------------------------------------
PATHS_TSV=${CHROMBPDIR}/modisco_${MODEL_HEAD}_h5_paths.tsv
: > ${PATHS_TSV}
for ct in "${CELLTYPES[@]}"; do
    h5=${CHROMBPDIR}/${ct}/no_bias_model/modisco_${MODEL_HEAD}/modisco_results_${MODEL_HEAD}.h5
    if [ ! -f "${h5}" ]; then
        echo "ERROR: missing modisco h5 for ${ct}: ${h5}"
        exit 1
    fi
    echo -e "${ct}\t${h5}" >> ${PATHS_TSV}
done
echo "Paths TSV written: ${PATHS_TSV}"

# ---- 2. Convert modisco patterns to PFMs ------------------------------------
echo "=== Step 1: modisco -> PFMs ==="
PYTHONPATH=${HDMA_UTILS} ${PYTHON} ${REPODIR}/utils/modisco_to_pfm.py \
    -c ${PATHS_TSV} \
    -o ${MERGED_DIR}/pfms_pos.txt \
    -p pos_patterns \
    -t 0.3 -ml 3

PYTHONPATH=${HDMA_UTILS} ${PYTHON} ${REPODIR}/utils/modisco_to_pfm.py \
    -c ${PATHS_TSV} \
    -o ${MERGED_DIR}/pfms_neg.txt \
    -p neg_patterns \
    -t 0.3 -ml 3

# ---- 3. Cluster PFMs with gimme ---------------------------------------------
echo "=== Step 2: gimme cluster (pos + neg) ==="
gimme cluster ${MERGED_DIR}/pfms_pos.txt ${MERGED_DIR}/cluster_pos -t 0.8 -N 8
gimme cluster ${MERGED_DIR}/pfms_neg.txt ${MERGED_DIR}/cluster_neg -t 0.8 -N 8

# Build combined cluster key with batch assignments (every 3 rows = 1 batch)
awk -F'\t' 'BEGIN{OFS="\t"} {n=int((NR-1)/3)+1; print $1,"pos_patterns",$2,n}' \
    ${MERGED_DIR}/cluster_pos/cluster_key.txt \
    > ${MERGED_DIR}/cluster_pos/cluster_key_modified_pos.txt

awk -F'\t' 'BEGIN{OFS="\t"} {n=int((NR-1)/3)+1; print $1,"neg_patterns",$2,n}' \
    ${MERGED_DIR}/cluster_neg/cluster_key.txt \
    > ${MERGED_DIR}/cluster_neg/cluster_key_modified_neg.txt

cat ${MERGED_DIR}/cluster_pos/cluster_key_modified_pos.txt \
    ${MERGED_DIR}/cluster_neg/cluster_key_modified_neg.txt \
    > ${MERGED_DIR}/cluster_key_combined.txt

CLUSTER_KEY=${MERGED_DIR}/cluster_key_combined.txt
echo "Cluster key: $(wc -l < ${CLUSTER_KEY}) entries"

# ---- 4. Submit merge jobs (one per batch) -----------------------------------
echo "=== Step 3: submitting modisco merge jobs ==="
mkdir -p ${MERGED_DIR}/merged_motifs

BATCHES=$(cut -f4 ${CLUSTER_KEY} | sort -n | uniq)
MERGE_JOB_IDS=""

for batch in ${BATCHES}; do
    job_id=$(bsub \
        -J modisco_merge_b${batch} \
        -P acc_Tsankov_Normal_Lung \
        -q premium \
        -n 8 -W 12:00 \
        -R "rusage[mem=80000]" -R "span[hosts=1]" \
        -o ${MERGED_DIR}/modisco_merge_${batch}.out \
        -e ${MERGED_DIR}/modisco_merge_${batch}.err \
        ${REPODIR}/utils/modisco_merge_job.sh \
            "${CHROMBPDIR}" "${MODEL_HEAD}" "${CLUSTER_KEY}" \
            "${CHROMBPDIR}" "${CHROMBPDIR}" "${batch}" \
        | awk '{print $2}' | tr -d '<>')
    echo "  Batch ${batch}: job ${job_id}"
    MERGE_JOB_IDS="${MERGE_JOB_IDS} done(${job_id}) &&"
done

# Strip trailing &&
MERGE_DEP="${MERGE_JOB_IDS% &&}"

# ---- 5. Submit compile + tomtom + finemo as post-merge job ------------------
echo "=== Submitting post-merge job (compile + tomtom + finemo) ==="

bsub \
    -J modisco_postmerge_${MODEL_HEAD} \
    -P acc_Tsankov_Normal_Lung \
    -q premium \
    -n 4 -W 24:00 \
    -R "rusage[mem=32000]" -R "span[hosts=1]" \
    -w "${MERGE_DEP}" \
    -o ${MERGED_DIR}/postmerge.out \
    -e ${MERGED_DIR}/postmerge.err \
    /bin/bash -c "
set -euo pipefail
ml anaconda3/2022.10
source activate chrombpnet

PYTHON=${PYTHON}
PYTHONPATH=${HDMA_UTILS}
export PYTHONPATH

MERGED_DIR=${MERGED_DIR}
CLUSTER_KEY=${CLUSTER_KEY}
CHROMBPDIR=${CHROMBPDIR}
MODEL_HEAD=${MODEL_HEAD}
REPODIR=${REPODIR}
MOTIF_DB=${MOTIF_DB}

# 4a. Compile all merged h5s into one
echo '=== Compile merged modisco ==='
mkdir -p \${MERGED_DIR}/compiled
\${PYTHON} \${REPODIR}/utils/04a-compile_modisco_obj.py \
    --output_dir \${MERGED_DIR} \
    --cluster_key_path \${CLUSTER_KEY} \
    --compiled_modisco_h5_path \${MERGED_DIR}/compiled/modisco_compiled.h5 \
    --compiled_modisco_tsv_path \${MERGED_DIR}/compiled/modisco_compiled.tsv \
    --modisco_merged_dir \${MERGED_DIR}/merged_motifs

# 4b. TOMTOM matching
echo '=== TOMTOM matching ==='
mkdir -p \${MERGED_DIR}/compiled_tomtom
\${PYTHON} \${REPODIR}/utils/04b-get_tomtom_matches.py \
    --modisco-h5 \${MERGED_DIR}/compiled/modisco_compiled.h5 \
    --out-dir \${MERGED_DIR}/compiled_tomtom \
    --meme-db \${MOTIF_DB} \
    --verbose True

# 5. Re-run FI-NeMo per cell type with compiled merged motifs
echo '=== FI-NeMo on merged motifs ==='
source activate finemo
for ct in ${CELLTYPES[*]}; do
    echo \"  finemo: \${ct}\"
    bash \${REPODIR}/utils/finemo_motif_calls_on_merged.sh \
        \${CHROMBPDIR} \${ct} \${MODEL_HEAD}
done

echo '=== Post-merge pipeline complete ==='
"

echo ""
echo "=== All jobs submitted. Monitor with: bjobs -u giottb01 ==="
echo "    Merge jobs will complete first, then post-merge (compile+tomtom+finemo) runs."
