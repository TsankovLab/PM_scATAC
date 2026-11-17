#!/bin/bash
#BSUB -J modisco_merge
#BSUB -P <LSF_ACCOUNT>
#BSUB -q gpu
#BSUB -n 8
#BSUB -W 48:00
#BSUB -gpu num=1
#BSUB -R v100
#BSUB -R rusage[mem=8000]
#BSUB -R span[hosts=1]
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

source "${MESO_CONFIG:-$(for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$PWD"; do while [ ! -f "$d/config.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; [ -f "$d/config.sh" ] && { echo "$d/config.sh"; break; }; done)}"   # data / results / reference roots: see config.sh at the repo root
ml anaconda3/2022.10

output_dir=${1}
model_head=${2}
cluster_key=${3}
modisco_dir=${4}
contribs_dir=${5}
batch=${6}

source activate chrombpnet
mkdir -p ${output_dir}/modisco_merged_${model_head}/merged_motifs

PYTHONPATH=${DATA_DIR}/HDMA/code/03-chrombpnet \
${CONDA_ENVS}/chrombpnet/bin/python -u ${REPO_DIR}/utils/03-merge_modisco.py \
                    --out-dir ${output_dir}/modisco_merged_${model_head}/merged_motifs \
                    --model-head ${model_head} \
                    --cluster-key ${cluster_key} \
                    --modisco-dir ${modisco_dir} \
                    --contribs-dir ${contribs_dir} \
                    --batch ${batch}
