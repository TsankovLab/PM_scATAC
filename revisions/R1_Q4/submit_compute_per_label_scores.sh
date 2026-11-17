#!/bin/bash
#BSUB -J per_label_scores_R1Q4
#BSUB -P <LSF_ACCOUNT>
#BSUB -q premium
#BSUB -n 1
#BSUB -R "rusage[mem=64000]"
#BSUB -W 2:00
#BSUB -o logs/per_label_scores_R1Q4_%J.out
#BSUB -e logs/per_label_scores_R1Q4_%J.err

source "${MESO_CONFIG:-$(for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$PWD"; do while [ ! -f "$d/config.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; [ -f "$d/config.sh" ] && { echo "$d/config.sh"; break; }; done)}"   # data / results / reference roots: see config.sh at the repo root
export OPENBLAS_NUM_THREADS=1

${CONDA_ENVS}/meso_scatac/bin/Rscript \
  ${REPO_DIR}/revisions/R1_Q4/compute_per_label_scores.R
