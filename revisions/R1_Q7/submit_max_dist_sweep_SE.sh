#!/bin/bash
#BSUB -J max_dist_SE
#BSUB -P <LSF_ACCOUNT>
#BSUB -q premium
#BSUB -n 8
#BSUB -R "rusage[mem=240000] span[hosts=1]"
#BSUB -W 24:00
#BSUB -o logs/max_dist_SE_%J.out
#BSUB -e logs/max_dist_SE_%J.err
source "${MESO_CONFIG:-$(for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$PWD"; do while [ ! -f "$d/config.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; [ -f "$d/config.sh" ] && { echo "$d/config.sh"; break; }; done)}"   # data / results / reference roots: see config.sh at the repo root
mkdir -p ${RESULTS_DIR}/logs
module purge
${CONDA_ENVS}/meso_scatac/bin/Rscript ${REPO_DIR}/revisions/R1_Q7/max_dist_sweep_SE.R
