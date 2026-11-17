#!/bin/bash
#BSUB -J celltype_barplots
#BSUB -P <LSF_ACCOUNT>
#BSUB -q premium
#BSUB -n 2
#BSUB -R "rusage[mem=32000]"
#BSUB -W 2:00
#BSUB -o logs/celltype_barplots_%J.out
#BSUB -e logs/celltype_barplots_%J.err

source "${MESO_CONFIG:-$(for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$PWD"; do while [ ! -f "$d/config.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; [ -f "$d/config.sh" ] && { echo "$d/config.sh"; break; }; done)}"   # data / results / reference roots: see config.sh at the repo root
mkdir -p ${RESULTS_DIR}/logs

module purge

${CONDA_ENVS}/meso_scatac/bin/Rscript \
  ${REPO_DIR}/revisions/R2_Q2/celltype_barplots.R

${CONDA_ENVS}/meso_scatac/bin/Rscript \
  ${REPO_DIR}/revisions/R2_Q2/celltype_barplots_axisbreak.R
