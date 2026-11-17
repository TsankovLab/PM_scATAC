#!/bin/bash
#BSUB -J add_myeloid_col_R1Q4
#BSUB -P <LSF_ACCOUNT>
#BSUB -q premium
#BSUB -n 1
#BSUB -R "rusage[mem=32000]"
#BSUB -W 1:00
#BSUB -o logs/add_myeloid_col_R1Q4_%J.out
#BSUB -e logs/add_myeloid_col_R1Q4_%J.err

source "${MESO_CONFIG:-$(for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$PWD"; do while [ ! -f "$d/config.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; [ -f "$d/config.sh" ] && { echo "$d/config.sh"; break; }; done)}"   # data / results / reference roots: see config.sh at the repo root
export OPENBLAS_NUM_THREADS=1

${CONDA_ENVS}/meso_scatac/bin/Rscript \
  ${REPO_DIR}/revisions/R1_Q4/add_myeloid_compartment_col.R
