#!/bin/bash
# Run one R2_Q10 script in the meso_scatac environment, with the working directory set to
# R2_Q10 itself -- R2_Q5/run_step.sh cds into R2_Q5, so scripts using relative output paths
# silently wrote there.
#   ./run_step.sh R2_Q10_sarc_TF_boxplots_paired.R
source "${MESO_CONFIG:-$(for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$PWD"; do while [ ! -f "$d/config.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; [ -f "$d/config.sh" ] && { echo "$d/config.sh"; break; }; done)}"   # data / results / reference roots: see config.sh at the repo root
set -euo pipefail
source ${CONDA_INIT}
conda activate meso_scatac
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$(dirname "$0")"
Rscript "$@"
