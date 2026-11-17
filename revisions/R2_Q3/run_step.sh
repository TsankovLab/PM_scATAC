#!/bin/bash
# Run one pipeline script in the meso_scatac environment.
#   ./run_step.sh 03_clone_calls.R
#   ./run_step.sh 02_run_epianeufinder.R P4
source "${MESO_CONFIG:-$(for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$PWD"; do while [ ! -f "$d/config.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; [ -f "$d/config.sh" ] && { echo "$d/config.sh"; break; }; done)}"   # data / results / reference roots: see config.sh at the repo root
set -euo pipefail
source ${CONDA_INIT}
conda activate meso_scatac
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$(dirname "$0")"
Rscript "$@"
