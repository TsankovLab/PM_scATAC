#!/bin/bash
# Run one pipeline script in the meso_scatac environment.
#   ./run_step.sh 03_clone_calls.R
#   ./run_step.sh 02_run_epianeufinder.R P4
set -euo pipefail
source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh
conda activate meso_scatac
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$(dirname "$0")"
Rscript "$@"
