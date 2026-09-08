#!/bin/bash
# Run one R2_Q5 script in the meso_scatac environment.
#   ./run_step.sh 01_annotations.R
set -euo pipefail
source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh
conda activate meso_scatac
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$(dirname "$0")"
Rscript "$@"
