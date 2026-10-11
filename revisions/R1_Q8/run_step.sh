#!/bin/bash
# Run one R1_Q8 script in the meso_scatac environment, with the working directory set to
# R1_Q8 itself -- R2_Q5/run_step.sh cds into R2_Q5, so scripts using relative output paths
# silently wrote there.
#   ./run_step.sh R1_Q8_sarc_TF_boxplots_paired.R
set -euo pipefail
source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh
conda activate meso_scatac
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$(dirname "$0")"
Rscript "$@"
