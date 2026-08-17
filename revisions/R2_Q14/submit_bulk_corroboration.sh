#!/bin/bash
set -euo pipefail
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WD=/sc/arion/scratch/giottb01/meso_R2Q14; mkdir -p "$WD"
bsub -J q14_bulk -P acc_Tsankov_Normal_Lung -q premium -n 2 -W 1:00 \
  -R "rusage[mem=32000]" -R "span[hosts=1]" -o "$WD/bulk.out" -e "$WD/bulk.err" \
  /bin/bash -c "export OPENBLAS_NUM_THREADS=2; source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh; conda activate meso_scatac; Rscript '$SD/R2_Q14_bulk_corroboration.R'"
echo "WD=$WD"
