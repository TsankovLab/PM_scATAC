#!/bin/bash
set -euo pipefail
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WD=/sc/arion/scratch/giottb01/meso_R1Q8; mkdir -p "$WD"
bsub -J q8_TFscatter -P acc_Tsankov_Normal_Lung -q premium -n 4 -W 3:00 \
  -R "rusage[mem=64000]" -R "span[hosts=1]" -o "$WD/o.out" -e "$WD/o.err" \
  /bin/bash -c "export OPENBLAS_NUM_THREADS=4; source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh; conda activate meso_scatac; Rscript '$SD/R1_Q8_discover_oncogenic_TFs_errorbars.R'"
echo "WD=$WD"
