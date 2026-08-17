#!/bin/bash
set -euo pipefail
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WD=/sc/arion/scratch/giottb01/meso_RUNXdimer_sarc; mkdir -p "$WD"
bsub -J runxdimer_sarc -P acc_Tsankov_Normal_Lung -q premium -n 4 -W 2:00 \
  -R "rusage[mem=48000]" -R "span[hosts=1]" -o "$WD/s.out" -e "$WD/s.err" \
  /bin/bash -c "export OPENBLAS_NUM_THREADS=4; source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh; conda activate meso_scatac; Rscript '$SD/R1_Q3_RUNXdimer_P2G_vs_sarc.R'"
