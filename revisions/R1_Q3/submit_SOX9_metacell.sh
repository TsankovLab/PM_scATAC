#!/bin/bash
set -euo pipefail
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WD=/sc/arion/scratch/giottb01/meso_SOX9_metacell; mkdir -p "$WD"
bsub -J sox9_metacell -P acc_Tsankov_Normal_Lung -q premium -n 4 -W 3:00 \
  -R "rusage[mem=24000]" -R "span[hosts=1]" -o "$WD/mc.out" -e "$WD/mc.err" \
  /bin/bash -c "source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh; conda activate meso_scatac; Rscript '$SD/R1_Q3_SOX9_metacell_expr_vs_module.R'"
