#!/bin/bash
set -euo pipefail
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WD=/sc/arion/scratch/giottb01/meso_allTF_cis; mkdir -p "$WD"
bsub -J alltf_cis -P acc_Tsankov_Normal_Lung -q premium -n 4 -W 3:00 \
  -R "rusage[mem=32000]" -R "span[hosts=1]" -o "$WD/a.out" -e "$WD/a.err" \
  /bin/bash -c "source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh; conda activate meso_scatac; Rscript '$SD/R1_Q3_allTF_cis_enrichment_SOX9.R'"
