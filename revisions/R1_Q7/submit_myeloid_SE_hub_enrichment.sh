#!/bin/bash
#BSUB -J myeloid_SE_hub_enrich
#BSUB -P <LSF_ACCOUNT>
#BSUB -q premium
#BSUB -n 8
#BSUB -R "rusage[mem=128000] span[hosts=1]"
#BSUB -W 12:00
#BSUB -o logs/myeloid_SE_hub_enrich_%J.out
#BSUB -e logs/myeloid_SE_hub_enrich_%J.err

source "${MESO_CONFIG:-$(for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$PWD"; do while [ ! -f "$d/config.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; [ -f "$d/config.sh" ] && { echo "$d/config.sh"; break; }; done)}"   # data / results / reference roots: see config.sh at the repo root
mkdir -p ${RESULTS_DIR}/logs

module purge

${CONDA_ENVS}/meso_scatac/bin/Rscript \
  ${REPO_DIR}/revisions/R1_Q7/myeloid_SE_hub_enrichment.R
