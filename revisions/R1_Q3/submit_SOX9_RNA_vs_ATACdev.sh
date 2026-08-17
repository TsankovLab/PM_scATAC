#!/bin/bash
# submit_SOX9_RNA_vs_ATACdev.sh — cross-patient SOX9 scRNA expr vs ATAC chromVAR deviation
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR=/sc/arion/scratch/giottb01/meso_SOX9_RNA_vs_ATAC
mkdir -p "${WORKDIR}"
bsub \
  -J sox9_rna_vs_atac \
  -P acc_Tsankov_Normal_Lung \
  -q premium \
  -n 4 -W 4:00 \
  -R "rusage[mem=24000]" -R "span[hosts=1]" \
  -o "${WORKDIR}/rna_vs_atac.out" \
  -e "${WORKDIR}/rna_vs_atac.err" \
  /bin/bash -c "
source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh
conda activate meso_scatac
Rscript '${SCRIPTDIR}/R1_Q3_SOX9_RNA_vs_ATACdeviation_persample.R'
"
echo "Submitted. logs: ${WORKDIR}/rna_vs_atac.{out,err}"
