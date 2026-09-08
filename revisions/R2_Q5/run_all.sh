#!/bin/bash
###############################################################################
# Reproduce R2_Q5.
#
# Step 2 is the expensive one (cross-modality label transfer, ~1-3 h, ~120 GB);
# submit it to LSF and wait for atac_label_transfer.csv before running steps 3-6:
#     ./submit_02_label_transfer.sh
#
# Steps 1 and 4 read the ArchR project and the Seurat object and need a few GB;
# steps 3, 5 and 6 are seconds and read only the CSVs in this folder.
###############################################################################
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
for s in 01_annotations.R \
         03_concordance.R \
         04_pseudobulk_concordance.R \
         05_composition.R \
         07_malignant_disagreement.R \
         09_figure_joint.R \
         06_figures.R ; do
  echo "=================================================================="
  echo ">>> $s   $(date '+%H:%M:%S')"
  echo "=================================================================="
  ./run_step.sh "$s" 2>&1 | tee "logs/${s%.R}.log"
done
echo; echo "ALL DONE  $(date '+%H:%M:%S')"
