#!/bin/bash
# Runs the CROP-seq validation figures end to end. Output folder: $CROP_OUT (default: <perturbation_panel>/figs_repro).
# Needs the module-score outputs screen_cm_diff_matrix.csv and cm_sample_scores_AddModuleScore_bulk_RUN1.csv (see profiles.R).
set -euo pipefail
module load R/4.2.0
cd "$(dirname "$0")"
for s in screen_export_figs.R bulk_cellcycle_RUN1.R fig1_variable_genes_heatmaps.R fig2_proliferation.R \
         fig3_KO_KO_correlation.R fig4_bulk_vs_CROPseq_correlation.R fig5_bestreplicate_correlation.R; do
  echo "=== $s ==="; Rscript "$s"
done
