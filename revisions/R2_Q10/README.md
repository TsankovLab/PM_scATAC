# R2_Q10 — keratin gene scores in epithelioid tumour cells vs normal mesothelium

Run in this order (`run_step.sh <script>` runs a script of this folder in the `meso_scatac` environment):

1. `R2_Q10_keratins_epithelioid_vs_mesothelium.R` — scRNA side (`main/scrna/srt.rds`): tumours called epithelioid from the
   scS-score (cNMF20 of the 25-programme malignant cNMF), keratins compared between epithelioid tumour cells and mesothelium
   at cell level and per-patient pseudobulk. Writes `R2_Q10_scS_per_sample.csv` and `R2_Q10_keratin_stats.csv`, used by step 2.
2. `R2_Q10_keratins_scATAC_genescore.R` — the same on the ArchR gene scores of `main/scatac_ArchR` (scS from gene scores,
   cut 0 on the tumour median), plus the agreement with the scRNA result.
   Output: `Plots/R2_Q10_keratins_scATAC_genescore.pdf`, `R2_Q10_scATAC_*.csv`, `R2_Q10_keratin_RNA_vs_ATAC.csv`.

Note: the mesothelium cells come from two patients only (P1, P13), neither of them epithelioid, so the comparison is between
patients.
