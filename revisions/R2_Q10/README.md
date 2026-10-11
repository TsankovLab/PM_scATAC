# R2_Q10 — keratin expression in epithelioid tumour cells vs normal mesothelium (scRNA-seq)

- `R2_Q10_keratins_epithelioid_vs_mesothelium.R` — malignant and mesothelium cells of `main/scrna/srt.rds`; tumours are called
  epithelioid from the scS-score (cNMF20 of the 25-programme malignant cNMF, top 20 genes, AddModuleScore; median over the
  malignant cells of each tumour, cut at the first clear break of the ordered medians, 0.18). Keratins detected in >= 5% of the
  cells of a group are compared between epithelioid tumour cells and mesothelium at cell level and as per-patient pseudobulk.
- `run_step.sh` — runs a script of this folder in the `meso_scatac` environment.

Epithelioid tumours by this rule: P8, P11, P4, P5, P12. Higher scS: P3, P13, P14, P1.
The 318 mesothelium cells come from two patients only (P1: 286, P13: 32), neither of them an epithelioid tumour, so the
comparison is between patients; with 5 vs 2 patients the smallest possible sample-level Wilcoxon p is 0.095.

Outputs: `R2_Q10_scS_per_sample.csv`, `R2_Q10_keratin_stats.csv`, `R2_Q10_keratin_pseudobulk_per_sample.csv`,
`Plots/R2_Q10_keratins_epithelioid_vs_mesothelium.pdf`.
