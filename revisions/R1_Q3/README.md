# R1_Q3 — SOX9 / RUNX2 / SNAI2 and survival, with and without scS adjustment

- `R1_Q3_F30b_survival_forest_compact.R` — Cox proportional hazards models of overall survival on SOX9, RUNX2 and SNAI2
  expression in the three bulk cohorts (TCGA, Bueno, MESOMICS; `bulkRNA_meso/bulk_RNA_studies*.rds`), alone and adjusted
  for the sarcomatoid score (scS). HR per 1 SD of expression; the median split is computed and kept in the table only.

Output: `Plots/F30b_survival_forest_compact.pdf`, `F30b_survival_forest_table.csv`.
