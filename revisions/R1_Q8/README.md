# R1_Q8 — tumour vs normal TF activity and expression, with error bars

- `R1_Q8_discover_oncogenic_TFs_errorbars.R` (`submit_discover_oncogenic_TFs.sh`) — re-implements the TF selection of
  `git_repo/tumor_analysis/discover_oncogenic_TFs.R` (chromVAR deviation per sample, tumour vs normal1, presto Wilcoxon per
  tumour, >= 5 tumours significant, expressed in the tumour scRNA) and draws the activity-difference vs RNA-difference scatter
  with +/- SEM across tumour samples on the labelled TFs.

Output: `Plots/Diff_normal_tumor_deviation_and_rna_scatterplot5_errorbars.pdf` (5 x 4 in) and `..._errorbars_large.pdf`,
`TF_diff_rna_with_errorbars.csv`, `selected_TF_table_with_errorbars.csv`.
