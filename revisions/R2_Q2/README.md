# R2_Q2 — cell type composition per sample

`submit_celltype_barplots.sh` runs both scripts:

- `celltype_barplots.R` — main ArchR project: proportions per sample (as in the main analysis) and raw cell counts per
  sample with the per-sample total. Output `celltype_barplots.pdf`.
- `celltype_barplots_axisbreak.R` — the same figure with a broken count axis, because P23 (18,228 of 49,849 cells) squeezes
  the other samples into a quarter of the panel. Cell counts are cached to `celltype_counts.csv` on the first run.
  Output `celltype_barplots_axisbreak.pdf`.
