# R2_Q9 — improved sample palette and regenerated per-sample figures

The reviewer asked for a more immediate sample colour scheme. The old `palette_sample` (`git_repo/utils/palettes.R`) used ten muted
brown/grey/blue tones for the tumours and the same red for every normal sample.

## New palette (`R2_Q9_palette_sample.R`, drop-in replacement)
- **Tumours:** 11 clearly different hues, assigned along the sarcomatoid score (median scS of malignant cells, scATAC):
  cool colours (green, teal, blue) = epithelioid-like, violet/purple = intermediate, warm colours (pink, amber, red) = higher scS.
- **Normal samples:** black/greys, one shade each; `normal1` (the one plotted with the tumours) is near-black.
- **P11_HOX:** light tint of the P11 colour.
- Smallest pairwise colour distance among the 11 tumours (OKLab x 100): 18.3 for normal vision and 10.1 under simulated red-green
  colour blindness, against 5.0 and 2.1 for the old palette.
- The mid greys of normal2/normal3/HU37/HU62 can come close to P10/P12 under colour-blind simulation; these samples are not plotted
  together with the tumours in any scATAC figure.

## Scripts (plotting code copied from the originals, only the palette changes; projects are read-only)
`R2_Q9_common.R` (sourced by every script) loads `git_repo/utils` and replaces `palette_sample` with `R2_Q9_palette_sample.R`.
`run_step.sh <script>` runs a script in the `meso_scatac` environment.

| script | original | output in `Plots/` |
|---|---|---|
| `R2_Q9_02_sample_umaps.R` | sample UMAP of the compartment scripts | `{main,TNK,myeloid,tumor,Bcells,stroma}_sample_umap.pdf` (page 2 = old palette) |
| `R2_Q9_03_main_QC.R` | `scatac_main_ArchR.R` QC plots | `main_QC_plots.pdf` |
| `R2_Q9_04_main_CNV_heatmap.R` | `scatac_main_ArchR.R` CNV heatmap | `main_GL_method_z_CNV_heatmap2.pdf` |
| `R2_Q9_05_tumor_plots.R` | `scatac_tumor_SOX9.R` scS violin, SOX6/9 scatter | `tumor_Sarcomatoid_signatures_boxplot.pdf`, `tumor_SOX6_9_scSscore_scatterplot.pdf` |
| `R2_Q9_06_myeloid_cnmf_heatmap.R` | `scatac_myeloid_ArchR.R` | `myeloid_cnmf_scatac_scaled_heatmap2.pdf`, `..._old_palette.pdf` |
| `R2_Q9_09_tumor_DAM_heatmap.R` | `scatac_tumor_SOX9.R` | `tumor_DAM_clusters_Clusters_sample_heatmaps.pdf` |
| `R2_Q9_10_endothelial_plots.R` | `scatac_endothelial_ArchR.R` | `endothelial_adult_vs_fetal_scatterplot_smoothed.pdf` (page 2 = old palette) |
