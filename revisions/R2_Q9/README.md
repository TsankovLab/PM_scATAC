# R2_Q9 — improved sample palette and regenerated per-sample figures

The reviewer asked for a more immediate sample colour scheme. The old `palette_sample` (`git_repo/utils/palettes.R`) used ten muted
brown/grey/blue tones for the tumours and the same red for every normal sample.

## New palette (`R2_Q9_palette_sample.R`, drop-in replacement)
- **Tumours:** 11 clearly different hues, assigned along the sarcomatoid score (median scS of malignant cells, scATAC):
  cool colours (green, teal, blue) = epithelioid-like, violet/purple = intermediate, warm colours (pink, amber, red) = higher scS.
- **Normal samples:** black/greys, one shade each; `normal1` (the one plotted with the tumours) is near-black.
- **P11_HOX:** light tint of the P11 colour.
- Smallest pairwise colour distance among the 11 tumours (OKLab x 100): 18.3 for normal vision and 10.1 under simulated red-green
  colour blindness, against 5.0 and 2.1 for the old palette (`R2_Q9_palette_distances.csv`, `Plots/R2_Q9_palette_old_vs_new.pdf`).
- The mid greys of normal2/normal3/HU37/HU62 can come close to P10/P12 under colour-blind simulation; these samples are not plotted
  together with the tumours in any scATAC figure.

## Scripts (plotting code copied from the originals, only the palette changes; projects are read-only)
| script | original | output in `Plots/` |
|---|---|---|
| `R2_Q9_00_palette_design.R` | — | chooses the 11 tumour colours |
| `R2_Q9_01_palette_validation.R` | — | `R2_Q9_palette_old_vs_new.pdf` |
| `R2_Q9_02_sample_umaps.R` | sample UMAP of the 7 compartment scripts | `<compartment>_sample_umap.pdf` |
| `R2_Q9_03_main_QC.R` | `scatac_main_ArchR.R` QC plots | `main_QC_plots.pdf` |
| `R2_Q9_04_main_CNV_heatmap.R` | `scatac_main_ArchR.R` CNV heatmap | `main_GL_method_z_CNV_heatmap2.pdf` |
| `R2_Q9_05_tumor_plots.R` | `scatac_tumor_SOX9.R` scS violin, SOX6/9 scatter | `tumor_Sarcomatoid_signatures_boxplot.pdf`, `tumor_SOX6_9_scSscore_scatterplot.pdf` |
| `R2_Q9_06_myeloid_cnmf_heatmap.R` | `scatac_myeloid_ArchR.R` | `myeloid_cnmf_scatac_scaled_heatmap2.pdf` |
| `R2_Q9_07_R2_Q5_joint_figure.R` | `revisions/R2_Q5/09_figure_joint.R` | `R2Q5_joint_*.pdf` |
| `R2_Q9_08_Bcells_TF_module_hist.R` | `scatac_Bcells_ArchR.R` | `Bcells_TF_modules_dist_sample_hist.pdf` |
| `R2_Q9_09_tumor_DAM_heatmap.R` | `scatac_tumor_SOX9.R` | `tumor_DAM_clusters_Clusters_sample_heatmaps.pdf` |
| `R2_Q9_10_endothelial_plots.R` | `scatac_endothelial_ArchR.R` | `endothelial_fetal_TFs_heatmap.pdf`, `endothelial_adult_vs_fetal_scatterplot_smoothed.pdf` |

Old-palette versions are kept for comparison (second page of a pdf, or `*_old_palette.pdf`).

Not regenerated: `scatac_stroma_hub_analysis.R` `scrna_region_boxplots.pdf` (the block uses objects that the script does not define:
`sample_sarc_order_levels` and a scRNA object with `sampleID3`).
