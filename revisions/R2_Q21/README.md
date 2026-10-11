# R2_Q21 — CD8 exhausted DAP: gene annotation and figures at FDR < 0.05

The peaks more accessible in CD8 exhausted than in CD8 T cells come from `git_repo/tnk_analysis/scatac_NKTcells_ArchR.R`
(`DAP_CD8_CD8_ext_pairwise.rds`: PeakMatrix, CD8_exhausted vs CD8, Wilcoxon). The original cut-off was FDR < 0.01 and
Log2FC > 0 (418 peaks); the labelled heatmap of `tnk_exhaustion_annotation` uses FDR < 0.05 (884 peaks). Everything
here uses FDR < 0.05 so that table and figures describe the same peak set.

- `R2_Q21_DAP_CD8_Ex_gene_annotation.R` — table of the 884 peaks with a gene per peak, by the rule of the labelled heatmap:
  `gene` = heatmap label if the peak is labelled there, else the Peak2Gene-linked gene with the highest correlation
  (>= 0.35; links stored in the T/NK project by `addPeak2GeneLinks`, 250 kb, filtered as `getPeak2GeneLinks(corCutOff = 0.2)`),
  else the ArchR nearest gene. `gene_source` says which.
  Output: `DAP_CD8_Ex_FDR0.05_gene_annotation.csv`.
- `R2_Q21_Ex_DAP_figures_FDR0.05.R` — the figures of the original script that depend on the cut-off, regenerated at
  FDR < 0.05 with unchanged code: MA plot, TF-module enrichment, TF-in-module enrichment, coverage density, coverage heatmap
  (`Plots/*_FDR0.05.pdf`, plus `TFmodule_enrichments_FDR0.05.csv` and `TF_in_module_enrichments_FDR0.05.csv` with both cut-offs).
  The GSEA plot ranks all peaks and does not depend on the cut-off.
- The labelled heatmaps already use FDR < 0.05: `tnk_exhaustion_annotation/Plots/CD8_ext_peakset_enrichments_annotated_FDR0.05.pdf`
  and `CD8_ext_peakset_focused_FGFBP2_CD8__KLRC1_CD8ex_FDR0.05.pdf`.
- `run_step.sh` — runs a script of this folder in the `meso_scatac` environment.
