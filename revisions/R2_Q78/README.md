# R2_Q78 — peak-to-gene links of cHub2711 (HAND2) and cHub3727 / cHub3728 (MEF2C)

Run in this order (`run_step.sh <script>` runs a script in the `meso_scatac` environment):

1. `R2_Q78_peak2gene.R` — ArchR `addPeak2GeneLinks` on `main/scatac_ArchR` (GeneScoreMatrix, IterativeLSI) with maxDist
   500 kb instead of the manuscript's 250 kb, because the MEF2C hubs lie 363-387 kb from the MEF2C TSS. The project's output
   directory is redirected to scratch (not modified); the links are cached to
   `/sc/arion/scratch/giottb01/R2_Q78_peak2gene/peak2gene_links_500kb.rds`. Tables of every hub peak vs the expected gene:
   `R2_Q78_peak2gene_hub_peaks.csv`, `R2_Q78_peak2gene_all_genes_for_hub_peaks.csv`, `R2_Q78_peak2gene_summary.csv`.
2. `R2_Q78_p2g_links_tracks.R` — browser tracks of the two loci side by side: coverage per cell type, all cHubs in the window
   (labelled), peak-to-gene links within 400 kb (correlation >= 0.2, FDR < 1e-4 recomputed over the 400 kb set), genes.
   Output: `Plots/R2_Q78_cHubs_p2g_links.pdf`, `R2_Q78_p2g400_hub_peak_links.csv`, `R2_Q78_p2g400_links_to_HAND2_MEF2C.csv`.

Hub set: `main/scatac_ArchR/hubs_obj_cor_0.3_md_12500_dgs_0_min_peaks_5`. HUB3722 (requested) is 14.5 Mb from MEF2C; the
second MEF2C hub next to HUB3728 is HUB3727, which is the one plotted.
