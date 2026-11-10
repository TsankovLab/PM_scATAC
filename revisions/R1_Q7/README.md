# R1_Q7 — benchmarking of the cHub definition against super-enhancer sets

Run in this order:

1. `myeloid_SE_hub_enrichment.R` (`submit_myeloid_SE_hub_enrichment.sh`) — hubs from co-accessibility (r 0.2-0.6) and from
   peak stitching (consensus peaks, Stitch-AP; per-cell-type peaks, Stitch CT) at 12.5 kb, for minimum peaks 3/4/5 and TSS
   filter 0/2 kb; differential hubs per cell type; enrichment of each cell type's hubs in every SE set (hypergeometric, BH).
   Output `myeloid_SE_hub_enrichment.csv`.
2. `max_dist_sweep_SE.R` (`submit_max_dist_sweep_SE.sh`) — the same at maximum distance 5, 10, 15 and 20 kb.
   Output `max_dist_sweep_SE.csv`.
3. `max_dist_sweep_summary.R` — SE recovery (matched cell type | SE pairs) and specificity (minus non-matching SE sets)
   per method, setting and distance. Output `plot_max_dist_sweep_table.csv`.
4. `plot_benchmark_summary.R` — A: SE recovery at 12.5 kb for the matched pairs; B: recovery and specificity over all
   distances and settings. Output `plot_benchmark_summary.pdf/.png`.
