# R1_Q4 — myeloid scATAC cells projected onto the PBMC / bone marrow reference

Run in this order (each R script has a matching LSF submit script):

1. `add_myeloid_compartment_col.R` (`submit_add_myeloid_compartment.sh`) — adds `celltype_main = "Myeloid"` to the myeloid
   ArchR project so that all cells can be projected in one call (saves the project).
2. `project_myeloid_on_pbmc.sh` — projects all myeloid cells onto the IterativeLSI space of the Greenleaf PBMC/BM reference
   (`Xmen/greenleaf_PBMC_BM`, label column `BioClassification`) with the projection wrapper
   `/sc/arion/scratch/leew17/MesoRNA/ProjectionCode/run_projection_default.sh`; writes `output/pbmc_projection/myeloid_all_lsi_projection.rds`.
3. `compute_per_label_scores.R` (`submit_compute_per_label_scores.sh`) — per cell, the fraction of its 25 nearest reference
   cells (LSI space) carrying each reference label; writes `output/pbmc_projection/per_label_scores.csv`.
4. Figures (read `per_label_scores.csv`):
   - `plot_all_label_scores_boxplots.R` → `output/pbmc_projection/all_label_scores_boxplots.pdf`
   - `plot_mono_score_boxplots.R` → `output/pbmc_projection/boxplot_mono_scores.pdf`
   - `plot_mod2_vs_cd14mono2.R` (also reads `mod_2` from the myeloid ArchR project) → `output/pbmc_projection/scatter_mod2_cd14mono2.pdf`
