# R2_Q5 — joint embedding of the scATAC and scRNA cells

Reviewer 2, Q5: does the cell type annotation from chromatin accessibility overlap with the annotation from expression?
The two annotations were made independently (`celltype_lv1` of `main/scatac_ArchR` and of `main/scrna/srt.rds`).
Method text: `METHODS_joint_embedding.md`.

`./run_all.sh` runs the three steps (step 8 needs ~1-2 h and 40 GB: use `submit_08_joint_embedding.sh` on LSF).

| step | script | what it does | output |
|---|---|---|---|
| — | `00_common.R` | paths, label vocabulary, colours | — |
| 1 | `01_annotations.R` | per-cell labels of both assays | `atac_cells.csv`, `rna_cells.csv`, `annotation_inventory.csv` |
| 8 | `08_joint_embedding.R` | CCA anchors, RNA expression imputed for every ATAC cell, shared PCA/UMAP; mixing, centroid and dispersion statistics | `joint_embedding.csv`, `joint_*.csv` |
| 9 | `09_figure_joint.R` | joint UMAP by assay (A) and by cell type (B); for each scATAC cell type, the label of the nearest scRNA cell in the joint UMAP | `Plots/R2Q5_joint_panel_A.pdf`, `Plots/R2Q5_joint_panel_B.pdf`, `Plots/scatac_celltype_comp_in_scrna.pdf` |

`run_step.sh` activates the `meso_scatac` environment and runs one script.
