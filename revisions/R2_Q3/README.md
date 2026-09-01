# R2_Q3 — subclonal structure in the scATAC tumours (epiAneufinder)

Reviewer 2, Question 3. Three things, in order:

1. **epiAneufinder subclones** in every scATAC tumour, from the malignant-cell fragments.
2. **TF activity between the clones**, across samples (chromVAR).
3. **Spatial validation of the P4 chr8q clone** in Visium.

Everything in this folder belongs to that chain. The earlier, hand-rolled arm-level CNV
pipeline (Granja–Lareau windows, cohort-adjusted multimodal arm scores, inferCNV
cross-modal trees) has been removed — epiAneufinder replaces it and reaches the same
conclusion about P4 with far less machinery.

---

## How to reproduce

```bash
cd .../git_repo_claude/R2_Q3
./run_all.sh                      # steps 3-9, ~7 min on one interactive core
```

Ask for a node with a few GB: step 3 holds an 11 331 x 11 331 distance matrix for P23
(~0.5 GB) and a login-node run of it was killed. `bsub -Is -P acc_Tsankov_Normal_Lung
-q premium -n 4 -R "rusage[mem=8000]" -W 2:00 /bin/bash` is plenty.

`run_all.sh` deliberately starts at **step 3**, because `out_5Mb/` already holds the
epiAneufinder output. To redo the two expensive steps from scratch:

```bash
./run_step.sh 01_extract_fragments.R          # ~10 min, writes ~11 GB into frags/
for S in P1 P4 P5 P8 P10 P11 P12 P14 P23; do
  bsub -P acc_Tsankov_Normal_Lung -q premium -n 8 -R "rusage[mem=8000] span[hosts=1]" \
       -W 6:00 -J epi_$S -o logs/epi_$S.out -e logs/epi_$S.err \
       "$PWD/run_step.sh 02_run_epianeufinder.R $S"
done
```

`run_step.sh` activates the `meso_scatac` conda environment and runs one script; use it
for any single step.

## Scripts

| | script | what it does | reads | writes |
|---|---|---|---|---|
| — | `00_common.R` | paths, epiAneufinder parameters, the clone-split function, palette. Sourced by every step, so the clone definition cannot drift between figures and statistics. | — | — |
| 1 | `01_extract_fragments.R` | per-cell fragments for each tumour out of the ArchR project | `tumor_compartment/scatac_ArchR` | `frags/<S>_fragments_archr.tsv` |
| 2 | `02_run_epianeufinder.R <S>` | epiAneufinder at 5 Mb, one sample | `frags/` | `out_5Mb/<S>/epiAneufinder_results/` |
| 3 | `03_clone_calls.R` | **defines the clones**; clone CNV profiles, driver arms | `out_5Mb/` | `epi_clone_labels.csv`, `epi_clone_profiles.rds`, `epi_leaves.csv`, `epi_subclone_summary.csv` |
| 4 | `04_figure_subclones_overview.R` | landscape of all 18 subclones + what separates each pair | step 3 | `Plots/epiAneufinder_subclones_overview.pdf` |
| 5 | `05_figure_circular_tree.R` | one circular tree over all 18 subclones | step 3 | `Plots/epiAneufinder_circular_tree.pdf` |
| 6 | `06_P4_chr8q_validation.R` | is P4's split really chr8q? chr8p as control | step 3 | `P4_chr8q_epiAneufinder.csv`, `Plots/P4_chr8q_epiAneufinder.pdf` |
| 7 | `07_chromvar_tf_variability.R` | TF activity between clones, per sample and across | step 3, ArchR MotifMatrix | `epi_chromvar_*.csv`, `chromvar_z_cache.rds` |
| 8 | `08_figure_tf_variability.R` | TF variability figure | step 7 | `Plots/epi_TF_variability.pdf` |
| 9 | `09_P4_visium_chr8q.R` | the P4 chr8q clone in space | Visium object | `P4_visium_chr8q_spatial.csv`, `Plots/R2_Q3_P4_visium_chr8q.pdf` |

Steps 4–9 depend only on step 3, so any one of them can be re-run on its own.

## Method, in brief

**Cells.** The ArchR project `tumor_compartment/scatac_ArchR` contains the malignant
cells only, so "every cell of sample *S* in that project" *is* the malignant population
— there is no separate barcode list that could fall out of sync. Autosomes only;
fragments are written as `chr / start / end / barcode`.

**CNV calling.** epiAneufinder 1.1.5, **5 Mb windows**, ENCODE hg38 v2 blacklist,
chrX/chrY/chrM excluded, GC correction on, `minFrags = 5000`, `k = 4`. Output is a
per-cell, per-window call: 0 = loss, 1 = normal, 2 = gain, on a 298-window genome grid
shared by all nine samples.

*Why 5 Mb.* At 10 Mb only ~100 windows survive the blacklist genome-wide (4–5 per
chromosome) and the caller errors out at gain/loss assignment. At 100 kb the data are
too sparse — 99.4 % of calls come back "normal". *Why `minFrags = 5000`.* The package
default of 20 000 would keep 154 of P4's 3 082 malignant cells.

**Clones.** epiAneufinder does not estimate a number of subclones: `split_subclones`
is literally `cutree(hclust(dist, "ward.D"), k = tree_depth)`. Each sample's own
dendrogram is therefore cut at **k = 2** — the depth that isolates P4's chr8q clone —
and the same depth is applied everywhere so the samples stay comparable. A clone must
hold ≥ 20 cells and ≥ 5 % of the sample. A clone's CNV profile is, per window,
(fraction of its cells called gain) − (fraction called loss).

**Silhouettes are ≈ 0 at this resolution**, so k is a fixed analysis choice, not a
data-driven estimate. What makes P4 credible is not the silhouette: it is that the
clone is arm-specific (chr8q, not chr8p), that it survives every k from 2 to 8, and
that it reappears in an independent modality (Visium).

**TF activity.** chromVAR z from the ArchR `MotifMatrix` (869 cisBP motifs), cached to
`chromvar_z_cache.rds`. Per sample: difference of mean z between clones, η², Wilcoxon +
BH, and the ARI of an *independent* chromVAR-only k = 2 clustering against the CNV
clones. Across samples: median |difference| with IQR — median, not mean, because samples
differ several-fold in how far apart their clones sit. The FDR column ranks motifs
*within* a sample only; it scales with cell count and must not be compared across
samples.

AP-1/bZIP + SMARCC1 are flagged **technical**: that axis tracks per-cell fragment count
in this dataset and floors in every comparison. Do not read it as clonal biology.

## Results

**Subclones.** All nine tumours split into two clones (`epi_subclone_summary.csv`).
P4 is the only split driven by a clean, focal, single-arm event.

| sample | clone 1 | clone 2 | top arm | Δ (clone2 − clone1) |
|---|---|---|---|---|
| P1 | 81 | 109 | 12q | 0.462 |
| **P4** | 2277 | **555** | **8q** | **0.759** |
| P5 | 1493 | 343 | 1q | 0.954 |
| P8 | 107 | 152 | 12p | 0.342 |
| P10 | 538 | 351 | 11q | 0.839 |
| P11 | 1154 | 787 | 15q | −0.318 |
| P12 | 374 | 312 | 10p | −0.341 |
| P14 | 44 | 68 | 12p | 0.316 |
| P23 | 10036 | 1295 | 1q | 0.898 |

⚠️ **P5 and P23 peak at the same chr1 bins** (Δ 0.95 and 0.90). A cross-tumour
signature at identical coordinates is what a gene-density / mappability confound looks
like; check those bins before reporting either as a chr1q subclone. P4's chr8q has no
such twin.

**P4 is a genuine chr8q clone** (`P4_chr8q_epiAneufinder.csv`). 555 of 2 832 cells
(19.6 %):

- chr8q gain fraction **0.688 vs 0.129** (Wilcoxon p = 2.6 × 10⁻²⁶²)
- chr8p, the control arm, shifts the same way but **~30× less** (0.048 vs 0.029)
- against a *targeted* split (k-means on the chr8q fraction alone): ARI 0.483,
  precision 0.939, recall 0.583 — i.e. the genome-wide clustering finds the chr8q
  clone **without being pointed at chr8**, and what it finds is almost pure

**Spatial validation** (`P4_visium_chr8q_spatial.csv`). chr8q-high malignant spots form
contiguous territories on both slides, not salt-and-pepper:

| slide | malignant spots | chr8q-amp | Moran's I | p |
|---|---|---|---|---|
| C1 | 1343 | 72 % | **0.760** | < 0.0033 |
| D1 | 1235 | 31 % | **0.664** | < 0.0033 |

(299 permutations, so p = 0 means below the resolution of the null.)

**TF activity does not track the clones.** ARI between an independent chromVAR
clustering and the CNV clones is ≈ 0 or negative in 8 of 9 tumours (only P11 = 0.100;
P4 = −0.015), and mean η² is 0.0002–0.03. The clones are CNV states, not distinct
regulatory states.

What *is* consistent across tumours is a set of motifs that move between clones
everywhere: **TEAD3** (median |Δz| 0.68, IQR 0.50–0.76 — the top non-technical motif in
P8 and P10 individually), then NFYA 0.48, CEBPZ 0.43, PITX2 0.42, NFYB/NFIC/NFYC ≈ 0.41,
RELA 0.37, TEAD4/TEAD1 0.34. Immune/interferon motifs are enriched among the top decile
but not significantly (OR 2.30, p = 0.056).

## Files kept

- `out_5Mb/<S>/epiAneufinder_results/` — the CNV calls (the analysis input; 44 MB)
- `frags/` — per-sample fragment files (11 GB; regenerable by step 1)
- `chromvar_z_cache.rds` — 869 motifs × 21 723 cells (128 MB; regenerable by step 7)
- `epi_*.csv`, `P4_*.csv` — all result tables
- `Plots/` — the five figures
- `logs/` — per-step run logs

## Environment

conda `meso_scatac` · R 4.3.3 · epiAneufinder 1.1.5 · ArchR 1.0.2 · chromVAR 1.24.0 ·
ComplexHeatmap 2.18.0 · ggtree 3.10.1 · Seurat 5.2.1 · data.table 1.15.4
