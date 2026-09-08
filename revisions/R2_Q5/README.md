# R2_Q5 — does the chromatin-based cell type annotation agree with the transcriptomic one?

Reviewer 2, Question 5: *"Does cell type annotation clustered with open chromatin
accessibility overlap with cell annotation using transcriptomic data?"*

---

## The two annotations being compared

| | scATAC | scRNA |
|---|---|---|
| object | `main/scatac_ArchR` | `main/scrna/srt.rds` |
| cells | 49,849 | 58,950 |
| patients | 11 (P1, P3, P4, P5, P8, P10–P14, P23) | 10 (the same, **without P23**) |
| label column | `celltype_lv1` | `celltype_lv1` |
| cell types | 12 | 13 (12 + Glia) |
| how it was made | IterativeLSI on the TileMatrix (25,000 variable features, LSIMethod 2) → `addClusters(resolution = 3)` → 59 clusters → each cluster named from the **gene score** of canonical markers | standard Seurat clustering → clusters named from **expression** of canonical markers |

The chromatin annotation was produced in
`git_repo/main_analysis/scatac_main_ArchR.R` and **never used the scRNA object**.
The two are therefore independent annotations of the same tumours, which is what
makes the comparison meaningful.

**These are not multiome data.** The two assays were run on separate cells from the
same patients; step 1 confirms it — 5 barcodes coincide out of 49,849 × 58,950, which
is what random 16-mer collision gives. There is no ground-truth per-cell pairing, so
the comparison is made three ways, with decreasing reliance on integration.

---

## Three tests, deliberately

| step | test | what it assumes | what it can show |
|---|---|---|---|
| 2–3 | **cell-by-cell label transfer** — CCA anchors from RNA expression onto ATAC gene score; every ATAC cell gets the RNA label it looks most like | that the anchor transfer works | per-cell agreement, per-cell-type recall, which types are confused |
| 4 | **pseudobulk / marker enrichment** — no integration at all; RNA marker sets scored against ATAC chromatin | only that a gene has the same name in both | whether the two annotations name the same biology |
| 5 | **composition** — per-patient cell type fractions, one annotation against the other | nothing beyond the labels | whether the annotations are interchangeable for downstream composition work |

Step 7 then asks what the one substantial disagreement is, and tests a falsifiable
prediction about which tumours should show it. Steps 8–9 put both assays in one CCA
embedding so the agreement can be seen as well as counted.

If the three tests agree, the answer does not depend on any one method. They do.

---

## How to reproduce

```bash
cd .../git_repo_claude/R2_Q5
./submit_02_label_transfer.sh    # LSF, ~1-3 h -- writes atac_label_transfer.csv
./run_all.sh                     # steps 1, 3, 4, 5, 6 once step 2 has finished
```

`run_step.sh` activates the `meso_scatac` conda environment and runs one script.
Steps 1 and 4 read the ArchR project and the Seurat object (a few GB); steps 3, 5
and 6 read only the CSVs in this folder and take seconds.

**Nothing is written back into `main/scatac_ArchR` or `main/scrna`.** Step 2 reads
the project and does the transfer in memory; `saveArchRProject()` is never called.
Step 4 redirects ArchR's output directory into `archr_scratch/` inside this folder
before calling `getGroupSE()`.

---

## Files

| file | contents |
|---|---|
| `atac_cells.csv` | every scATAC cell: patient, LSI cluster, chromatin label, nFrags, TSS, UMAP |
| `rna_cells.csv` | every scRNA cell: patient, level-1 and fine transcriptomic label |
| `annotation_inventory.csv` | cells per patient per modality |
| `atac_label_transfer.csv` | transferred RNA label + confidence per ATAC cell (step 2) |
| `atac_joined.csv` | the two labels side by side, per cell |
| `concordance_summary.csv` | agreement / kappa / ARI, overall and by subset |
| `confusion_matrix.csv` | chromatin label × transferred label, counts |
| `concordance_by_celltype.csv` | recall, precision, F1, commonest disagreement per type |
| `concordance_by_score.csv` | agreement as a function of transfer confidence |
| `concordance_by_cluster.csv` | the 59 LSI clusters: majority label in each modality |
| `pseudobulk_atac.csv`, `pseudobulk_rna.csv` | mean profile per cell type per modality |
| `pseudobulk_correlation.csv` | Pearson of z-scored profiles (Spearman version alongside) |
| `pseudobulk_bestmatch.csv`, `marker_bestmatch.csv` | best RNA match per chromatin cell type |
| `marker_enrichment.csv`, `rna_markers.csv` | RNA marker sets scored in ATAC chromatin |
| `marker_genescore_by_atac_label.csv` | the 35 canonical markers the annotation was made from |
| `composition_by_patient.csv`, `composition_concordance.csv` | per-patient cell type fractions |
| `malignant_disagreement_by_sample.csv` | where each tumour's malignant cells are sent, vs its sarcomatoid score |
| `joint_embedding.csv` | joint CCA UMAP coordinates for all 108,799 cells, with assay, cell type, patient |
| `joint_mixing_by_celltype.csv` | cross-assay neighbour fraction per cell type per assay |
| `joint_label_agreement.csv` | each ATAC cell's chromatin label vs the majority label of its nearest RNA neighbours |
| `joint_centroid_distance.csv` | per cell type, assay-to-assay centroid distance over distance to other cell types |
| `joint_dispersion.csv` | spread of each assay around its own centroid (the imputation shrinkage) |
| `Plots/R2Q5_joint_embedding.pdf` | the joint embedding figure, panels A–F |
| `Plots/R2Q5_annotation_concordance.pdf` | the figure, panels A–H (panels also written separately) |

---

## Results

### Short answer

**Yes.** Across 49,849 scATAC cells the two annotations give the same label to
**81.1 %** of cells (Cohen's kappa **0.76**, ARI **0.63**); restricted to the 10
patients that have matched scRNA it is **88.7 %** (kappa **0.86**, ARI **0.75**).
At the level the chromatin annotation was actually made — the 59 LSI clusters — the
majority label agrees in **51 of 59 clusters, covering 90.0 % of cells**.
Every immune and stromal cell type agrees at 92.6–99.9 % recall. Essentially all of
the disagreement is one cell type, and it has a specific, testable cause.

| subset | cells | agreement | kappa | ARI | compartment |
|---|---|---|---|---|---|
| all ATAC cells | 49,849 | 0.811 | 0.760 | 0.630 | 0.857 |
| patients with matched scRNA | 31,621 | **0.887** | **0.860** | 0.755 | 0.940 |
| P23 (ATAC only, no matched RNA) | 18,228 | 0.678 | 0.559 | 0.471 | 0.713 |
| transfer confidence ≥ 0.5 | 40,915 | 0.875 | 0.845 | 0.762 | 0.913 |

### Per cell type

| chromatin label | n | recall | precision | median confidence | commonest disagreement |
|---|---|---|---|---|---|
| Myeloid | 6,394 | **0.999** | 0.985 | 0.99 | pDCs |
| B_cells | 2,802 | 0.998 | 0.877 | 1.00 | T_cells |
| Endothelial | 1,707 | 0.994 | 0.801 | 1.00 | T_cells |
| pDCs | 165 | 0.988 | 0.872 | 0.98 | B_cells |
| Fibroblasts | 1,953 | 0.988 | 0.275 | 1.00 | SmoothMuscle |
| NK | 1,948 | 0.984 | 0.456 | 0.98 | Fibroblasts |
| Plasma | 141 | 0.979 | 0.780 | 0.75 | B_cells |
| Alveolar | 195 | 0.964 | 0.787 | 0.96 | T_cells |
| SmoothMuscle | 336 | 0.926 | 0.915 | 0.87 | Fibroblasts |
| Mesothelium | 171 | 0.889 | 0.393 | 0.78 | Fibroblasts |
| T_cells | 12,726 | 0.829 | 0.938 | 0.69 | NK |
| **Malignant** | 21,311 | **0.666** | 0.997 | 0.63 | **Fibroblasts (23.5 %)** |

Low **precision** for Fibroblasts (0.275), NK (0.456) and Mesothelium (0.393) is the
same fact seen from the other side: those labels absorb the malignant cells the
transfer cannot place, not a failure of the chromatin annotation of the cells that
actually are fibroblasts (recall 0.988).

### The one substantial disagreement, and what it is

All 8 LSI clusters whose majority labels disagree are **Malignant** clusters, and
23.5 % of chromatin-defined malignant cells are transferred as **Fibroblasts**.

That is not noise. Sarcomatoid mesothelioma is spindle-cell and mesenchymal by
definition, and an RNA reference whose "Fibroblasts" are normal stroma has no other
label to offer such a cell. Step 7 tests the prediction that follows: the per-tumour
rate of Malignant → Fibroblasts should track the tumour's sarcomatoid score, computed
independently from the scATAC data in `R2_Q14`.

It does — **Spearman rho = 0.82, p = 0.007** over the 9 tumours with a score
(rho = 0.75, p = 0.052 over the 7 with ≥ 50 malignant cells):

| tumour | malignant cells | kept Malignant | → Fibroblasts | sarcomatoid score |
|---|---|---|---|---|
| P1 | 192 | 43.8 % | 52.1 % | **0.94** |
| P14 | 131 | 42.0 % | 56.5 % | **0.68** |
| P13 | 18 | 0.0 % | 94.4 % | **0.40** |
| P3 | 21 | 71.4 % | 19.0 % | 0.23 |
| P12 | 832 | 55.0 % | 24.0 % | 0.12 |
| P5 | 1,920 | 98.3 % | 0.8 % | 0.10 |
| P11 | 1,980 | 60.8 % | 13.0 % | 0.09 |
| P4 | 3,082 | 89.6 % | 6.0 % | −0.04 |
| P8 | 472 | 97.0 % | 1.3 % | −0.05 |
| P10 | 976 | 96.7 % | 2.7 % | n/a |
| P23 | 11,687 | 54.1 % | 35.3 % | n/a |

The epithelioid tumours (P5, P8, P4, P10) keep 90–98 % of their malignant cells; the
sarcomatoid ones (P1, P14) lose about half. The disagreeing cells are still counted
as disagreements in every number above — this explains which tumours disagree, it
does not repair the agreement statistic.

### Joint CCA embedding (steps 8–9)

CCA anchors → `TransferData` imputes an expression profile for every ATAC cell from
the RNA reference → merge with all 58,950 real RNA cells → centre → PCA → UMAP.
108,799 cells in one space; cell type labels were assigned independently in each
assay and never matched.

**The space is organised by cell type, not by assay — 12 of 12.** For every type the
distance between the two assays' centroids is a small fraction of the distance to
other cell types:

| cell type | ratio | cell type | ratio |
|---|---|---|---|
| Endothelial | 0.046 | Myeloid | 0.114 |
| B_cells | 0.059 | Plasma | 0.115 |
| pDCs | 0.067 | T_cells | 0.130 |
| Fibroblasts | 0.078 | Alveolar | 0.140 |
| NK | 0.078 | SmoothMuscle | 0.171 |
| Mesothelium | 0.225 | **Malignant** | **0.402** |

Asking each ATAC cell for the majority label among its nearest *RNA* neighbours
reproduces the cell-by-cell result: **86.9 % agreement** (kappa 0.83, ARI 0.72), with
Malignant at 75.9 % and losing 14.6 % to Fibroblasts — the same sarcomatoid signal.

**But the two assays do not interleave, and this must not be over-read.** ATAC cells
draw 0.000–0.002 of their 30 nearest neighbours from the other assay, against 0.542
expected; overall mixing is 0.119 where a mixed space gives 0.497. Restricting to a
50/50 balanced subsample (99,602 cells) gives 0.117, so cohort composition is not the
cause.

The cause is the imputation. `TransferData` returns an anchor-weighted average, so
imputed profiles are shrunk toward the local reference mean: the ATAC cloud is about
**half as wide** as the RNA one (median spread ratio **0.484**, from 0.24 for pDCs to
0.73 for Malignant). It is not in the *wrong place* — the median ATAC cell is only
1.27× further from its nearest RNA cell than an RNA cell is from its nearest RNA
neighbour — it is simply denser, so k-nearest-neighbour mixing picks up only ATAC
cells. Panel D is a density artefact; panel E is the result, and panel F is the
explanation.

**How to use this.** The joint embedding is evidence that the two annotations name
the same populations (panels B, C, E). It is **not** evidence of integration or
batch-correction quality, and overlap in panel A is expected by construction even for
a bad transfer. Note also that the joint-space agreement (86.9 %) exceeds the
transfer agreement (81.1 %) because imputation flatters itself; **quote 81.1 % / 88.7 %,
not 86.9 %**.

If genuine cross-assay interleaving is wanted, that is a different procedure —
`FindIntegrationAnchors` + `IntegrateData` treating gene score and expression as two
batches, which corrects both directions rather than projecting ATAC onto RNA. It
answers "can the assays be batch-corrected together", not "do the annotations agree",
and is not implemented here.

### The disagreement is concentrated at low transfer confidence

| predictedScore | cells | % of cells | agreement |
|---|---|---|---|
| < 0.2 | 29 | 0.1 % | 0.24 |
| 0.2–0.4 | 4,331 | 8.7 % | 0.43 |
| 0.4–0.5 | 4,574 | 9.2 % | 0.60 |
| 0.5–0.6 | 5,548 | 11.1 % | 0.68 |
| 0.6–0.8 | 10,930 | 21.9 % | 0.82 |
| ≥ 0.8 | 24,437 | **49.0 %** | **0.94** |

Half the cells sit at confidence ≥ 0.8 and agree 94 % of the time. So a large share
of the disagreement is the transfer being unsure, not the two annotations making
confident opposite claims.

### Confirmation without any integration (step 4)

Scoring the **top 100 RNA marker genes of each cell type against ATAC chromatin**
puts the matching cell type first for **11 of 12** types, and **10 of 12** exceed the
95th percentile of a 1,000-permutation null over random gene sets of the same size:

```
Alveolar 2.13 | SmoothMuscle 1.64 | pDCs 1.57 | Mesothelium 1.49 | Endothelial 1.48
Fibroblasts 1.27 | NK 1.19 | Plasma 1.10 | B_cells 0.92 | Malignant 0.38
Myeloid 0.28 | T_cells -0.04            (null 95th pct = 0.295)
```

The single miss is **Malignant**, whose best match is **Mesothelium** (0.53 vs 0.38)
— mesothelioma's tissue of origin, so the two chromatin profiles genuinely overlap.
Profile correlation (Pearson on z-scored pseudobulk, 2,000 genes) agrees at
**10 of 12**, missing Mesothelium → Malignant and Plasma → B_cells (141 and 179 cells).

### Composition (step 5)

Over 130 patient × cell-type pairs from the 10 shared patients, the two annotations
give cell type fractions correlated at **r = 0.946** (Spearman 0.881). Per patient
r runs 0.55 (P1) to 1.00 (P8). Two systematic offsets, neither surviving FDR:

| cell type | mean fraction, ATAC | mean fraction, RNA | difference | FDR |
|---|---|---|---|---|
| Myeloid | 0.109 | 0.170 | **−0.061** | 0.24 |
| Malignant | 0.346 | 0.297 | **+0.049** | 0.24 |

Both are in the direction expected from nucleus extraction rather than annotation:
scATAC under-recovers myeloid cells and over-represents tumour. Note this test cannot
separate annotation disagreement from assay bias, and is reported as the weakest of
the three.

---

## Caveats

- **P23 has no matched scRNA.** Its 18,228 cells (37 % of the ATAC cohort) can only
  borrow a label from another patient's tumour, and they agree at 67.8 % against
  88.7 % for the matched patients. P23 alone accounts for most of the gap between the
  cohort-wide and matched-patient numbers. Quote the matched-patient figure for the
  annotation, the cohort figure for the dataset as delivered.
- **Malignant cells are patient-specific.** A transfer that must map one tumour's
  malignant cells onto a reference built from other tumours is the hardest case in
  this dataset, and it is where all the disagreement is.
- **Glia has no chromatin counterpart.** The RNA object has 55 Glia cells; the ATAC
  annotation has no Glia cluster, and 29 ATAC cells are transferred to it. Too few
  cells in either modality to say whether this is a real miss.
- **Rare types are near the floor.** Plasma (141 ATAC / 179 RNA), Alveolar (195/147)
  and Mesothelium (171/318) have noisy pseudobulk profiles; their step 4 results
  should be read with that in mind.
- **Step 5 confounds annotation with assay.** Composition differences can come from
  the labels or from dissociation and nucleus-extraction bias, and these data cannot
  separate them.

---

## Software notes

- `ArchR 1.0.2`, `Seurat 5.2.1`, `SeuratObject 5.0.2`, `Matrix 1.6.5`, R 4.3.3,
  conda environment `meso_scatac`.
- **ArchR's `addGeneIntegrationMatrix()` does not run against Seurat 5** in this
  environment: it mixes `CreateSeuratObject()` with the v3-only
  `CreateAssayObject()`, and the block loop dies with
  `slot(object, "features")[[layer]] <- features : more elements supplied than there
  are to replace`. Setting `options(Seurat.object.assay.version = "v3")` is not
  sufficient. Step 2 therefore implements the same procedure directly in Seurat —
  gene score as the query feature space, CCA anchors, ArchR's IterativeLSI as the
  transfer weight reduction — which also puts every parameter in the open.
- **Do not rebuild ArchR's gene score `dgCMatrix`** with `Matrix::sparseMatrix()` to
  attach gene symbols to its rows; the result no longer dispatches `[`
  (`object of type 'S4' is not subsettable`). Index the original by position and name
  the small per-block matrix instead.
- **The joint embedding is imputation-based** and therefore cannot be used to argue
  integration quality; see the joint-embedding section. Its agreement number is
  optimistic relative to step 3 and should not be the one quoted.
- **Step 4 uses Pearson, not Spearman.** Ranking 2,000 mostly-flat z-profiles lets
  noise dominate: Spearman puts the diagonal first for only 6 of 12 cell types and
  makes Plasma the spurious best match for three unrelated types, while Pearson gets
  10 of 12. The Spearman matrix is written to
  `pseudobulk_correlation_spearman.csv` so the choice can be checked.
