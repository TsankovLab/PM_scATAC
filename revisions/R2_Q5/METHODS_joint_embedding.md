# Joint embedding of the scATAC and scRNA cells (R2_Q5, steps 8–9)

Two parts: a paragraph for the paper, then the step-by-step with the reasoning and
the numbers behind each choice.

---

## Part 1 — for the paper

> **Joint scATAC/scRNA embedding.** scATAC and scRNA profiles were generated from
> separate cells of the same tumours (10 shared patients; the assays are not
> multiome and share no cell barcodes), so the two modalities were co-embedded by
> canonical correlation analysis. An scRNA reference was built by subsampling up to
> 3,000 cells per level-1 cell type (21,165 cells), log-normalising, and taking the
> 2,000 most variable genes, of which 1,694 were present in the ArchR gene score
> matrix and defined the shared feature space. scATAC cells were processed in blocks
> of 10,000: gene scores were log-normalised and scaled over the same features, and
> transfer anchors were identified against the reference by CCA on 30 dimensions
> (Seurat `FindTransferAnchors`, 14,528–16,623 anchors per block). Reference
> expression was then transferred onto every scATAC cell with `TransferData`,
> weighting anchors by the cell's ArchR IterativeLSI coordinates (30 dimensions,
> correlation-to-depth cutoff 0.75), giving each scATAC cell an imputed expression
> profile over the shared genes. Imputed scATAC profiles and the log-normalised
> profiles of all 58,950 scRNA cells were concatenated (108,799 cells × 1,694
> genes), gene-centred without variance scaling, reduced to 30 principal components
> (`irlba`) and embedded by UMAP (`uwot`, 30 neighbours, min_dist 0.3, cosine
> metric). Cell type labels were assigned independently within each assay and were
> never used by the embedding. Agreement was quantified in the shared principal
> component space as (i) the majority scRNA label among each scATAC cell's 30
> nearest scRNA neighbours, scored against its chromatin label, and (ii) the
> distance between the two assays' centroids for each cell type, expressed as a
> fraction of that cell type's median distance to the other cell types. Because
> imputation is an anchor-weighted average and therefore contracts the scATAC cells
> toward the local reference mean, the spread of each assay around its own centroid
> was also measured, and nearest-neighbour modality mixing is reported only
> alongside that contraction rather than as an integration-quality statistic.

~250 words. If space is tight, the last sentence is the one to cut — but then the
mixing statistic must be cut with it, not reported alone.

---

## Part 2 — the steps, and why each one

### Step 0 — why a joint embedding at all, and why CCA

The two assays were run on **separate cells**. Five barcodes coincide between 49,849
scATAC and 58,950 scRNA cells, which is what random 16-mer collision gives, so there
is no per-cell pairing to exploit and no ground truth. Any statement that the two
annotations describe the same populations therefore has to be built from shared
*gene-level* signal.

CCA is the right tool because it does not require the two matrices to be on a common
scale. It finds the linear combinations of genes that are maximally correlated
*between* the two datasets, which is what is needed when one matrix is log-normalised
UMI counts and the other is an ArchR gene score — a fundamentally different quantity
(a distance-weighted sum of Tn5 insertions around a gene body). PCA on the two
stacked directly would be dominated by the difference between the assays.

### Step 1 — the shared feature space

**What:** subsample up to 3,000 scRNA cells per level-1 cell type (21,165 cells),
`NormalizeData`, `FindVariableFeatures(nfeatures = 2000)`, intersect with the gene
score matrix → **1,694 genes**. `ScaleData` on those genes.

**Why subsample, and why stratified:** anchor finding scales badly with reference
size, and the cohort is dominated by three types (Malignant 21,882; T_cells 15,232;
Myeloid 8,594). An unstratified subsample would leave Glia (55 cells), Alveolar (147)
and Plasma (179) below the point where they can anchor anything. Capping the large
types and keeping every cell of the small ones costs nothing and preserves the rare
types.

**Necessary?** Yes for the rare types; the cap on the large ones is for runtime.

### Step 2 — anchors, in blocks

**What:** for each block of 10,000 scATAC cells — build a Seurat object from the gene
score matrix, `NormalizeData`, `ScaleData` on the 1,694 shared genes,
`FindTransferAnchors(reduction = "cca", dims = 1:30)`.

**Why blocks:** peak memory. Five blocks gave 14,528 / 16,623 / 15,827 / 15,527 /
15,287 anchors — consistent across blocks, so the blocking is not creating
block-specific structure.

**Why gene score at all:** it is the only representation of an ATAC cell that lives
in gene space, which is the only space the two assays share. Its known weakness —
it inherits copy number, so a gene on an amplified arm scores higher for reasons
unrelated to regulation — is a real limitation and is why the embedding is used for
cell type identity and not for quantitative comparison.

### Step 3 — transfer, weighted by chromatin structure

**What:** `TransferData(refdata = <reference log-normalised expression on the 1,694
genes>, weight.reduction = <ArchR IterativeLSI, 30 dims, corCutOff 0.75>)`. Each
scATAC cell receives an imputed expression profile.

**Why weight by LSI and not by the CCA space:** the anchors say which reference cells
a query cell resembles; the weights say how much each anchor should count for a given
query cell, and those weights should come from a space that describes the query
cells' *own* structure. LSI is the scATAC manifold ArchR clustered on. Weighting by
the CCA space instead would let the transfer reinforce its own alignment.

**`corCutOff = 0.75`** drops LSI components correlated with sequencing depth, so
depth does not drive the anchor weighting.

### Step 4 — concatenate, centre, PCA, UMAP

**What:** `cbind` the imputed scATAC profiles with the log-normalised profiles of
**all** 58,950 scRNA cells → 1,694 × 108,799. Subtract each gene's mean.
`irlba::prcomp_irlba(n = 30, center = FALSE, scale. = FALSE)`. `uwot::umap`
(30 neighbours, min_dist 0.3, cosine).

**Why centre but not scale:** imputed values are already on the reference's
log-normalised scale, so the two halves are commensurate. Unit-variance scaling would
divide by the standard deviation of each gene — and imputation leaves some genes
nearly constant across the scATAC cells, so scaling would multiply exactly those
genes' residual noise up to the same weight as real signal.

**Why all the RNA cells, not the subsample:** the reference subsample exists to make
anchor finding tractable. The embedding has no such constraint, and using every cell
keeps the RNA manifold at full density so the comparison of local structure is fair.

### Step 5 — what is measured, and what is deliberately not claimed

Three statistics in the 30-dimensional PC space, not the UMAP (UMAP distances are not
metric and should never be measured on):

**(a) Joint-space label agreement.** For each scATAC cell, the majority cell type
among its 30 nearest **scRNA** neighbours, scored against its chromatin label.
Result: **86.9 %** (Cohen's kappa 0.83, ARI 0.72).

**(b) Centroid geometry.** Per cell type, the distance between the scATAC and scRNA
centroids divided by that cell type's median distance to the other cell types. Below
1 means the space is organised by cell type rather than by assay. Result: **12 of 12
cell types**, ratios 0.046 (Endothelial) to 0.402 (Malignant).

**(c) Dispersion, and the mixing statistic it explains.** scATAC cells draw
0.000–0.002 of their 30 nearest neighbours from the other assay against 0.542
expected from cohort composition; overall mixing is 0.119 where a mixed space gives
0.497, and a 50/50 balanced subsample (99,602 cells) gives 0.117, so composition is
not the cause.

The cause is the imputation. `TransferData` returns an anchor-weighted average, so
imputed profiles are contracted toward the local reference mean: the scATAC cloud is
about **half as wide** as the scRNA one (median spread ratio **0.484**, from 0.24 for
pDCs to 0.73 for Malignant). It is not in the wrong *place* — the median scATAC cell
is only **1.27×** further from its nearest scRNA cell than an scRNA cell is from its
own nearest neighbour — it is simply denser, so k-nearest-neighbour mixing returns
only scATAC cells.

**What must not be claimed.** Imputation pulls scATAC cells toward the reference by
construction, so visual overlap of the two assays in the UMAP is expected **even for a
poor transfer** and is not evidence on its own. This embedding supports the statement
that the two annotations name the same populations (b, and the side-by-side panels);
it does **not** support any statement about batch-correction or integration quality.
The joint-space agreement (86.9 %) also exceeds the direct label-transfer agreement
(81.1 % cohort-wide, 88.7 % in the patients with matched scRNA) because imputation
flatters itself — **the transfer numbers are the ones to quote.**

Genuine cross-assay interleaving would need a different procedure
(`FindIntegrationAnchors` + `IntegrateData` treating gene score and expression as two
batches, correcting both directions instead of projecting scATAC onto scRNA). That
answers "can these assays be batch-corrected together", which is not the question
here, and it was not run.

---

## Software

ArchR 1.0.2 (gene score matrix, IterativeLSI), Seurat 5.2.1 / SeuratObject 5.0.2
(`FindTransferAnchors`, `TransferData`), irlba (PCA), uwot (UMAP), FNN (neighbour
statistics), R 4.3.3, conda environment `meso_scatac`.

**Implementation note.** ArchR's own `addGeneIntegrationMatrix()` does not run against
SeuratObject 5.x in this environment — it mixes `CreateSeuratObject()` with the
v3-only `CreateAssayObject()` and the block loop fails with
`slot(object, "features")[[layer]] <- features : more elements supplied than there
are to replace`; setting `options(Seurat.object.assay.version = "v3")` is not
sufficient. The procedure above is therefore written out directly in Seurat, which is
the same algorithm with every parameter visible rather than in package defaults.

Code: `08_joint_embedding.R` (embedding and statistics), `09_figure_joint.R`
(figures).
