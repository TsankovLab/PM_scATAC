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
| 3 | `03_clone_calls.R` | **defines the clones**; runs the split test, clone CNV profiles, driver arms | `out_5Mb/` | `epi_clone_labels.csv`, `epi_clone_decision.csv`, `epi_clone_profiles.rds`, `epi_leaves.csv`, `epi_subclone_summary.csv` |
| 4 | `04_figure_subclones_overview.R` | landscape of all 14 clones + what separates each pair | step 3 | `Plots/epiAneufinder_subclones_overview.pdf` |
| 5 | `05_figure_circular_tree.R` | one circular tree over all 14 clones | step 3 | `Plots/epiAneufinder_circular_tree.pdf` |
| 6 | `06_P4_chr8q_validation.R` | is P4's split really chr8q? chr8p as control | step 3 | `P4_chr8q_epiAneufinder.csv`, `Plots/P4_chr8q_epiAneufinder.pdf` |
| 7 | `07_chromvar_tf_variability.R` | TF variability between clones, per sample and across; plus the mean ArchR gene score of each TF gene over the same cells | step 3, ArchR MotifMatrix + GeneScoreMatrix | `epi_chromvar_*.csv`, `chromvar_z_cache.rds`, `archr_tf_genescore.csv` |
| 8 | `08_figure_tf_variability.R` | three figures: ARI bars + top-22 forest + all 869 motifs ranked per tumour; the square motif-vs-locus variability panel with its top-right quadrant named; and the per-tumour / immune supporting panels | step 7 | `Plots/epi_TF_variability.pdf`, `Plots/epi_TF_genescore.pdf`, `Plots/epi_TF_genescore_supporting.pdf`, `epi_TF_genescore_quadrant.csv` |
| 9 | `09_P4_visium_chr8q.R` | the P4 chr8q clone in space, all four Visium sections | Visium object | `P4_visium_chr8q_spatial.csv`, `Plots/R2_Q3_P4_visium_chr8q.pdf` |

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

**Clones — one distance threshold for the whole cohort.** epiAneufinder does not
estimate a number of subclones: `split_subclones` is literally
`cutree(hclust(dist), k = tree_depth)`, so asking for two gives two in every tumour
whether or not there is anything there. Instead:

1. cluster the cells of one tumour on the Euclidean distance between their 0/1/2 call
   vectors, normalised to a **per-bin RMS** so the scale is the same in every tumour
   (`ward.D2` — Ward's criterion on unsquared distances, which is what `dist()` returns);
2. take the **primary branch point** (k = 2);
3. accept it as a clone boundary only if both branches hold ≥ 20 cells and ≥ 5 % of the
   sample **and** the two branches differ by at least **`ARM_MIN` = 0.55** in mean CNV
   profile on some chromosome arm. Otherwise the tumour is reported as **one clone**.

No cell is discarded either way, and the same threshold is applied to all nine tumours,
so "P8 has one clone" and "P4 has two" are statements on one scale. `c1` is always the
major clone and `c2` the minor one, so every `c2 − c1` reads as *minor clone relative to
the bulk of the tumour*. A clone's CNV profile is, per window, (fraction of its cells
called gain) − (fraction called loss); the separation is that profile difference averaged
over the arm where it is largest, so 0.55 means *"on this arm, 55 % more of one clone's
cells carry the change than of the other's"*.

*Why 0.55.* It is the separation of **P4's chr8q split (0.571)** — the one subclone with
orthogonal validation (step 9, Visium). Setting the bar there means nothing is called
that is less clear-cut than the clone we can independently see in tissue. The next
tumour down is P8 at 0.465, so the threshold sits in a real gap rather than on top of a
sample.

*Why arm level and not genome-wide distance.* Genome-wide distance does not discriminate
at all. The RMS profile difference of the primary split is 0.12–0.19 in **all nine**
tumours and P4 sits at the *bottom* of that range (0.117, the smallest of the nine) —
because P4's difference is a single arm out of 38 while the others are spread thinly
over the genome. Subtracting each tumour's own noise floor does not help either (the
noise-corrected values are 0.112–0.188, P4 still near the bottom). What separates P4 is
*focality*, so the criterion is focal. Averaging within an arm is also what suppresses
the per-bin sampling noise that shallow cells generate — see *Depth* below. Arms with
fewer than 5 informative bins (19p, 22q, …) are excluded, because their mean is one or
two windows and they otherwise dominate the maximum.

*Reported alongside.* `epi_clone_decision.csv` gives, per tumour, the separation, the arm
it sits on, and the mean separation of 20 **random** splits of the same two sizes — the
noise floor for that tumour. Every accepted split clears its own noise floor by more than
tenfold; the four rejected tumours are 2–6× above theirs, i.e. not noise, just not big
enough to call.

*Depth.* epiAneufinder's stage 1 divides each cell by its own mean coverage, so the
per-cell *mean* is normalised — but the *variance* is not, and a call requires a segment
mean to deviate ≥ 50 % from the cell mean (`round()` at 1.5 / 0.5). In a shallow cell the
5 Mb segment means are imprecise, so more segments cross that fixed relative threshold:
Spearman(depth, fraction of non-normal calls) is negative in all nine tumours (−0.12 to
−0.50) and within all clones. This inflates genome-wide distance, not arm-level
coherence, which is a second reason the criterion is arm level. P4 runs *against* the
noise direction — its chr8q clone is the deeper and less-called of the two (median 5 003
vs 4 226 fragments; call rate 0.066 vs 0.073), and its chr8q gain is flat across depth
quartiles (0.718 → 0.666).

*Deeper structure exists and is deliberately not reported.* Applying the same 0.55
threshold recursively (test each branch, then each of its branches) finds further
divisions in P4 (9p, within the chr8q clone), P10 (11p, 8q) and P23 (9q, 4p). The
pipeline stops at the primary split so that every tumour is described on the same
footing and every downstream contrast is a clean two-group comparison.

**TF activity.** chromVAR z from the ArchR `MotifMatrix` (869 cisBP motifs), cached to
`chromvar_z_cache.rds`. Only the five tumours that split have two clones to compare; the
other four are skipped. Per sample: **|Δz|**, the unsigned difference of mean z between
the clones, η², Wilcoxon + BH, and the ARI of an *independent* chromVAR-only k = 2
clustering against the CNV clones. Across samples: median |Δz| with IQR — median, not
mean, because samples differ several-fold in how far apart their clones sit. The FDR
column ranks motifs *within* a sample only; it scales with cell count and must not be
compared across samples.

*Unsigned on purpose.* Which clone is called `c1` is an arbitrary label from `cutree`, so
a signed difference would carry no meaning. Everything reported and plotted is |Δz|.

**TF gene score — the same contrast on a second readout.** Mean ArchR `GeneScoreMatrix`
value of each motif's TF gene, log2(score + 1), computed **per clone** so that
**|Δ gene score| = |c2 − c1|** is the gene-level twin of |Δz|. Cached per (tumour, clone,
gene) to `archr_tf_genescore.csv`, with the between-clone differences in
`archr_tf_genescore_clonediff.csv`. 772 of the 869 motifs have a matching gene, and
because the gene score comes from the same ArchR project as the clone calls, **all five
tumours that split are covered** (P4, P5, P10, P11, P23) over all 20 076 cells.

⚠️ **Two reasons this is a consistency check and not independent evidence.** First, gene
score and chromVAR z are computed from the same fragments. Second, and more concretely,
**a clone defined by a copy-number gain has more fragments across that entire arm**, so
every gene there gains score for a reason that has nothing to do with regulation. That is
measured, not assumed: the median |Δ gene score| of genes **on the arm that defines the
split** against genes elsewhere is 0.094 vs 0.021 (P4, 8q), 0.148 vs 0.020 (P5, 1q),
0.066 vs 0.024 (P10, 11q), 0.189 vs 0.076 (P11, 15q) and 0.014 vs 0.010 (P23, 1q) — a
3–7× inflation in four of the five. Those genes are drawn ringed in panel D, and they are
excluded from the gene-score median in panel E.

AP-1/bZIP + SMARCC1 are flagged **technical**: that axis tracks per-cell fragment count
in this dataset and floors in every comparison. Do not read it as clonal biology.

## Results

**Clones.** 14 clones in 9 tumours (`epi_clone_decision.csv`, `epi_subclone_summary.csv`).
Five tumours pass the 0.55 arm-level threshold and are split in two; four do not and are
reported as a single clone. All nine had a large enough minor branch — the four
rejections are on separation alone, not on cell counts.

| sample | cells | clones | major / minor | arm | separation | noise floor |
|---|---|---|---|---|---|---|
| P1 | 190 | 1 | 105 / 85 | 5q | 0.369 | 0.086 |
| **P4** | 2 832 | **2** | 2 220 / **612** | **8q** | **0.571** | 0.025 |
| P5 | 1 836 | 2 | 1 489 / 347 | 1q | 0.826 | 0.027 |
| P8 | 259 | 1 | 219 / 40 | 6q | −0.465 | 0.060 |
| P10 | 889 | 2 | 599 / 290 | 11q | 0.695 | 0.053 |
| P11 | 1 941 | 2 | 1 481 / 460 | 15q | 0.598 | 0.030 |
| P12 | 686 | 1 | 520 / 166 | 10p | −0.261 | 0.060 |
| P14 | 112 | 1 | 69 / 43 | 4q | −0.230 | 0.118 |
| P23 | 11 331 | 2 | 9 858 / 1 473 | 1q | 0.583 | 0.018 |

"Noise floor" is the mean separation of 20 random splits of the same two sizes. P4 is the
weakest split accepted, by construction.

⚠️ **P5 and P23 both peak at the same chr1 bins** (separation 0.83 and 0.58). A
cross-tumour signature at identical coordinates is what a gene-density / mappability
confound looks like; check those bins before reporting either as a chr1q subclone. P4's
chr8q has no such twin.

**P4 is a genuine chr8q clone** (`P4_chr8q_epiAneufinder.csv`). 612 of 2 832 cells
(21.6 %):

- chr8q gain fraction **0.683 vs 0.116** (Wilcoxon p = 1.3 × 10⁻²⁸⁶)
- chr8p, the control arm, shifts the same way but **~17× less** (0.059 vs 0.026,
  p = 3.0 × 10⁻⁹) — a whole-chromosome or global artefact would move both arms equally
- against a *targeted* split (k-means on the chr8q fraction alone): ARI 0.521,
  precision 0.925, recall 0.634 — i.e. the genome-wide clustering finds the chr8q
  clone **without being pointed at chr8**, and what it finds is almost pure

**Spatial validation** (`P4_visium_chr8q_spatial.csv`). Independent of the scATAC calls,
and run on **all four sections in the Visium object**. All four are the same patient —
the cellranger runs are `37_ST_meso_PT811_Fresh_A1`, `38_..._fresh_B1`,
`39_..._snapFrozen_C1`, `40_..._snapFrozen_D1`, and PT811 = P4 (11 618 malignant scRNA
cells under `sampleID2 == "p811"`, `sampleID == "P4"`). So A1/B1 and C1/D1 are two fresh
and two snap-frozen sections of one tumour, and this is a four-section replication that
also crosses the preservation method.

chr8q-high malignant spots form contiguous territories on every section, not
salt-and-pepper:

| slide | prep | malignant spots | chr8q-amp | Moran's I | p |
|---|---|---|---|---|---|
| A1 | fresh | 1 288 | 38 % | **0.902** | < 0.0033 |
| B1 | fresh | 1 136 | 34 % | **0.794** | < 0.0033 |
| C1 | snap-frozen | 1 343 | 72 % | **0.760** | < 0.0033 |
| D1 | snap-frozen | 1 235 | 31 % | **0.664** | < 0.0033 |

(299 permutations, so p = 0 means below the resolution of the null.) Moran's I 0.664–0.902,
median 0.777; the two **fresh** sections score highest, so the result is not an artefact of
one preservation method.

⚠️ **C1 is the weakest section of the four for the targeted split.** Its k-means puts 72 %
of malignant spots in the "amp" group with the smallest separation of the four
(0.044 vs −0.043), against 31–38 % on the others. The spatial coherence is real on C1, but
if a single representative section is wanted, A1 or B1 is the better choice.

**TF activity does not track the clones.** Testable in the five tumours that split. The
ARI between an independent chromVAR clustering and the CNV clones is ≈ 0 in four of them
(P4 −0.011, P5 −0.007, P10 −0.001, P23 0.007) and 0.114 in P11; mean η² is 0.0003–0.039.
The clones are CNV states, not distinct regulatory states.

Panel C of `Plots/epi_TF_variability.pdf` shows this without any selection: all 869
motifs ranked within each tumour. The curve falls away immediately — max |Δz| is
0.32–1.03 in P4, P5, P10 and P23, and the "most variable" motifs are the handful at the
very start of a curve that is otherwise flat. In P5 and P23 the entire top 5 is the
AP-1/SMARCC1 technical axis. **P11 is the one exception**: max |Δz| 4.47 and mean
η² 0.039, driven by the posterior *HOX* paralogues (HOXC13, HOXD13, HOXB13, HOXA13,
HOXC10) — the same tumour that carries the only non-zero ARI. See the expression check
below before reading much into it.

**Does the TF's own locus move between the clones as well?** Barely
(`Plots/epi_TF_genescore.pdf`, panel D). The Spearman correlation between |Δ gene score|
and |Δz| within a tumour is **−0.03 to 0.11**, and across motifs at the median level
(on-arm genes excluded) it is **0.09**. The motifs that move most are not, in general,
motifs whose TF locus moves. Restricting to the **193 motifs whose locus moves most**
(top quartile of median |Δ gene score|), the ones whose motif also moves most are
**TEAD3, CTCF, TEAD1, NFYC, CEBPD, FOS, JUND, JUN** — TEAD3 and TEAD1 again, alongside the
AP-1 technical axis.

`Plots/epi_TF_genescore.pdf` is the square headline panel: one point per motif, median
over the five tumours, motif variability against locus variability, with **IQR bars over
the tumours** on the named motifs. The shaded top-right quadrant (above the 85th
percentile on **both** axes) holds **28 motifs**, all named and listed in
`epi_TF_genescore_quadrant.csv`: **CTCF, TEAD1, NFYC, CEBPD**, then the AP-1 set
(JUND, JUN, FOSL2, JUNB, FOSB, NFE2, NFE2L2), **TEAD4, ELK4, SNAI1, CEBPB, ID3, ID4, SOX6,
ETV2, ETV3, TCF12, SIX5, ELF3, HIVEP3, BCL6, TFAP4, GLI2, ZNF524**.

⚠️ **The error bars are wide and that is the real message.** Tumours differ several-fold
in how far apart their clones sit, so the IQR across five tumours is large — for most
quadrant members it straddles the cut. Only **8 of the 28** have an entire interquartile
range above the |Δz| threshold (led by TEAD1, TEAD4, SNAI1, CTCF, CEBPB, CEBPD) and only
**12 of 28** above the gene-score threshold. Quadrant membership is a ranking of medians,
not a set of motifs that behave consistently in every tumour; treat the named list as a
shortlist to check, not as a result.

Note what is *not* in the quadrant: **TEAD3**, the single most variable motif overall
(median |Δz| 0.62), sits far left — its locus barely moves between the clones. The two
readouts pick out different TFs, which is the same point the ρ = 0.09 makes.

The supporting file `Plots/epi_TF_genescore_supporting.pdf` carries the per-tumour version
and the immune comparison. It also shows the copy-number effect directly: the ringed
on-arm genes sit clearly to the right (higher |Δ gene score|) but not higher up (their
motifs do not move more). The clones differ in copy number, and that propagates into gene
score without propagating into motif activity.

**The immune/interferon motifs need two statements, not one** (`Plots/epi_TF_median_scatter.pdf`).
Collapsing the five tumours to one median per motif:

- as a **distribution**, the 41 immune/interferon motifs sit slightly higher than the
  other non-technical motifs — median |Δz| **0.101 vs 0.076, Wilcoxon p = 0.0016**
  (all 869 motifs, all five tumours that split);
- as an **enrichment among the leaders**, they are not there — top-decile Fisher
  **OR 1.26, p = 0.59** (it was OR 2.30, p = 0.056 when all nine tumours were forced to
  split, so the earlier trend did not survive).

Both are true: the whole immune distribution is shifted up by about 0.03 z, and no immune
motif is anywhere near the top of the ranking. The highest are RELA, NFKB1/2 and ETS1,
and panel E shows them sitting in the low-expression, low-variability corner with
everything else. Report the shift as what it is — small and distributional — not as
"immune programmes distinguish the clones".

## Files kept

- `out_5Mb/<S>/epiAneufinder_results/` — the CNV calls (the analysis input; 44 MB)
- `frags/` — per-sample fragment files (11 GB; regenerable by step 1)
- `chromvar_z_cache.rds` — 869 motifs × 21 723 cells (128 MB; regenerable by step 7)
- `epi_*.csv`, `P4_*.csv` — all result tables
- `Plots/` — the seven figures
- `archr_tf_genescore.csv` — mean/percent ArchR gene score of 772 TF genes per clone, and
  `archr_tf_genescore_clonediff.csv` — the between-clone differences with each gene's arm
  and whether it sits on the arm that defines the split (both regenerable by step 7; the
  GeneScoreMatrix is pulled from the project once)
- `logs/` — per-step run logs

## Environment

conda `meso_scatac` · R 4.3.3 · epiAneufinder 1.1.5 · ArchR 1.0.2 · chromVAR 1.24.0 ·
ComplexHeatmap 2.18.0 · ggtree 3.10.1 · Seurat 5.2.1 · data.table 1.15.4
