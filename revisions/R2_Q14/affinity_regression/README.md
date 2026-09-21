# Reproducing the TCGA MESO BAP1 inferred-TF-activity analysis (affinity regression)

Reviewer context: R2_Q14. Target of the reproduction is **supplementary table S2B** of
Hmeljak *et al.*, *Cancer Discovery* 2018 (TCGA MESO; file
`../21598290cd180804-sup-205173_2_supp_5073240_pg14sl.xlsx`, tab `2B_bap1_inferred_TF_activity`)
and its main-text Figure 2D–G: TF activities inferred with the method of
Osmanbeyoglu *et al.*, *Nat Commun* 2017 (PMID 28139702), compared between BAP1-inactivated and
BAP1 wild-type tumours by t-test.

---

## 1. The method being reproduced

Affinity regression (AR; Pelossof *et al.* 2015) models the tumour × gene expression matrix as a
bilinear function of TF binding and protein levels:

```
Y  ≈  D · W · Pᵀ
Y  genes × tumours     mean-centred log10(RSEM + 1) expression
D  genes × TFs         binary TF → target-gene prior (motif hits in promoters)
P  tumours × proteins  mean-centred RPPA
W  TFs × proteins      learned interaction matrix
```

- **TF activity per tumour** = `W · Pᵀ` (TFs × tumours).
- W is fitted on the transformed system `Yᵀ Y ≈ (Yᵀ D) W P`, reduced by truncated SVDs of `YᵀD`
  and `Pᵀ`, and solved with the SLEP `LeastR` solver (accelerated proximal gradient with an
  L1 penalty λ; the published text calls the fit "ridge", the released code applies an L1 term
  with `rsL2 = 0`).
- Hyperparameters (λ, spectrum fractions for the two SVDs) were chosen in the original work by
  10-fold cross-validation on expression reconstruction.
- The MESO paper then compared inferred activity between BAP1-inactivated and wild-type tumours by
  t-test with Benjamini–Hochberg correction (FDR < 0.01 called significant).

## 2. What was recovered exactly from the published material

| component | recovered as | evidence |
|---|---|---|
| TF → target prior | **MSigDB v6.1 `c3.tft` (TRANSFAC)** | the paper's published YY1 (tab 2C) and IRF8 (tab 2D) target lists are reproduced **exactly** by `GCCATNTTG_YY1_Q6` (427/427) and `ICSBP_Q6` (248/248); v2023.2 and v6.0 do not match |
| motif chosen per TF | **the largest target set** among a TF's motifs | the rule that yields both published lists (e.g. YY1 has four motifs; the largest is the one published) |
| TF panel | the paper's **141 TFs** (tab 2B), used as given | — |
| BAP1 labels | the paper's own integrated calls (tab 2A) | 26 inactivated, 23 no inactivation, 5 possibly inactivated among tumours with RPPA |
| data | TCGA MESO from cBioPortal: PanCancer Atlas and Firehose legacy releases | RPPA available for 63 tumours in both |
| algorithm | line-by-line Python port of the Osmanbeyoglu lab MATLAB code (`code/Affreg`: `ar_train.m`, `ar_model2w.m`, `ar_predict.m`, `ar_reconstruction.m`, `SLEP/LeastR.m`) | see §3 |

### TF-to-motif crosswalk (`prepared/tf_motif_crosswalk.csv`)

TRANSFAC matrix names are not gene symbols (IRF8 = `ICSBP`, NFE2L1 = `TCF11`, RELA =
`NFKAPPAB65`, SPI1 = `PU1`). The paper's names come from an automated alias lookup, including
some that are biologically wrong but must be kept to reproduce the table (POSTN from the `OSF2`
motif, ACTR1A from `ARP1`, ADD1 as an SREBF1 alias, SF1 and NF1).

- **121 / 141** TFs map automatically by HGNC symbol / alias / previous symbol.
- **19** are assigned by hand from TRANSFAC nomenclature, each with a confidence label
  (high / medium / low). Low-confidence assignments: CREBBP (shares `P300_01` with EP300), NFYB,
  NFATC4, RUNX2 (generic `AML_Q6`), MEF2D.
- **RARA** could not be mapped — v6.1 has no RAR motif — and is absent from the reproduction
  (140 TFs compared).

## 3. Implementation and its validation

MATLAB could not be used: the cluster licence does not include this user. `affreg.py` is a direct
port; every column-major MATLAB operation (`vec`, `reshape`, the diagonal-equation removal) is
reproduced with `order='F'`, and `LeastR` keeps the lab copy's `+0.001` line-search tolerance.

Validation on noiseless simulated data with a known W (`python affreg.py`, mirroring the lab's
`run_helloworld.m`):

| λ | spectrum A/B | corr(W, Ŵ) | held-out Y reconstruction | TF activity |
|---|---|---|---|---|
| 0.01 | 1 / 1 | 0.983 | 1.000 | 0.982 |
| 0.1 | 0.95 / 0.9 | 0.843 | 0.880 | 0.868 |

## 4. Choices the methods do not state, and how they were made

| choice | options tested | how decided |
|---|---|---|
| gene universe for Y | top 5,000 most variable genes (2017 methods) / all expressed genes with ≥1 motif hit | CV |
| data release | PanCancer Atlas / Firehose legacy | CV |
| D normalisation | binary / column unit-norm | CV |
| Y normalisation | none / per-tumour unit-norm (applied in the lab's run scripts) | CV |
| training tumours | all 63 with RPPA / the 54 with a tab-2A label | CV |
| λ, spectrum A, spectrum B | 8 λ values × 5–6 spectrum pairs | CV |
| t-test | Welch / Student | makes no difference (checked) |
| "possibly inactivated" | excluded / grouped with inactivated | makes little difference (checked); excluded, as in the paper's two-group comparison |

**Selection rule, fixed before the complete grids were inspected (`select_final.py`):** the
configuration with the highest 10-fold CV reconstruction correlation — the original authors'
criterion. **Agreement with tab 2B is never used to choose.** Also reported: the best-CV
configuration under the most literal reading of the published methods (binary D, no
normalisation, top 5,000 genes, all tumours), and the spread across the 20 best CV configurations.

## 5. Results

Both grids completed without error (760 configurations with `rsL2 = 0`). Cross-validation is
flat: held-out reconstruction ranges only r = 0.299–0.304 across the best 20 configurations, so the
selection separates candidates weakly. Comparison is 26 BAP1-inactivated vs 23 wild-type
tumours, 140 TFs (RARA unmapped), Welch t-test, BH.

| | **primary (best CV)** | literal methods |
|---|---|---|
| configuration | PanCan release, top 5,000 genes, all 63 tumours, **column-normalised D**, no Y normalisation, spectrum 0.8/0.8 (λ inert) | PanCan release, top 5,000 genes, all 63 tumours, binary D, no Y normalisation, spectrum 0.9/0.9 (λ inert) |
| CV reconstruction r | 0.304 | 0.298 |
| Spearman, effect estimates vs tab 2B | **0.54** | 0.51 |
| Spearman, signed −log p | 0.56 | 0.54 |
| same direction | 67 % | 66 % |
| significant, FDR < 0.01: ours / paper | 35 / 28 | 20 / 28 |
| **paper-significant TFs recovered** (FDR < 0.01 in both, same sign) | **12 / 28** | 8 / 28 |
| EGR2 (paper: up, #2) | **up, #2**, FDR 2 × 10⁻⁵ | up, #12 |
| IRF8 (paper: up, #1) | up, #26, FDR 0.002 | up, #14 |
| YY1 (paper: down, FDR 3 × 10⁻⁵) | down, #64, FDR 0.17 (ns) | up, #117, FDR 0.8 (ns) |
| MAX (paper: down, #6) | down, #18, FDR 4 × 10⁻⁴ | down, #22, FDR 0.017 |

Recovered in the primary configuration (FDR < 0.01 in both, same direction): IRF1, EGR2, PAX8, TFAP2A, FOXO1, SP3, IRF8 (up in inactivated); MAX, STAT4, ATF3, CEBPG, SMAD4 (down). A further 2 (NR3C1, TTF1) are significant in both but in opposite directions; the volcano's label count (14) includes them.

**Spread across the 20 best CV configurations** (min / median / max): Spearman 0.40 / 0.50 / 0.54;
paper-significant TFs recovered 4 / 7 / 12; **IRF8 rank 1 / 3.5 / 26; EGR2 rank 1 / 2.5 / 26;
YY1 rank 31 / 64 / 118**.

Reading: the reproduction recovers the rank structure of tab 2B and its two headline TFs.
EGR2 and IRF8 are higher in BAP1-inactivated tumours in 20/20 and 20/20 of the best-CV
configurations (and in 750/760 and 720/760 of all grid points), and near the top in most. MAX, CEBPG, ATF3 and STAT4 reproduce as lower. **The YY1 decrease, emphasised in the paper, does not reproduce reliably:** YY1 is lower in only
13 of the 20 best-CV configurations (ranks 31–118), is not significant in the CV-selected fit, and
is slightly higher (ns) under the literal configuration. It is the least
robust of the named findings. (It does reach p = 0.0017 under per-tumour Y normalisation,
`results/provisional/`, a configuration chosen for inspection and not by CV.)

Discordant TFs: **E2F1** is the strongest hit in the reproduction but not significant in the paper.
**NR3C1** and **TTF1** (TRANSFAC's NKX2-1 motif) go the opposite way.

## 6. Label-permutation null

Agreement could in principle arise from structure in the inferred activities alone. BAP1 labels
are permuted among the 49 labelled tumours (2,000 times) and the t-test and agreement statistics
are recomputed on the *same* activity matrix (`permutation_null.py`).

| statistic | primary: observed | null mean (SD) | null max | p | literal: observed | null max | p |
|---|---|---|---|---|---|---|---|
| Spearman, estimates | 0.542 | 0.00 (0.21) | 0.516 | 0.0005 | 0.513 | 0.508 | 0.0005 |
| Spearman, signed −log p | 0.564 | 0.00 (0.22) | 0.521 | 0.0005 | 0.536 | 0.509 | 0.0005 |
| paper-significant TFs recovered | 12 | 0.0 (0.1) | 4 | 0.0005 | 8 | 3 | 0.0005 |
| same direction | 67 % | 50 % (8) | 72 % | 0.016 | 66 % | 71 % | 0.012 |

The agreement is driven by BAP1 status. The null is nonetheless wide (SD 0.2): permuted labels
occasionally reach ρ ≈ 0.5, because the 140 inferred activities are strongly correlated with each
other. A single agreement coefficient should therefore not be over-read. The recovery of 12
significant TFs (null maximum 4) is the more discriminating statistic.

## 7. Caveats

- Without the original code run, the original parameter values and the paper's supplementary
  methods (not retrievable — PMC serves a download gate), an exact numerical match is not expected;
  the question is whether the ranking and the named findings reproduce.
- The crosswalk is reconstructed; low-confidence TFs should not be over-interpreted individually.
- Nine tumours with RPPA have no tab-2A label; they contribute to fitting W (in the `all63`
  configurations) but not to the t-test.
- Activity scale differs from the paper's (units depend on normalisation); estimates are compared
  by rank, not magnitude.
- The published supplementary table lists 141 TFs but reports no sample numbers; the 26 vs 23
  split is inferred from tab 2A ∩ tumours with RPPA.

## 8. Files

| file | role |
|---|---|
| `affreg.py` | Python port of affinity regression + simulation test |
| `build_inputs.py` | crosswalk, D, Y, P for each release × gene universe |
| `mesobap1.py` | shared loading, fitting, t-test, comparison with tab 2B, CV |
| `reproduce_quick.py`, `reproduce_sensitivity.py`, `scan_ynorm.py` | exploratory diagnostics (not used for selection) |
| `cv_grid.py`, `cv_grid2.py` | cross-validation grids (LSF) |
| `select_final.py` | pre-specified configuration selection |
| `report_tfactivity.py` + `plot_tfactivity.R` | final activity matrix, 2B-equivalent table, figures |
| `permutation_null.py` | label-permutation null |
| `prepared/` | inputs, crosswalk, grid results, selection |
| `results/<tag>/` | `tf_activity.csv`, `bap1_tf_ttest.csv`, `summary.json`, `permutation_null.json`, PDFs |
| `code/` | cloned Osmanbeyoglu lab repositories (reference) |
| `data/` | cBioPortal TCGA MESO study bundles |
| `c3.tft.v6.1.symbols.gmt`, `hgnc_complete_set.txt` | TF prior and symbol table |
