# R2_Q14 — TF activity vs BAP1 status (malignant compartment)

**Reviewer 2, Q14:** systematically explore whether TF activity (chromVAR)
differs between mesothelioma tumors **with vs without BAP1 loss**, and whether
TF activity increases/decreases in relation to low BAP1 (gene score). Reviewer
noted a trend for higher RUNX1/RUNX2 activity in BAP1-retained samples.

## Data
- ArchR object: `../tumor_compartment/scatac_ArchR`
- TF activity: `MotifMatrix` chromVAR **deviations** (`assays(mSE)[[1]]`),
  aggregated by mean per sample — the same metric plotted as `dev_diff` in
  `discover_oncogenic_TFs.R` (raw deviations, unscaled). Base-R `scale()` is
  applied only in the per-cell `wilcoxauc` supplement, mirroring the original.
- BAP1 accessibility: `GeneScoreMatrix` (BAP1 row)
- Sample label: `Sample` cellColData column
- **BAP1 status (per-sample genetic annotation):**
  - LOST (n=5): P3, P4, P5, P8, P10
  - RETAINED (n=6): P1, P11, P12, P13, P14, P23

## Approach (systematic)
Cell counts per tumor are highly unbalanced (P23 = 11,687 vs P13 = 18), so a
per-cell test is pseudoreplicated and dominated by one sample. **The primary
analysis uses the sample as the unit of replication (pseudobulk).**

- **(A) Group test** — per-sample mean TF activity, LOST vs RETAINED, Wilcoxon
  rank-sum + Welch t-test, BH-adjusted. → `group_test_lost_vs_retained.csv`,
  `Plots/volcano_group_lost_vs_retained.pdf`
- **(B) Continuous** — Spearman correlation of each TF's per-sample activity
  with the sample's BAP1 gene score across the 11 tumors (rho < 0 ⇒ activity up
  as BAP1 down). → `correlation_TFactivity_vs_BAP1genescore.csv`,
  `Plots/correlation_TFactivity_vs_BAP1.pdf`
- **Supplement** — per-cell presto Wilcoxon (flagged as pseudoreplicated),
  columns `*_cell` in the master table.
- **Validation** — BAP1 gene score is lower in the genetically BAP1-lost
  samples. → `Plots/BAP1_genescore_by_status.pdf`

## Key findings
- **BAP1 gene score validates the genetic annotation:** BAP1-lost tumors have
  lower BAP1 accessibility than BAP1-retained.
- **RUNX confirmed and generalized:** RUNX1/RUNX2/RUNX3 (and cofactor CBFB) TF
  activity is **higher in BAP1-retained** tumors — i.e. RUNX activity *decreases*
  with BAP1 loss. Sample-level t-test p: RUNX3 = 0.0027, RUNX2 = 0.0053,
  RUNX1 = 0.085; Wilcoxon p = 0.017 for RUNX2/RUNX3. Spearman rho vs BAP1 gene
  score positive (RUNX3 0.61, RUNX2 0.47, RUNX1 0.27).
- **Other TFs up in BAP1-retained:** SOX2/SOX21, POU5F1, TCF21, HOXB9, ZBTB4.
- **TFs up in BAP1-lost / rising as BAP1 falls:** CTCF/CTCFL, RXRA/RXRB/RARG,
  NR1I3, NR4A1/NR4A3, TP73, ZNF784, TBR1, EOMES.

## Bulk RNA-seq corroboration (Bueno, TCGA, MESOMICS)
Independent validation in the three bulk cohorts used elsewhere in the study
(loaded from `bulkRNA_meso/bulk_RNA_studies*.rds`, built by
`git_repo/bulk_analysis.R`). Bulk has expression, not chromVAR activity, so the
analog is the **RUNX–BAP1 relationship**: scATAC predicts RUNX co-varies
*positively* with BAP1 (RUNX down when BAP1 down/lost).

- **RUNX ~ BAP1 expression correlation (Spearman):** positive in **11/12**
  gene×cohort cells (only MESOMICS CBFB −0.07, ns). Strongest in Bueno (n=211):
  RUNX1 rho=0.29 (p=2.7e-5), RUNX2 rho=0.30 (p=9e-6), RUNX3 0.14 (p=0.036),
  CBFB 0.25 (p=3e-4). TCGA: RUNX2 0.24 (p=0.026), CBFB 0.23 (p=0.032).
- **Subtype-adjusted `lm(RUNX ~ BAP1 + subtype)`:** the positive BAP1 term is
  significant for RUNX1 and/or RUNX3 in **all three** cohorts (e.g. MESOMICS
  RUNX1 p=0.009, RUNX3 p=0.011; TCGA RUNX1 p=0.017, RUNX2 p=0.010), i.e. the
  association is not merely driven by histology.
- **MESOMICS BAP1 IHC status** (real BAP1 loss; direction verified against BAP1
  mRNA, retained=6.07 vs lost=5.34, p=4e-6): RUNX1/RUNX3 trend lower in
  BAP1-lost (p=0.11 / 0.13), RUNX2 flat.
- **TCGA/Bueno BAP1-mRNA tertiles (low = loss proxy):** RUNX lower in BAP1-low in
  **all 8** gene×cohort tests (every delta negative); Bueno RUNX1 p=1.7e-4,
  RUNX2 p=4e-4, CBFB p=0.011.

**Conclusion (unadjusted):** across three independent bulk cohorts the raw
RUNX–BAP1 direction matches the scATAC — RUNX (esp. RUNX1/RUNX2) is reduced with
BAP1 loss/low BAP1. BUT this is substantially confounded by histology — see next
section, which is the important qualification.

## Histology confound — is the RUNX–BAP1 link real or histology-driven?
(`R2_Q14_histology_disentangle.R`) In MPM, BAP1 loss is enriched in
**epithelioid** tumors while RUNX1/RUNX2 are mesenchymal/EMT TFs higher in
**sarcomatoid** — so "RUNX low when BAP1 lost" could be pure histology. Tested
directly, using the study's continuous sarcomatoid score (cNMF20) as the
histology axis (available for both bulk metadata and the scATAC samples).

**(1) BAP1 loss DOES track histology — strongly.**
- BAP1 mRNA is lowest in epithelioid, higher in sarcomatoid: Bueno KW p=3.5e-9
  (epi 3.15 vs sarc 4.1, p=2.9e-6); MESOMICS KW p=0.01. BAP1~sarc-score rho
  positive in all cohorts (Bueno 0.25 p=2e-4).
- MESOMICS real IHC status × subtype: **52/63 BAP1-lost cases are epithelioid**
  (only 2 sarcomatoid); retained is balanced (18 epi / 11 sarc / 14 biphasic).
  Fisher p = 2.3e-5.
- scATAC (9 samples with a sarc score): BAP1 gene score vs sarc score rho=0.48;
  mean sarc score lost=0.06 vs retained=0.45 — the BAP1-retained scATAC samples
  are markedly more sarcomatoid (P1, P14, P13 drive it).

**(2) Does RUNX~BAP1 survive histology adjustment? Partly — weaker & gene/cohort-inconsistent.**
- **Epithelioid-only:** in Bueno the RUNX1/RUNX2 correlations collapse to ~0
  (−0.02, −0.03) — Bueno's raw signal was largely histology. They persist in
  **TCGA** (RUNX1 0.25, RUNX2 0.24, p≈0.07) and partly MESOMICS (RUNX1 0.23).
- **Partial Spearman | sarc score:** attenuated but still positive/sig for
  RUNX1/RUNX2 in Bueno (0.17*, 0.20**) and RUNX2 in TCGA (0.23*).
- **`lm(RUNX ~ BAP1 + sarc_score)`:** BAP1 term still significant for RUNX1/RUNX2/
  CBFB in Bueno (p=0.011/0.007/0.006), RUNX2/CBFB in TCGA (0.013/0.025), RUNX3 in
  MESOMICS (0.016). No BAP1×sarc interaction.
- **scATAC (n=9):** raw RUNX3 activity ~ BAP1 gs rho=0.62 (p=0.09) drops to
  partial 0.44 (p=0.24); RUNX1/RUNX2 raw 0.35 flip to −0.19 after adjustment.
  With 9 samples BAP1 and histology are not separable.

**Bottom line for the reviewer:** BAP1 loss is strongly confounded with
epithelioid histology, and a large part of the raw RUNX–BAP1 signal (especially
in Bueno and in the 9 scATAC samples) is explained by the sarcomatoid skew of the
BAP1-retained group. A **weaker, histology-independent residual** association
survives sarc-score adjustment for RUNX1/RUNX2 in the larger bulk cohorts
(Bueno, TCGA), so the effect is not entirely histology — but it is smaller and
less consistent than the unadjusted numbers suggest. The honest framing is:
"the RUNX trend is real at the descriptive level but is largely, though not
wholly, attributable to the epithelioid/sarcomatoid axis; it cannot be cleanly
separated from histology in the 9-sample scATAC set."

### Synthesis figure (`R2_Q14_histology_synthesis.R`)
One 6-panel figure that makes the confound visible directly, in both modalities.
Top row = bulk RNA (BAP1 & RUNX = expression); bottom row = scATAC (BAP1 = gene
score, RUNX = mean RUNX1/2/3 chromVAR activity).
- **A/D** BAP1 rises along the sarcomatoid-score axis (BAP1 mRNA lowest in
  epithelioid; scATAC BAP1 gene score higher in the more sarcomatoid samples).
- **B/E** RUNX rises along the *same* axis — very strongly (bulk RUNX2~sarc
  ρ=0.57/0.66/0.41; **scATAC RUNX activity~sarc ρ=0.82, p=0.011**).
- **C/F** the BAP1↔RUNX scatter is coloured by the continuous sarc score, so the
  "positive BAP1–RUNX correlation" is seen to be a sarcomatoid-score gradient
  (points sort by colour along the trend line).
Numbers behind the panels: `histology_synthesis_correlations.csv`. Plots:
`Plots/R2_Q14_histology_synthesis.pdf` (full 6-panel) and
`Plots/R2_Q14_confound_BAP1_RUNX_by_sarcscore.pdf` (just panels C+F).

Disentangle files: `histology_BAP1_association.csv`,
`RUNX_BAP1_histology_adjusted.csv`, `scATAC_histology_disentangle.csv`,
`scATAC_per_sample_RUNX_BAP1_sarc.csv`; plots `Plots/BAP1_by_subtype_bulk.pdf`,
`Plots/RUNX_vs_BAP1_epithelioid_only.pdf`,
`Plots/RUNX_BAP1_histology_adjusted_heatmap.pdf`,
`Plots/scATAC_BAP1_vs_sarcscore.pdf`.

## Within-tumor & genetic analysis — is BAP1 status linked to ANY TF, independent of histology?
(reviewer follow-up) The between-tumor RUNX–BAP1 link is largely histology
(above). To test BAP1 **independent of histology** we (1) compare cells **within a
single tumor** (same histology) — BAP1-high vs BAP1-low, across **all** TFs; (2)
use **inferCNV** on the matched scRNA to ask whether BAP1-retained subclones even
exist (reviewer suggested P8); (3) **match** the scRNA subclones to scATAC; and
(4) use **genetic** (not just mRNA) BAP1-loss calls in the bulk cohorts.

### Does a BAP1-retained subclone exist? (inferCNV, reused existing run)
Reused `tumor_compartment/scrna/infercnv/run.final.infercnv_obj` (P8 = 4,109
cells + 1,124 normal-pleura reference). Per-cell chr3p / BAP1-window CNV,
normalized to each cell's genome-wide median (`persample_chr3p_bimodality.csv`,
`scRNA_persample_clonality.csv`):
- **P8 is CLONALLY BAP1-deleted** — chr3p median 0.92, **73 % of cells "lost" at
  the BAP1 window, only 3 % retained**, one coarse subcluster. In the few
  CNV-"retained" P8 cells BAP1 mRNA is not elevated (p=0.25). **So the
  BAP1-retained subclone the reviewer hoped to find in P8 does not exist — P8's
  BAP1 loss is truncal/clonal.**
- **P11 is genuinely SUBCLONAL** — 47 % chr3p-retained / 29 % lost (bimodal), and
  the CNV calls are **validated**: CNV-retained P11 cells express more BAP1 mRNA
  (0.140 vs 0.102, Wilcoxon **p=0.012**). 1,021 scRNA + 1,980 scATAC cells.
- P5 = clonal-lost; P4 = mixed/intermediate (chr3p arm retained but BAP1 focal —
  its annotated loss is likely mutational, invisible to CNV); P1/P12/P13/P14/P3 =
  clonal-retained.

### Within-tumor TF activity (scATAC, all TFs) — BAP1-high vs BAP1-low
Per sample, split cells BAP1-accessible (gene score>0) vs BAP1-silent (=0),
presto `wilcoxauc` on scaled chromVAR deviations across **all ~870 motifs**
(`within_sample_TFactivity_BAP1_high_vs_low.csv`, `within_sample_TF_meta.csv`):
- **RUNX shows NO consistent within-tumor direction.** RUNX1 is actually *lower*
  in BAP1-high cells in P5 (AUC 0.42, padj 3.8e-6) and P11 (AUC 0.41, padj
  3.4e-10), but RUNX3/CBFB are slightly *higher* in BAP1-high in P8/P12. Removing
  between-tumor histology removes the RUNX–BAP1 relationship.
- The TFs most **consistently** tracking within-tumor BAP1 accessibility are
  **IRF3/5/9, CEBPE/G/Z, ETS (ELF4/EHF/ERF/SPDEF), and generic YY1/NFYA/NFYB/
  CREB1/ATF1/YBX1** — but with **small effects (AUC−0.5 ≈ 0.05)**, consistent with
  a global chromatin-openness gradient rather than a BAP1-specific program. **Not
  RUNX.**
- In the one clean subclonal tumor (**P11**), the BAP1-retained subclone has
  higher GATA4/2/3/6, TEAD1/4, CEBP, STAT1/3/5 and **lower RUNX1/PKNOX1/TGIF2** —
  a broad regulatory difference between subclones, again not a "RUNX-up" program.

### Cross-modal match
Per-sample scRNA CNV-retained fraction vs scATAC BAP1-accessible fraction:
Spearman **rho=0.48** (positive, ns at n=9) — the two modalities give a broadly
concordant picture of subclonal BAP1 status (`crossmodal_scRNA_vs_scATAC_BAP1.csv`).

### Genetic BAP1 loss in the bulk cohorts
Genetic (DNA/protein) BAP1-loss calls **do exist** and were used
(`bulk_genetic_*` files):
- **Bueno `FISH.chrom3`** (chr3p deletion) — but 44/45 called samples are "del 3p/3"
  (near-universal), so it cannot support a Loss-vs-Retained contrast.
- **TCGA `STATUS_3P`** (arm CNV) — only 5 "Lost" (underpowered) + per-sample BAP1
  segment in `TCGA_CNV_hg38.rds`.
- **MESOMICS `IHC.BAP1`** — 63 loss / 43 retained (best-powered; direction VERIFIED
  here: "NO"=loss, BAP1 mRNA 5.34 vs 6.07 — note the old `R1_Q3_F26` header had this
  reversed).
- **Genetic BAP1 loss strongly tracks epithelioid** (MESOMICS Fisher **p=5e-4**,
  83 % of losses epithelioid vs 42 % of retained) — the histology confound holds
  with genetic calls, not just mRNA.
- **RUNX by genetic BAP1:** RUNX1 trends down in loss overall (all 3 cohorts,
  ns) but **collapses within epithelioid** (MESOMICS delta_epi −0.12, p=0.68).
- **MESOMICS within-epithelioid genome-wide TF scan** (52 loss vs 18 retained):
  785 TFs, 95 nominal p<0.05, but **none survive FDR (padj<0.1)** and **RUNX is ns**;
  top nominal hits are FOXB1/E2F8/HOX (down in loss) and ZNF219/IRF7/SOX11 (up) —
  not RUNX.

### Full TF differential EXPRESSION scan by BAP1 status, all cohorts
`R2_Q14_bulk_TF_DE_by_BAP1.R` runs a Wilcoxon DE (BAP1-Loss vs Retained) for
**every one of the 869 chromVAR-motif TFs**, per cohort, both **overall** and
**within epithelioid**, plus a sarc-score-adjusted lm
(`bulk_TF_DE_by_BAP1_percohort.csv`, `bulk_TF_DE_by_BAP1_crosscohort.csv`;
volcano `Plots/bulk_TF_DE_by_BAP1_mesomics_volcano.pdf`, RUNX panel
`Plots/bulk_TF_DE_RUNX_overall_vs_epi.pdf`):
- **Bueno** (44 loss / 1 retained) and **TCGA** (5 loss / 44 retained) are too
  degenerate/underpowered — 0 TFs survive FDR in either; TCGA's top nominal hits
  are KLF12/NR2F2/HOXB6 (up) and ETV7/BHLHE40 (down).
- **MESOMICS** (63 loss / 43 retained, the powered cohort) yields **28 TFs at
  BH<0.05 overall**: up in loss = **BARX2, SOX11, POU3F3, SOX6, VTN, HMX1**;
  down in loss = **FOXB1, MYBL1, NR3C1, ARNTL2, GLIS3, NFE2L3**. These are
  neural-crest / mesenchymal-patterning TFs — i.e. the **sarcomatoid axis**, not a
  BAP1-repression signature.
- **Every one of the 28 MESOMICS hits loses FDR significance within epithelioid
  (0 survive)** → the differential expression is the histology confound, not
  BAP1 itself.
- **RUNX1/2/3/CBFB are non-significant in all three cohorts**, overall and
  within-epithelioid (MESOMICS RUNX1 p=0.11 overall, p=0.68 within-epithelioid).
- Cross-cohort, **only NPAS3 and PRRX2** are consistent-direction and nominally
  significant in ≥2 cohorts (both up in loss, tcga+mesomics) — neither survives FDR
  nor is a canonical BAP1 target.

### BAP1 EXPRESSION as proxy, within epithelioid (bypassing genetic calls)
Because the genetic calls are degenerate/underpowered, `R2_Q14_BAP1expr_proxy_epithelioid.R`
uses **continuous BAP1 expression** to split BAP1-high vs -low **within epithelioid
tumors only** (histology held ~constant), in bulk RNA and in low-scS scATAC samples.
- **Bulk RNA (BAP1 mRNA proxy, epithelioid-only):** unlike the genetic scan, this IS
  powered — TCGA 42 / MESOMICS 15 TFs at BH<0.05 (Bueno 0), and there is a
  **cross-cohort-consistent BAP1-proxy TF signature** (same-sign in all 3 cohorts,
  sig in ≥2; `epi_proxy_bulk_consistentTFs.csv`):
  - *positively* correlated with BAP1 (down in BAP1-low): **SMARCC1** (SWI/SNF, a BAP1
    chromatin partner), the **HOXA1-5 cluster**, ZNF589, MITF, TCF4, E2F8, NR1D1, PPARD.
  - *negatively* (up in BAP1-low): **SOX11, BARX2, POU3F3**, SIX4, ZNF219, FOXD2,
    PRDM16, ID1, MAFA — note **SOX11/BARX2/POU3F3 replicate the direction of the
    genetic-loss scan** (up in loss), so they are the most robust BAP1-low TFs.
  - **RUNX1/2/3/CBFB are non-significant in every cohort** (all padj>0.28) — RUNX is
    absent from the BAP1-proxy signature even with this more powered design.
  Files: `epi_proxy_bulk_TFvsBAP1_spearman.csv`, `epi_proxy_bulk_TF_DE_lowVShigh.csv`,
  plot `Plots/epi_proxy_bulk_spearman_volcano.pdf`.
- **scATAC (BAP1 gene-score proxy, epithelioid/low-scS samples P4,P5,P8,P11):** the
  BAP1-high vs -low TF-activity meta is dominated by a **generic ETS/NFY/YY1 openness
  axis** (NFYA/B, ELF3/4, EHF, ETS1/2, ELK1/3, YY1, SPDEF — all 4/4 samples, small
  AUC−0.5≈0.07), consistent with BAP1 gene score tracking overall accessibility;
  **RUNX1 is mildly LOWER in BAP1-high cells** (flips vs the between-tumor trend).
  This scATAC proxy does **not** corroborate the bulk signature (cross-modal Spearman
  **rho=0.04, p=0.29**) — i.e. the single-cell BAP1-accessibility split is confounded
  by chromatin openness and is uninformative for BAP1-specific TF biology.
  Files: `epi_proxy_scATAC_TFactivity_meta.csv`, plots
  `Plots/epi_proxy_scATAC_meta_bars.pdf`, `Plots/epi_proxy_crossmodal_bulk_vs_scATAC.pdf`.
- **Caveat:** within-epithelioid BAP1 mRNA may still track a residual
  sarcomatoid-within-epithelioid gradient (SOX11/BARX2/ID1/PRDM16 are mesenchymal);
  the SMARCC1/HOXA arm is the more BAP1-specific-looking axis. Either way the key point
  holds: **a reproducible non-RUNX BAP1-proxy TF signature exists in bulk, RUNX is not
  part of it, and it does not reproduce in the single-cell openness-confounded scATAC.**

#### Does the bulk signature survive removing the residual sarc gradient? YES.
`R2_Q14_BAP1proxy_partial_sarc.R` computes the **partial Spearman** of each TF vs
BAP1 mRNA **controlling for the cNMF20 sarcomatoid score**, within epithelioid only,
per cohort (`epi_proxy_partial_sarc_percohort.csv`, `_crosscohort.csv`; plots
`Plots/epi_proxy_partial_marg_vs_partial.pdf`, `Plots/epi_proxy_partial_signature_dumbbell.pdf`):
- **The residual sarc confound within epithelioid is weak** — BAP1 mRNA ~ sarc-score
  rho is only 0.26 (Bueno, p=0.06), 0.14 (TCGA), 0.05 (MESOMICS). So there is little
  histology left to remove once epithelioid is imposed.
- **The signature is essentially unchanged after partialling.** Marginal-significant
  hits that survive sarc-adjustment: **TCGA 42→41 (98%)**, **MESOMICS 15→15 (100%)**;
  the marg-vs-partial rho scatter sits **on the diagonal** for every signature TF.
- **The cross-cohort survivors** (same-sign partial rho, padj_part<0.05 in ≥2 cohorts)
  are exactly the reported signature: **ZNF589, SMARCC1, HOXA1-5, CPEB1** (up with BAP1)
  and **SIX4** (up in BAP1-low) — e.g. MESOMICS SMARCC1 marg 0.545 → partial 0.543,
  HOXA4 0.652 → 0.653; TCGA HOXA1 0.612 → 0.606.
- **RUNX stays null after partialling** in all cohorts (partial p: MESOMICS RUNX1 0.03
  nominal only, TCGA RUNX1/2 ~0.12-0.14, Bueno ns; no cross-cohort survival).
- **Conclusion:** the **SMARCC1 + HOXA + ZNF589** BAP1-proxy program is a **genuinely
  sarc-independent** BAP1-associated signature, not residual histology — this resolves
  the caveat above and is the substantive "any TF beyond RUNX linked to BAP1?" answer.
  RUNX is specifically not part of it.

**Bottom line for the reviewer:** the specific P8-subclone strategy is not viable
because P8's BAP1 loss is **clonal/truncal** (inferCNV: no retained subclone).
Using the one genuinely subclonal tumor (P11), a systematic within-tumor scan
across every TF, cross-modal validation, and **genetic** BAP1 calls in the bulk
cohorts, **BAP1 status does not drive a specific, reproducible TF-activity program
once histology is held constant** — and RUNX in particular shows **no** within-tumor
or within-epithelioid association. This is the strongest evidence that the
descriptive RUNX–BAP1 trend is an epithelioid/sarcomatoid histology effect, not a
cell-intrinsic consequence of BAP1 loss. The only within-tumor TF signals are
small/generic (chromatin openness) or subclone-specific (GATA/TEAD in P11) and are
neither FDR-robust nor cross-cohort consistent.

Files: `R2_Q14_within_sample_TFactivity.R` (scATAC within-sample, all TFs),
`R2_Q14_scRNA_subclones.R` (inferCNV subclones + BAP1-mRNA validation + cross-modal),
`R2_Q14_bulk_genetic_BAP1.R` (genetic BAP1 calls), `R2_Q14_within_sample_figs.R`
(headline figure), `probe_p8*.R` (recon). Key plots:
`Plots/R2_Q14_within_sample_headline.pdf` (A clonality violin, B cross-modal,
C within-P11 volcano, D RUNX by tumor), `Plots/within_P11_TF_volcano.pdf`,
`Plots/scRNA_CNV_BAP1mrna_validation.pdf`, `Plots/bulk_genetic_RUNX_epithelioid.pdf`.
Tables: `within_sample_TFactivity_BAP1_high_vs_low.csv`, `within_sample_TF_meta.csv`,
`scRNA_persample_clonality.csv`, `scRNA_CNV_BAP1mrna_validation.csv`,
`persample_chr3p_bimodality.csv`, `crossmodal_scRNA_vs_scATAC_BAP1.csv`,
`bulk_genetic_BAP1_histology.csv`, `bulk_genetic_RUNX_by_BAP1.csv`,
`bulk_genetic_epithelioid_DE_allgenes.csv`.

Bulk files: `R2_Q14_bulk_corroboration.R` (run via `submit_bulk_corroboration.sh`),
`bulk_RUNX_BAP1_correlation.csv`, `mesomics_IHC_BAP1_RUNX_grouptest.csv`,
`tcga_bueno_BAP1tertile_RUNX_grouptest.csv`; plots
`Plots/bulk_RUNX_vs_BAP1_scatter.pdf`, `Plots/mesomics_IHC_BAP1_RUNX_boxplots.pdf`,
`Plots/bulk_RUNX_BAP1_correlation_heatmap.pdf`.

## Files
- `R2_Q14_BAP1_TF_activity.R` — scATAC analysis (run via `submit_BAP1_TF_activity.sh`)
- `R2_Q14_bulk_corroboration.R` — bulk RNA-seq corroboration
- `probe.R`, `probe_bulk.R`, `probe_bulk2.R` — object inspection
- `BAP1_TF_activity_master_table.csv` — per-TF: group deltas + p/padj,
  Spearman/Pearson vs BAP1, per-cell supplement, RUNX flag
- `TFactivity_deviation_per_sample.csv`, `BAP1_genescore_per_sample.csv` — plot inputs
- `Plots/` — scATAC: volcano (group), correlation (vs BAP1), RUNX by status,
  RUNX vs BAP1 gene score, BAP1 gene score validation; bulk: RUNX-vs-BAP1
  scatter, MESOMICS IHC boxplots, correlation heatmap

## Caveats
- 11 samples (5 vs 6); small-n non-parametric tests. Low-cell tumors (P13=18,
  P3=21, P14=131) give noisier pseudobulk means.
- padj (BH over 869 TFs) leaves most individual TFs non-significant after
  correction; the RUNX signal is a consistent *trend* across two independent
  read-outs (group test + BAP1-gene-score correlation), not a single test.
