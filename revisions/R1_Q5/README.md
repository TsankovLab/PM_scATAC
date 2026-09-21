# R1_Q5 — Functional validation of the NR4A2 distal enhancer (CRISPRa, Jurkat)

Reviewer 1, Q5 asks for functional validation of the NR4A2 distal enhancer. A CRISPRa
(dCas9-activator) assay targeting the element was run in Jurkat cells and NR4A2 mRNA measured
by RT-qPCR against two normalisers.

## Source data

`../NR4A2_CRISPRa_Lillian/7_15_26_qPCR_JK_CRISPRa.xlsx` (read-only; never modified)

| | |
|---|---|
| NR4A2-targeting gRNA | 5 biological replicates, 7/15/26 |
| Scramble gRNA | 3 biological replicates, 7/15/26 + 1 from 6/24/26 |
| Genes | NR4A2, beta-actin, 18S — technical triplicates |

## Analysis

`R1_Q5_NR4A2_CRISPRa_qPCR.R` recomputes everything from the raw Ct values rather than reusing
the spreadsheet's derived columns. Standard delta-delta-Ct, with the test performed on dCt
(already log2 units) rather than on fold changes, which are log-normal:

    dCt  = mean Ct(NR4A2) - mean Ct(reference)      per biological replicate
    ddCt = dCt - mean dCt(scramble)
    fold = 2^(-ddCt)                                Welch t-test on dCt

Three normalisers are reported: beta-actin, 18S, and their mean Ct (geometric mean of the two).
Technical-triplicate means reproduce the spreadsheet's own averaged-Ct block exactly
(max difference 3e-5 Ct), so the parsing is verified.

## Result — the direction is right, the significance depends on which replicates are kept

Fold change (geometric mean normaliser), Welch t-test on dCt:

| replicates used | n | fold change (95% CI) | p |
|---|---|---|---|
| all | 5 vs 3 | 2.2 (0.2–21.2) | 0.44 |
| all + 6/24 scramble | 5 vs 4 | 2.7 (0.4–19.8) | 0.27 |
| drop NR4A2 1 (technical QC) | 4 vs 4 | 4.7 (0.8–27.1) | 0.077 |
| spreadsheet subset | 3 vs 2 | 10.7 (2.7–43.0) | 0.015 |
| spreadsheet subset + 6/24 | 3 vs 3 | 11.9 (3.0–47.4) | 0.013 |

Both normalisers agree throughout, so the conclusion does not depend on the reference gene.

### Per-replicate fold change (vs scramble mean)

| replicate | fold | dropped in source file |
|---|---|---|
| NR4A2 1 | 0.32 | yes |
| NR4A2 2 | 6.57 | |
| NR4A2 3 | 1.21 | yes |
| NR4A2 4 | 14.34 | |
| NR4A2 5 | 4.10 | |
| Scramble 1 | 0.78 | |
| Scramble 2 | 4.39 | yes |
| Scramble 3 | 0.59 | |
| Scramble 6/24 | 0.49 | |

**3 of 5 targeting replicates show clear induction (4–14x); 2 show none.** Scramble is tightly
clustered (0.49–0.78) except Scramble 2 (4.39).

### The exclusions need a stated justification

The source file's significant result drops exactly the three replicates that most oppose the
hypothesis: the two non-responding targeting replicates and the one high scramble. p moves from
0.44 to 0.015 on that basis alone.

Of the three, only **NR4A2 1** has support in the assay itself: its technical triplicates
disagree badly (SD 0.70 Ct for NR4A2, ~3x every other well). **NR4A2 3** (SD 0.061) and
**Scramble 2** (SD 0.177) replicate cleanly and have no technical grounds for removal in the
data provided. Dropping only the QC failure gives 4.7-fold, p = 0.077.

If there is an independent reason to exclude NR4A2 3 and Scramble 2 — transduction efficiency,
selection, cell viability, dCas9 or gRNA expression — it should be recorded and applied as a
pre-specified rule; the analysis then stands. Without one, the defensible statement is a
consistent direction of effect that does not reach significance at n = 5 vs 4.

A rank test is uninformative at this n (minimum attainable p with 5 vs 4 is 0.016) and is
reported only for completeness.

## Outputs

| file | contents |
|---|---|
| `NR4A2_CRISPRa_replicate_level.csv` | per-replicate Ct, dCt, fold change |
| `NR4A2_CRISPRa_summary_stats.csv` | every variant x normaliser, fold change, CI, p |
| `NR4A2_CRISPRa_technical_QC.csv` | technical triplicate SD per well |
| `Plots/NR4A2_CRISPRa_qPCR.pdf` | dCt and fold change by group and normaliser |
| `Plots/NR4A2_CRISPRa_sensitivity.pdf` | how the result moves with replicate inclusion |

## Open items

- Enhancer coordinates and gRNA sequences are not in the spreadsheet; needed for the methods
  (the element analysed elsewhere in the revision is chr2:156,480,366-156,480,866).
- No transduction-efficiency or dCas9-expression readout accompanies the qPCR.
- The 6/24/26 scramble has no matching targeting sample from that date.
- n = 5 vs 4 is low for an assay with this spread; 2-3 further targeting replicates would
  settle it.
