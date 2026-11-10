# R1_Q5 — functional validation of the NR4A2 distal enhancer (CRISPRa, Jurkat)

A CRISPRa (dCas9-activator) assay targeting the element was run in Jurkat cells and NR4A2 mRNA measured by RT-qPCR
against two normalisers (beta-actin, 18S).

- `R1_Q5_NR4A2_CRISPRa_barplot.R` — recomputes dCt, ddCt and fold change from the raw Ct values of
  `../NR4A2_CRISPRa_Lillian/7_15_26_qPCR_JK_CRISPRa.xlsx` (read-only), for the analysis set of the spreadsheet's final
  columns (targeting NR4A2 2, 4, 5; scramble 1, 3, 6/24/26), checks them against the spreadsheet, and draws the barplot
  (mean fold change, SEM, replicate points; equal-variance t-test on dCt).

Output: `Plots/NR4A2_CRISPRa_barplot.pdf`, `NR4A2_CRISPRa_barplot_values.csv`.
