# R2_Q14 — TF activity vs BAP1 status (malignant compartment)

Reviewer 2, Q14: does TF activity differ between tumours with and without BAP1 loss? BAP1 status per sample (genetic):
lost P3, P4, P5, P8, P10; retained P1, P11, P12, P13, P14, P23.

Run in this order (LSF submit scripts alongside):

1. `R2_Q14_BAP1_TF_activity.R` (`submit_BAP1_TF_activity.sh`) — `tumor_compartment/scatac_ArchR`: per-sample mean chromVAR
   deviation and BAP1 gene score; lost vs retained (t-test / Wilcoxon, sample as unit), Spearman of activity with the BAP1
   gene score, per-cell presto test as a supplement. Tables `BAP1_TF_activity_master_table.csv`,
   `group_test_lost_vs_retained.csv`, `correlation_TFactivity_vs_BAP1genescore.csv`, `TFactivity_deviation_per_sample.csv`,
   `BAP1_genescore_per_sample.csv`. Figure `Plots/BAP1_genescore_by_status.pdf` (BAP1 gene score by genetic status).
2. `R2_Q14_histology_disentangle.R` (`submit_histology_disentangle.sh`) — BAP1 vs histology in the bulk cohorts
   (`bulkRNA_meso`), RUNX ~ BAP1 adjusted for the sarcomatoid score (`RUNX_BAP1_histology_adjusted.csv`), and the same in
   the scATAC samples (needs the two per-sample tables of step 1 and
   `tumor_compartment/scrna/cnmf20_sarcomatoid_sample_order.csv`). Figure `Plots/scATAC_BAP1_vs_sarcscore.pdf`.
