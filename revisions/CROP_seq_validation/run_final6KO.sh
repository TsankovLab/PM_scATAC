#!/bin/bash
# Reproduces the final bulk RNA-seq KO figures of spec195 (output folder $CROP_OUT, default: <perturbation_panel>/figs_repro):
#   final6KO_all.pdf (TMM log-CPM; every page also as final6KO_p*.png/pdf), including final6KO_p4_sample_identity and final6KO_p6_CNA_vs_CCLE
#   (these two pages do not use the normalised values: identical to the final6KO_TPM_p4/p6 pages of the earlier NORM=tpm run)
# Only the scripts whose outputs spec195 reads, directly or through another script, in dependency order:
#   spec71   sample audit (spec71_sample_audit.csv, read by load_audit() in agree_common.R)
#   spec72   bulk profiles of the S70 sample set (spec72_profiles_S70.rds: sample list)
#   spec75   screen module scores per cell (spec75_screen_scores_by_cell.rds, needed by spec128)
#   spec118  reads at the gRNA target sites in the BAMs; spec119 guide reads per sample; spec121 indel and guide-read matrices
#   spec155a non-targeting guide reads; spec156 indels at the NTC guide matches
#   spec128  per-sample annotation SA (genotype line, corrected KO label, indel %) in spec128_inputs.rds
#   spec195  the figures
#   final6KO_normalise.R  the same expression processing as spec195, written as tables ($CROP_OUT/final6KO_normalised/: counts, TMM log-CPM,
#            cell-cycle regressed, module scores, KO profiles, screen KO means/profiles); its check lines reproduce spec195 page 1 C and D
# Helpers sourced by these scripts: repro_common.R, spec_common.R, profiles.R, agree_common.R.
# External inputs (not in this folder): STAR counts and logs, coordinate-sorted indexed STAR BAMs (BAM_DIR), the CROP-seq screen object (load_screen()),
#   git_repo/files/scRNA_PM_modules.rds, "Oligos NR4A2+CROPseq.xlsx" (guide sequences, next to the scripts; read by spec118, spec119, spec121, spec155a),
#   CCLE RNA-seq RSEM TPM, GENCODE v46 GTF, hg38 genome. spec118/spec119/spec155a/spec156 call samtools.
# The genotype-based cell-line corrections are registered in repro_common.R (they come from the SNP check of spec71b_*, kept in git_repo_claude).
set -euo pipefail
module load R/4.2.0 samtools
cd "$(dirname "$0")"
run() { echo "=== $* ==="; env "$@"; }
run Rscript spec71_sample_audit_76.R
run SETS=S70 Rscript spec72_build_profiles.R
run SET=S70 Rscript spec75_Cm_preserved.R
for s in spec118_gRNA_site_reads.R spec119_guide_reads_in_samples.R spec121_indel_heatmap.R spec155a_count_NTC_guides.R spec156_NTC_guides_genome_matches_and_indels.R \
         spec128_final_grid_build.R spec195_final_6KO_figures.R; do run Rscript "$s"; done
run Rscript final6KO_normalise.R
