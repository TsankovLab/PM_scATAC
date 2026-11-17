#!/bin/bash
# Reproduces the final bulk RNA-seq KO figures of spec195 (output folder $CROP_OUT, default: $RESULTS_DIR/CROP_seq_validation):
#   final6KO_all.pdf (TMM log-CPM; every page also as final6KO_p*.png/pdf); NORM=tpm Rscript spec195_final_6KO_figures.R gives the _TPM versions
# Samples: the bulk KO cohort defined in repro_common.R (COHORT: 35 KO samples of 6 KOs + 8 NTC = 43 samples, MSTO-211H, NCI-H2452, NCI-H28).
# Scripts, in dependency order:
#   spec118  reads at the gRNA target sites in the BAMs; spec119 guide reads per sample; spec121 indel and guide-read matrices
#   spec155a non-targeting guide reads; spec156 indels at the NTC guide matches
#   spec195  the figures
#   final6KO_normalise.R  the same expression processing as spec195, written as tables ($CROP_OUT/final6KO_normalised/: counts, TMM log-CPM,
#            cell-cycle regressed, module scores, KO profiles, screen KO means/profiles); its check lines reproduce spec195 page 1 C and D
# Helpers sourced by these scripts: repro_common.R, spec_common.R, profiles.R, agree_common.R.
# External inputs (not in this folder): STAR counts (PROC/star_counts: 04_processed_counts of the data folder), coordinate-sorted indexed STAR BAMs with logs (BAM_DIR),
#   the CROP-seq screen object (load_screen()), git_repo/files/scRNA_PM_modules.rds, "Oligos NR4A2+CROPseq.xlsx" (guide sequences; OLIGO_XLSX, default next to the scripts;
#   read by spec118, spec119), CCLE RNA-seq RSEM TPM, GENCODE v46 GTF, hg38 genome. spec118/spec119/spec155a/spec156 call samtools.
set -euo pipefail
module load R/4.2.0 samtools
cd "$(dirname "$0")"
source "${MESO_CONFIG:-$(for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$PWD"; do while [ ! -f "$d/config.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; [ -f "$d/config.sh" ] && { echo "$d/config.sh"; break; }; done)}"   # data / results / reference roots: see config.sh at the repo root
run() { echo "=== $* ==="; env "$@"; }
for s in spec118_gRNA_site_reads.R spec119_guide_reads_in_samples.R spec121_indel_heatmap.R spec155a_count_NTC_guides.R spec156_NTC_guides_genome_matches_and_indels.R \
         spec195_final_6KO_figures.R; do run Rscript "$s"; done
run Rscript final6KO_normalise.R
