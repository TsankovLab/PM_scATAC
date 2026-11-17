#!/bin/bash
# Bulk RNA-seq KO cohort, 35 transcription-factor knockout samples (6 KOs x MSTO-211H, NCI-H2452, NCI-H28) + 8 non-targeting controls: from the sequencing run to the processed count matrices.
# 1. raw data -> per-sample gene counts (LSF jobs; paths are inside the scripts):
#      tar -xf <run>.BCL.tar -C <scratch>              BCL run folder (RUN_DIR of raw/02)
#      bsub < raw/01_star_index.lsf                    STAR 2.7.11b index, hg38 + GENCODE v46
#      bsub < raw/02_bclconvert_demux.lsf              bcl-convert 4.4.6 -> FASTQ
#      bsub < raw/03_star_align.lsf                    FastQC + STAR alignment, --quantMode GeneCounts -> <sample>_ReadsPerGene.out.tab
#      raw/04_gene_id_to_symbol.sh                     gene_id_to_symbol.tsv from the GTF
# 2. counts -> processed matrices (this script): raw counts, TMM log2 CPM and sample annotation (OUT_DIR, default 04_processed_counts/KO35)
set -euo pipefail
cd "$(dirname "$0")"
source "${MESO_CONFIG:-$(for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$PWD"; do while [ ! -f "$d/config.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; [ -f "$d/config.sh" ] && { echo "$d/config.sh"; break; }; done)}"   # data / results / reference roots: see config.sh at the repo root
module load R/4.4.1
Rscript 00_ko35_count_matrix.R
