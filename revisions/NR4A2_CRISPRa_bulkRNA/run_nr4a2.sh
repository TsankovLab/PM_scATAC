#!/bin/bash
# NR4A2 enhancer CRISPRa vs SCRAMBLE (Jurkat) bulk RNA-seq: from the sequencing run to 05_volcano_top_pathway_and_dotplot.{png,pdf} and the GEO matrices.
# 1. raw data -> per-sample gene counts (LSF jobs; paths are inside the scripts):
#      tar -xf <run>.BCL.tar -C <scratch>              BCL run folder (RUN_DIR of raw/02)
#      bsub < raw/01_star_index.lsf                    STAR 2.7.11b index, hg38 + GENCODE v46
#      bsub < raw/02_bclconvert_demux.lsf              bcl-convert 4.4.6 -> FASTQ
#      bsub < raw/03_star_align.lsf                    FastQC + STAR alignment, --quantMode GeneCounts -> <sample>_ReadsPerGene.out.tab
#      raw/04_gene_id_to_symbol.sh                     gene_id_to_symbol.tsv from the GTF
# 2. counts -> results (this script): count matrix, DE (R 4.4.1, DESeq2 1.46.0), GSEA (system Python 3.9, gseapy 1.3.1), figure and GEO matrices (R 4.4.1).
#    Outputs in DE_DIR (default processed/de/NR4A2_1_2_vs_SCRAMBLE).
set -euo pipefail
cd "$(dirname "$0")"
source "${MESO_CONFIG:-$(for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$PWD"; do while [ ! -f "$d/config.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; [ -f "$d/config.sh" ] && { echo "$d/config.sh"; break; }; done)}"   # data / results / reference roots: see config.sh at the repo root
module load R/4.4.1
Rscript 00_nr4a2_count_matrix.R
Rscript 01_nr4a2_deseq2.R
PYTHONPATH= /usr/bin/python3 02_nr4a2_gsea.py   # system Python 3.9 + gseapy 1.3.1 (pip --user); module R/4.4.1 puts Python 3.12 / gseapy 1.1.3 first,
                                                 # which gives different permutation results
Rscript nr4a2_volcano_pathway_dotplot.R
Rscript 03_nr4a2_geo_matrices.R   # GEO processed data: raw and DESeq2-normalised counts of the 5 samples (DE_DIR/GEO)
