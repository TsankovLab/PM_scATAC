#!/bin/bash
# 04_gene_id_to_symbol.sh: Ensembl gene id (with version) -> gene name table used by the DE step, from the GENCODE v46 GTF used for the STAR index
# (gene records only, file order). Reproduces 04_processed_counts/star_counts/gene_id_to_symbol.tsv byte for byte.
# Usage: ./04_gene_id_to_symbol.sh [GTF] [output]
source "${MESO_CONFIG:-$(for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$PWD"; do while [ ! -f "$d/config.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; [ -f "$d/config.sh" ] && { echo "$d/config.sh"; break; }; done)}"   # data / results / reference roots: see config.sh at the repo root
GTF=${1:-${GENOME_DIR}/human/gencode.v46.annotation.gtf}
OUT=${2:-${DATA_DIR}/CROPseq_validation_bulkRNA/04_processed_counts/star_counts/gene_id_to_symbol.tsv}
awk -F'\t' '$3=="gene"{match($9,/gene_id "[^"]+"/); id=substr($9,RSTART+9,RLENGTH-10); match($9,/gene_name "[^"]+"/); nm=substr($9,RSTART+11,RLENGTH-12); print id"\t"nm}' "$GTF" > "$OUT"
