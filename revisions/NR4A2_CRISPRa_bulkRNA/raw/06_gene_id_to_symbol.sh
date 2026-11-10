#!/bin/bash
# 06_gene_id_to_symbol.sh: Ensembl gene id (with version) -> gene name table used by the DE step, from the GENCODE v46 GTF used for the STAR index
# (gene records only, file order). Reproduces processed/counts/gene_id_to_symbol.tsv byte for byte.
# Usage: ./06_gene_id_to_symbol.sh [GTF] [output]
GTF=${1:-/sc/arion/projects/Tsankov_Normal_Lung/genome_assemblies/human/gencode.v46.annotation.gtf}
OUT=${2:-/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CROPseq_validation_bulkRNA/processed/counts/gene_id_to_symbol.tsv}
awk -F'\t' '$3=="gene"{match($9,/gene_id "[^"]+"/); id=substr($9,RSTART+9,RLENGTH-10); match($9,/gene_name "[^"]+"/); nm=substr($9,RSTART+11,RLENGTH-12); print id"\t"nm}' "$GTF" > "$OUT"
