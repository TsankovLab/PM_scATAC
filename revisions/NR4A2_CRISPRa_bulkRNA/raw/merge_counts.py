#!/usr/bin/env python3
"""Merge STAR ReadsPerGene.out.tab files into one gene x sample count matrix,
and report per-sample strandedness diagnostics (col2=unstranded, col3=fwd, col4=rev)."""
import glob
import os
import sys
import csv

COUNTS_DIR = "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CROPseq_validation_bulkRNA/processed/counts"
OUT_MATRIX_UNSTRANDED = os.path.join(COUNTS_DIR, "gene_counts_matrix_unstranded.tsv")
OUT_MATRIX_STRANDED = os.path.join(COUNTS_DIR, "gene_counts_matrix_stranded_reverse.tsv")
OUT_STRAND_QC = os.path.join(COUNTS_DIR, "strandedness_qc.tsv")

files = sorted(glob.glob(os.path.join(COUNTS_DIR, "*_ReadsPerGene.out.tab")))
if not files:
    sys.exit("No *_ReadsPerGene.out.tab files found in " + COUNTS_DIR)

samples = []
gene_ids = None
unstranded = {}
rev_stranded = {}
strand_totals = {}

for fpath in files:
    sid = os.path.basename(fpath).replace("_ReadsPerGene.out.tab", "")
    samples.append(sid)
    genes = []
    col2, col3, col4 = [], [], []
    with open(fpath) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0].startswith("N_"):
                continue  # skip N_unmapped / N_multimapping / N_noFeature / N_ambiguous summary rows
            genes.append(parts[0])
            col2.append(int(parts[1]))
            col3.append(int(parts[2]))
            col4.append(int(parts[3]))
    if gene_ids is None:
        gene_ids = genes
    elif gene_ids != genes:
        sys.exit(f"Gene order mismatch in {fpath}")
    unstranded[sid] = col2
    rev_stranded[sid] = col4
    strand_totals[sid] = (sum(col2), sum(col3), sum(col4))

with open(OUT_MATRIX_UNSTRANDED, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene_id"] + samples)
    for i, gid in enumerate(gene_ids):
        w.writerow([gid] + [unstranded[s][i] for s in samples])

with open(OUT_MATRIX_STRANDED, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene_id"] + samples)
    for i, gid in enumerate(gene_ids):
        w.writerow([gid] + [rev_stranded[s][i] for s in samples])

with open(OUT_STRAND_QC, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["Sample_ID", "unstranded_total", "fwd_stranded_total", "rev_stranded_total", "likely_protocol"])
    for sid in samples:
        u, fw, rv = strand_totals[sid]
        if fw > 0.8 * u and rv < 0.2 * u:
            proto = "stranded_fwd (use col3)"
        elif rv > 0.8 * u and fw < 0.2 * u:
            proto = "stranded_rev (use col4)"
        else:
            proto = "unstranded (use col2)"
        w.writerow([sid, u, fw, rv, proto])

print(f"Wrote unstranded count matrix: {OUT_MATRIX_UNSTRANDED} ({len(gene_ids)} genes x {len(samples)} samples)")
print(f"Wrote stranded-reverse count matrix: {OUT_MATRIX_STRANDED} ({len(gene_ids)} genes x {len(samples)} samples)")
print(f"Wrote strandedness QC: {OUT_STRAND_QC}")
