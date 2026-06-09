"""
Parse SCENIC motifs.csv to extract TF → target gene table.
Also identifies:
  - Whether SOX9/RUNX2/SNAI2 appear as targets in any regulon
  - The full SOX9 regulon target list
  - Overlap between SOX9 targets and RUNX2/SNAI2 gene programs
"""

import csv
import ast
import re
import os
import sys
from collections import defaultdict

scenic_dir = (
    '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
    '/tumor_compartment/scrna/SCENIC/vg_5000_mw_tss500bp/tumor_programs'
)
outdir = (
    '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
    '/git_repo_claude/R1_Q3'
)
motifs_csv = os.path.join(scenic_dir, 'motifs.csv')
auc_csv    = os.path.join(scenic_dir, 'auc_mtx.csv')

anchor_genes = {'SOX9', 'RUNX2', 'SNAI2'}

# ---------------------------------------------------------------------------
# Parse motifs.csv
# ---------------------------------------------------------------------------
# Format: rows 0-1 are headers; data starts at row 2
# Columns: TF, MotifID, AUC, NES, MotifSimilarityQvalue, OrthologousIdentity,
#          Annotation, Context, TargetGenes, RankAtMax

def parse_target_genes(raw):
    """Extract gene set from a string like frozenset({'GENE1', 'GENE2', ...})"""
    # strip frozenset wrapper
    inner = re.sub(r"^frozenset\(\{", '', raw.strip())
    inner = re.sub(r"\}\)$", '', inner)
    genes = set()
    for token in inner.split(','):
        token = token.strip().strip("'\"")
        if token and not token.startswith('('):
            genes.add(token)
    return genes

tf_regulons   = defaultdict(lambda: {'auc': 0, 'nes': 0, 'targets': set(),
                                      'motif': '', 'n_motifs': 0})
all_targets_by_tf = defaultdict(set)   # union of all motif target sets per TF

with open(motifs_csv, 'r') as f:
    lines = f.readlines()

# skip first 2 header rows, then parse
for line in lines[2:]:
    # Use csv reader on each line to handle quoted fields
    try:
        row = next(csv.reader([line]))
    except StopIteration:
        continue
    if len(row) < 9:
        continue
    tf       = row[0].strip()
    motif    = row[1].strip()
    try:
        auc  = float(row[2])
        nes  = float(row[3])
    except ValueError:
        continue
    raw_tg   = row[8] if len(row) > 8 else ''
    targets  = parse_target_genes(raw_tg) if raw_tg else set()

    all_targets_by_tf[tf] |= targets
    tf_regulons[tf]['n_motifs'] += 1

    # Keep the highest-AUC motif as the representative regulon
    if auc > tf_regulons[tf]['auc']:
        tf_regulons[tf].update({'auc': auc, 'nes': nes,
                                 'motif': motif, 'targets': targets})

# ---------------------------------------------------------------------------
# Write full TF × target gene table (long format)
# ---------------------------------------------------------------------------
long_table_path = os.path.join(outdir, 'SCENIC_TF_target_genes_long.csv')
with open(long_table_path, 'w') as f:
    f.write('TF,target_gene,regulon_AUC,regulon_NES,best_motif,n_motifs\n')
    for tf in sorted(tf_regulons.keys()):
        reg = tf_regulons[tf]
        for gene in sorted(reg['targets']):
            f.write(f"{tf},{gene},{reg['auc']:.4f},{reg['nes']:.4f},"
                    f"{reg['motif']},{reg['n_motifs']}\n")

print(f"Long table written: {long_table_path}")

# Summary table: one row per TF
summary_path = os.path.join(outdir, 'SCENIC_TF_regulon_summary.csv')
with open(summary_path, 'w') as f:
    f.write('TF,n_target_genes,regulon_AUC,regulon_NES,best_motif,n_motifs,'
            'SOX9_is_target,RUNX2_is_target,SNAI2_is_target\n')
    for tf in sorted(tf_regulons.keys()):
        reg     = tf_regulons[tf]
        targets = all_targets_by_tf[tf]
        f.write(f"{tf},{len(targets)},{reg['auc']:.4f},{reg['nes']:.4f},"
                f"{reg['motif']},{reg['n_motifs']},"
                f"{'TRUE' if 'SOX9'  in targets else 'FALSE'},"
                f"{'TRUE' if 'RUNX2' in targets else 'FALSE'},"
                f"{'TRUE' if 'SNAI2' in targets else 'FALSE'}\n")

print(f"Summary table written: {summary_path}")

# ---------------------------------------------------------------------------
# Focused report: which TFs have SOX9 / RUNX2 / SNAI2 as targets?
# ---------------------------------------------------------------------------
print('\n=== TFs that REGULATE (target) SOX9 / RUNX2 / SNAI2 ===')
for anchor in sorted(anchor_genes):
    regulators = [tf for tf, targets in all_targets_by_tf.items()
                  if anchor in targets]
    print(f"\n{anchor} is a target of: {', '.join(sorted(regulators)) or 'none'}")

# SOX9 regulon details
print('\n=== SOX9 regulon targets ===')
sox9_targets = sorted(all_targets_by_tf.get('SOX9', set()))
print(f"N targets: {len(sox9_targets)}")
print(f"Targets: {', '.join(sox9_targets)}")
print(f"RUNX2 in SOX9 targets: {'RUNX2' in sox9_targets}")
print(f"SNAI2 in SOX9 targets: {'SNAI2' in sox9_targets}")

# Check which anchor genes are co-targets (appear in same regulon)
print('\n=== Co-targeting: TFs that regulate ≥2 anchor genes ===')
for tf, targets in all_targets_by_tf.items():
    shared = anchor_genes & targets
    if len(shared) >= 2:
        print(f"  {tf} → {', '.join(sorted(shared))}")

# ---------------------------------------------------------------------------
# Also extract AUC matrix TF names (active regulons in cells)
# ---------------------------------------------------------------------------
with open(auc_csv, 'r') as f:
    header = f.readline().strip().split(',')
regulon_names = [h.replace('(+)', '').replace('(-)', '').strip()
                 for h in header[1:]]   # skip 'Cell' column
print(f'\nN active regulons in AUC matrix: {len(regulon_names)}')
print('Active regulons containing SOX9/RUNX2/SNAI2 family:')
for r in regulon_names:
    if any(a in r for a in anchor_genes):
        print(f'  {r}')

print('\nDone.')
