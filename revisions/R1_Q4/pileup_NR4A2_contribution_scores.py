#!/usr/bin/env python3
"""
pileup_NR4A2_contribution_scores.py  (R1_Q4)

Computes average ChromBPNet contribution score pileup at NR4A2 and FOSL2 (AP-1)
motif sites for CD8_exhausted and KLRC1_NK.

If NR4A2 sites show AP-1-like contribution patterns but no NR4A2 footprint,
that supports the hypothesis that AP-1 co-binding dominates ChromBPNet signals
at those loci, explaining the absent NR4A2 MoDISco pattern.

Inputs:
  NR4A2_sites_<celltype>.bed    exported by analyze_NR4A2_cooccurrence.R
  FOSL2_sites_<celltype>.bed    exported by analyze_NR4A2_cooccurrence.R
  <celltype>_averaged_contribution_scores_counts.bw

Output:
  NR4A2_contribution_pileup.pdf   metaplot: avg contribution ±300bp around motif
"""

import os
import sys
import numpy as np
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

HALF_WIN  = 300   # bp each side of motif centre
BW_DIR    = "/sc/arion/scratch/giottb01/Xmen/meso"
BED_DIR   = "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/NR4A2_analysis"
OUT_FILE  = os.path.join(BED_DIR, "NR4A2_contribution_pileup.pdf")

CELLTYPES = {
    "CD8_exhausted": "CD8_exhausted",
    "NK_KLRC1":      "KLRC1_NK",       # bw filename uses KLRC1_NK
}

MOTIFS = ["NR4A2", "FOSL2"]

COLORS = {
    "CD8_exhausted_NR4A2": "#D73027",
    "CD8_exhausted_FOSL2": "#FC8D59",
    "NK_KLRC1_NR4A2":      "#4575B4",
    "NK_KLRC1_FOSL2":      "#74ADD1",
}

# ── helpers ────────────────────────────────────────────────────────────────────

def load_bed(path):
    """Load BED6 → list of (chr, center, strand)."""
    sites = []
    if not os.path.exists(path):
        print(f"  [WARN] BED not found: {path}")
        return sites
    with open(path) as f:
        for line in f:
            cols = line.rstrip().split("\t")
            if len(cols) < 6:
                continue
            chrom  = cols[0]
            start  = int(cols[1])
            end    = int(cols[2])
            strand = cols[5]
            center = (start + end) // 2
            sites.append((chrom, center, strand))
    return sites

def fetch_window(bw, chrom, center, half_win):
    """Extract contribution scores ±half_win around centre. Returns array or None."""
    s = center - half_win
    e = center + half_win
    if s < 0:
        return None
    try:
        vals = bw.values(chrom, s, e, numpy=True)
    except Exception:
        return None
    if vals is None or len(vals) != 2 * half_win:
        return None
    vals = np.array(vals, dtype=float)
    vals[np.isnan(vals)] = 0.0
    return vals

def pileup(bw_path, sites, half_win=HALF_WIN, max_sites=5000):
    """Average contribution scores over all sites. Flips RC strands."""
    if not os.path.exists(bw_path):
        print(f"  [WARN] bw not found: {bw_path}")
        return None, 0
    bw  = pyBigWig.open(bw_path)
    mat = []
    # sample sites if too many
    if len(sites) > max_sites:
        np.random.seed(42)
        idx   = np.random.choice(len(sites), max_sites, replace=False)
        sites = [sites[i] for i in idx]
    for chrom, center, strand in sites:
        arr = fetch_window(bw, chrom, center, half_win)
        if arr is None:
            continue
        if strand == "-":
            arr = arr[::-1]   # reverse complement: flip the window
        mat.append(arr)
    bw.close()
    if not mat:
        return None, 0
    return np.nanmean(mat, axis=0), len(mat)

# ── main ───────────────────────────────────────────────────────────────────────

x_axis = np.arange(-HALF_WIN, HALF_WIN)

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
ax_nr4a2, ax_fosl2 = axes

for ct_label, bw_prefix in CELLTYPES.items():
    bw_path = os.path.join(BW_DIR, bw_prefix, "no_bias_model",
                           f"{bw_prefix}_averaged_contribution_scores_counts.bw")

    for motif in MOTIFS:
        bed_path = os.path.join(BED_DIR, f"{motif}_sites_{ct_label}.bed")
        sites    = load_bed(bed_path)
        print(f"{ct_label} | {motif}: {len(sites)} sites", flush=True)

        avg, n_used = pileup(bw_path, sites)
        if avg is None:
            print(f"  -> skipped (no data)")
            continue

        label = f"{ct_label} (n={n_used})"
        color = COLORS.get(f"{ct_label}_{motif}", "grey")
        ax = ax_nr4a2 if motif == "NR4A2" else ax_fosl2
        ax.plot(x_axis, avg, color=color, linewidth=1.0, label=label, alpha=0.9)

for ax, title in zip(axes, ["NR4A2 motif sites", "FOSL2 (AP-1) sites — positive control"]):
    ax.axvline(0, color="black", linewidth=0.6, linestyle="--", alpha=0.5)
    ax.set_xlabel("Distance from motif centre (bp)", fontsize=9)
    ax.set_ylabel("Avg. contribution score", fontsize=9)
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.legend(fontsize=7, frameon=False)
    ax.tick_params(labelsize=8)
    ax.spines[["top","right"]].set_visible(False)

fig.suptitle("ChromBPNet contribution score pileup at NR4A2 vs AP-1 sites\n"
             "(CD8_exhausted & NK_KLRC1 — counts model)",
             fontsize=11, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(OUT_FILE, bbox_inches="tight")
print(f"\nSaved: {OUT_FILE}")
