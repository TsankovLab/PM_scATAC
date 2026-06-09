#!/usr/bin/env python3
"""
pileup_SOX9_contribution_scores.py  (R1_Q4)

Average ChromBPNet contribution score pileup at SOX9 and TEAD motif sites
within SOX9_regulon_high peaks.

If TEAD sites show strong contribution signal but SOX9 sites show
TEAD-like (not SOX9-like) patterns, that supports the YAP/TAZ-TEAD
co-binding hypothesis explaining the absent SOX9 MoDISco motif.

Inputs:
  SOX9_sites_SOX9_regulon_high.bed   exported by analyze_SOX9_cooccurrence.R
  TEAD_sites_SOX9_regulon_high.bed   exported by analyze_SOX9_cooccurrence.R
  SOX9_regulon_high_averaged_contribution_scores_counts.bw

Output:
  R1_Q4/output/SOX9_analysis/SOX9_contribution_pileup.pdf
"""

import os
import sys
import numpy as np
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

HALF_WIN = 300
BW_PATH  = ("/sc/arion/scratch/giottb01/SOX9_chromBPnet/"
            "SOX9_regulon_high/no_bias_model/"
            "SOX9_regulon_high_averaged_contribution_scores_counts.bw")
BED_DIR  = ("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/"
            "scATAC_PM/git_repo_claude/R1_Q4/output/SOX9_analysis")
OUT_FILE = os.path.join(BED_DIR, "SOX9_contribution_pileup.pdf")

BEDS = {
    "SOX9 motif sites":         os.path.join(BED_DIR, "SOX9_sites_SOX9_regulon_high.bed"),
    "TEAD motif sites\n(+ctrl)": os.path.join(BED_DIR, "TEAD_sites_SOX9_regulon_high.bed"),
}
COLORS = {"SOX9 motif sites": "#D73027", "TEAD motif sites\n(+ctrl)": "#4575B4"}

# ── helpers ────────────────────────────────────────────────────────────────────

def load_bed(path):
    sites = []
    if not os.path.exists(path):
        print(f"[WARN] not found: {path}")
        return sites
    with open(path) as f:
        for line in f:
            c = line.rstrip().split("\t")
            if len(c) < 6: continue
            center = (int(c[1]) + int(c[2])) // 2
            sites.append((c[0], center, c[5]))
    return sites

def pileup(bw_path, sites, half_win=HALF_WIN, max_sites=5000):
    bw  = pyBigWig.open(bw_path)
    mat = []
    if len(sites) > max_sites:
        np.random.seed(42)
        idx   = np.random.choice(len(sites), max_sites, replace=False)
        sites = [sites[i] for i in idx]
    for chrom, center, strand in sites:
        s, e = center - half_win, center + half_win
        if s < 0: continue
        try:
            v = bw.values(chrom, s, e, numpy=True)
        except Exception:
            continue
        if v is None or len(v) != 2*half_win: continue
        v = np.array(v, dtype=float); v[np.isnan(v)] = 0.0
        if strand == "-": v = v[::-1]
        mat.append(v)
    bw.close()
    if not mat: return None, 0
    return np.nanmean(mat, axis=0), len(mat)

# ── main ───────────────────────────────────────────────────────────────────────

x = np.arange(-HALF_WIN, HALF_WIN)

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=False)

for ax, (label, bed_path) in zip(axes, BEDS.items()):
    sites = load_bed(bed_path)
    print(f"{label.replace(chr(10),' ')}: {len(sites)} sites", flush=True)
    avg, n = pileup(BW_PATH, sites)
    if avg is None:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
    else:
        color = COLORS[label]
        ax.fill_between(x, avg, alpha=0.25, color=color)
        ax.plot(x, avg, color=color, linewidth=1.2,
                label=f"n={n} sites")
        ax.axvline(0, color="black", linewidth=0.6, linestyle="--", alpha=0.5)
        ax.legend(fontsize=8, frameon=False)

    ax.set_xlabel("Distance from motif centre (bp)", fontsize=9)
    ax.set_ylabel("Avg. contribution score (counts model)", fontsize=9)
    ax.set_title(label, fontsize=10, fontweight="bold")
    ax.tick_params(labelsize=8)
    ax.spines[["top","right"]].set_visible(False)

fig.suptitle("ChromBPNet contribution score pileup — SOX9_regulon_high model\n"
             "SOX9 motif sites vs TEAD sites (positive control)",
             fontsize=11, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(OUT_FILE, bbox_inches="tight")
print(f"\nSaved: {OUT_FILE}")
