#!/usr/bin/env python3
"""
pileup_pattern24_nkt_celltypes.py  (R1_Q4)

Extracts genomic positions of pos_patterns.pattern_24 seqlets from the
KLRC1_NK MoDISco-counts H5, then computes and overlays the average
ChromBPNet contribution score profile at those sites for all 6 NKT cell
types (CD4, CD8, CD8_exhausted, FGFBP2_NK, KLRC1_NK, Tregs).

pattern_24 was identified as a composite AP-1 + NR4A2 motif by visual
inspection of the KLRC1_NK modisco report.

Coordinate mapping
------------------
ChromBPNet uses 2114-bp windows centred on peak summits.
MoDISco (`-w 400`) processes the central 400 bp (chr-BPNet positions
857-1256), so seqlet positions within the H5 are 0-based offsets into
that 400-bp sub-window:
    genomic_seqlet_center = (peak_start + summit_offset)  # summit
                            + seqlet_centre_in_400bp_win  # 0-399
                            - 200                          # centre of 400bp win

Output
------
  R1_Q4/output/NR4A2_analysis/pattern24_footprint_nkt_celltypes.pdf
"""

import os
import sys
import numpy as np
import h5py
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ──────────────────────────────────────────────────────────────────────
MESO_DIR   = "/sc/arion/scratch/giottb01/Xmen/meso"
MODISCO_H5 = f"{MESO_DIR}/KLRC1_NK/no_bias_model/modisco_counts/modisco_results_counts.h5"
PEAKS_BED  = f"{MESO_DIR}/KLRC1_NK/KLRC1_NK_peakset_all_no_blacklist.bed"
OUT_FILE   = ("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/"
              "scATAC_PM/git_repo_claude/R1_Q4/output/NR4A2_analysis/"
              "pattern24_footprint_nkt_celltypes.pdf")

CELLTYPES = {
    "CD4":          "CD4",
    "CD8":          "CD8",
    "CD8_exhausted":"CD8_exhausted",
    "FGFBP2_NK":   "FGFBP2_NK",
    "KLRC1_NK":    "KLRC1_NK",
    "Tregs":       "Tregs",
}

COLORS = {
    "CD4":           "#4DAF4A",
    "CD8":           "#377EB8",
    "CD8_exhausted": "#E41A1C",
    "FGFBP2_NK":    "#FF7F00",
    "KLRC1_NK":     "#984EA3",
    "Tregs":        "#A65628",
}

PATTERN   = "pos_patterns/pattern_24"
HALF_WIN  = 300   # bp each side of seqlet centre
MODISCO_HALF = 200  # half of the 400-bp modisco window

# ── 1. Load seqlet positions from modisco H5 ──────────────────────────────────
print(f"Loading seqlets for {PATTERN} from {MODISCO_H5}")
with h5py.File(MODISCO_H5, "r") as f:
    grp        = f[PATTERN]["seqlets"]
    ex_idx     = grp["example_idx"][()]   # peak row index
    sl_start   = grp["start"][()]         # start in 400-bp modisco window (0-based)
    sl_end     = grp["end"][()]           # end   in 400-bp modisco window
    is_rc      = grp["is_revcomp"][()]    # bool
print(f"  {len(ex_idx)} seqlets")

# ── 2. Load peak BED → summit genomic position ────────────────────────────────
print(f"Loading peaks from {PEAKS_BED}")
peaks = []
with open(PEAKS_BED) as fh:
    for line in fh:
        c = line.rstrip().split("\t")
        chrom        = c[0]
        peak_start   = int(c[1])       # 0-based BED start
        summit_off   = int(c[9])       # summit offset from peak_start
        summit_pos   = peak_start + summit_off  # 0-based genomic summit
        peaks.append((chrom, summit_pos))
peaks = np.array(peaks, dtype=object)
print(f"  {len(peaks)} peaks loaded")

# ── 3. Compute seqlet genomic centres ─────────────────────────────────────────
sl_centre_in_win = (sl_start + sl_end) / 2.0   # centre within 400-bp window
sl_offset_from_summit = sl_centre_in_win - MODISCO_HALF  # signed offset

seqlet_sites = []
for i, (eidx, off, rc) in enumerate(zip(ex_idx, sl_offset_from_summit, is_rc)):
    chrom, summit = peaks[eidx]
    centre = int(round(summit + off))
    seqlet_sites.append((chrom, centre, "-" if rc else "+"))

print(f"  Seqlet genomic centres computed (first 3): {seqlet_sites[:3]}")

# ── 4. Pileup at seqlet positions for each cell type ─────────────────────────
def pileup(bw_path, sites, half_win=HALF_WIN, max_sites=5000):
    if not os.path.isfile(bw_path):
        print(f"  [WARN] bigwig not found: {bw_path}", file=sys.stderr)
        return None, 0
    bw = pyBigWig.open(bw_path)
    mat = []
    if len(sites) > max_sites:
        np.random.seed(42)
        idx   = np.random.choice(len(sites), max_sites, replace=False)
        sites = [sites[i] for i in idx]
    for chrom, centre, strand in sites:
        s, e = centre - half_win, centre + half_win
        if s < 0:
            continue
        try:
            v = bw.values(chrom, s, e, numpy=True)
        except Exception:
            continue
        if v is None or len(v) != 2 * half_win:
            continue
        v = np.array(v, dtype=float)
        v[np.isnan(v)] = 0.0
        if strand == "-":
            v = v[::-1]
        mat.append(v)
    bw.close()
    if not mat:
        return None, 0
    return np.nanmean(mat, axis=0), len(mat)

x = np.arange(-HALF_WIN, HALF_WIN)

fig, ax = plt.subplots(figsize=(8, 4))

for ct_label, bw_prefix in CELLTYPES.items():
    bw_path = os.path.join(MESO_DIR, bw_prefix, "no_bias_model",
                           f"{bw_prefix}_averaged_contribution_scores_counts.bw")
    avg, n = pileup(bw_path, seqlet_sites)
    if avg is None:
        print(f"  [{ct_label}] skipped — no data")
        continue
    print(f"  [{ct_label}] n={n} sites used")
    color = COLORS[ct_label]
    ax.plot(x, avg, color=color, linewidth=1.4,
            label=f"{ct_label} (n={n})", alpha=0.9)

ax.axvline(0, color="black", linewidth=0.7, linestyle="--", alpha=0.4,
           label="seqlet centre")

# Shade the seqlet window (approx 30bp, the pattern width)
half_pat = 15
ax.axvspan(-half_pat, half_pat, alpha=0.07, color="grey")

ax.set_xlabel("Distance from pattern_24 seqlet centre (bp)", fontsize=10)
ax.set_ylabel("Avg. contribution score (counts model)", fontsize=10)
ax.set_title(
    "pos_patterns.pattern_24 (KLRC1_NK) — AP-1 + NR4A2 composite\n"
    "ChromBPNet contribution score pileup across NKT cell types",
    fontsize=10, fontweight="bold"
)
ax.legend(fontsize=8, frameon=False, loc="upper right")
ax.tick_params(labelsize=9)
ax.spines[["top", "right"]].set_visible(False)

fig.tight_layout()
fig.savefig(OUT_FILE, bbox_inches="tight", dpi=150)
print(f"\nSaved: {OUT_FILE}")
