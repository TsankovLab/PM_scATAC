#!/usr/bin/env python3
"""
plot_footprint_generic.py  (R1_Q4)
Averages ChromBPNet fold footprint H5s for one pattern and overlays all
6 NKT cell types on one plot. Called by run_footprint_generic.sh.
"""

import argparse, glob, os, sys
import h5py, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

CELLTYPES = ["CD4", "CD8", "CD8_exhausted", "FGFBP2_NK", "KLRC1_NK", "Tregs"]
COLORS    = {"CD4":"#4DAF4A","CD8":"#377EB8","CD8_exhausted":"#E41A1C",
             "FGFBP2_NK":"#FF7F00","KLRC1_NK":"#984EA3","Tregs":"#A65628"}
OUT_BASE  = ("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/"
             "scATAC_PM/git_repo_claude/R1_Q4/output/NKT_footprints")

def average_folds(out_dir, celltype, motif_group):
    files = sorted(glob.glob(os.path.join(out_dir, f"{celltype}_fold*_footprints.h5")))
    if not files:
        return None
    print(f"  [{celltype}] {len(files)} folds")
    arrays = []
    for fn in files:
        with h5py.File(fn, "r") as f:
            if motif_group not in f:
                print(f"    [WARN] group '{motif_group}' not in {fn}", file=sys.stderr)
                continue
            arrays.append(f[motif_group]["i0"][()])
    if not arrays:
        return None
    avg = np.mean(arrays, axis=0)
    avg_path = os.path.join(out_dir, f"{celltype}_average_footprints.h5")
    with h5py.File(avg_path, "w") as f:
        g = f.create_group(motif_group)
        g.create_dataset("i0", data=avg)
        g.create_dataset("i1", data=np.array([len(arrays)]))
    return avg

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out_dir",     required=True)
    p.add_argument("--pattern_id",  required=True)
    p.add_argument("--motif_group", required=True)
    args = p.parse_args()

    fig, ax = plt.subplots(figsize=(8, 4))
    for ct in CELLTYPES:
        avg = average_folds(args.out_dir, ct, args.motif_group)
        if avg is None:
            continue
        center  = avg.shape[0] // 2
        window  = avg[center - 100: center + 100]
        ax.plot(np.arange(-100, 100), window, color=COLORS[ct],
                linewidth=1.5, label=ct, alpha=0.9)

    ax.axvline(0, color="black", linewidth=0.7, linestyle="--", alpha=0.4)
    ax.axvspan(-4, 4, alpha=0.07, color="grey")
    ax.set_xlabel("Distance from motif centre (bp)", fontsize=10)
    ax.set_ylabel("Predicted accessibility (probability)", fontsize=10)
    pattern_display = args.pattern_id.replace("__", ".")
    ax.set_title(f"{pattern_display}\nChromBPNet footprint — NKT cell types (KLRC1_NK peaks)",
                 fontsize=10, fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    ax.tick_params(labelsize=9)
    ax.spines[["top", "right"]].set_visible(False)

    os.makedirs(OUT_BASE, exist_ok=True)
    out_pdf = os.path.join(OUT_BASE, f"{args.pattern_id}_footprints.pdf")
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight", dpi=150)
    print(f"Saved: {out_pdf}")

if __name__ == "__main__":
    main()
