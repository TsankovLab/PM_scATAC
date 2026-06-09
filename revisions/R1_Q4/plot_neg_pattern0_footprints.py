#!/usr/bin/env python3
"""
plot_neg_pattern0_footprints.py  (R1_Q4)

Averages ChromBPNet footprint H5 files across 5 folds per cell type,
then overlays the neg_pattern0 footprint profile for all 6 NKT
cell types in a single plot.

Called by run_pattern24_footprints.sh after all fold jobs complete.
"""

import argparse
import glob
import os
import sys
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CELLTYPES = ["CD4", "CD8", "CD8_exhausted", "FGFBP2_NK", "KLRC1_NK", "Tregs"]
MOTIF_GROUP = "neg_pattern0"

COLORS = {
    "CD4":           "#4DAF4A",
    "CD8":           "#377EB8",
    "CD8_exhausted": "#E41A1C",
    "FGFBP2_NK":    "#FF7F00",
    "KLRC1_NK":     "#984EA3",
    "Tregs":        "#A65628",
}

OUT_PDF = ("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/"
           "scATAC_PM/git_repo_claude/R1_Q4/output/NR4A2_analysis/"
           "neg_pattern0_chrombpnet_footprints.pdf")


def average_folds(out_dir, celltype, motif_group):
    """Average i0 across fold H5 files for one cell type."""
    files = sorted(glob.glob(os.path.join(out_dir, f"{celltype}_fold*_footprints.h5")))
    if not files:
        print(f"  [{celltype}] no fold files found — skipping", file=sys.stderr)
        return None
    print(f"  [{celltype}] averaging {len(files)} folds: {[os.path.basename(f) for f in files]}")
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

    # Save averaged H5
    avg_path = os.path.join(out_dir, f"{celltype}_average_footprints.h5")
    with h5py.File(avg_path, "w") as f:
        grp = f.create_group(motif_group)
        grp.create_dataset("i0", data=avg)
        grp.create_dataset("i1", data=np.array([len(arrays)]))
    print(f"    saved → {avg_path}  shape={avg.shape}")
    return avg


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()

    fig, ax = plt.subplots(figsize=(8, 4))

    for ct in CELLTYPES:
        avg = average_folds(args.out_dir, ct, MOTIF_GROUP)
        if avg is None:
            continue

        # Show ±100bp around the centre of the 1000bp footprint window
        center = avg.shape[0] // 2
        window = avg[center - 100: center + 100]
        x = np.arange(-100, 100)

        ax.plot(x, window, color=COLORS[ct], linewidth=1.5,
                label=ct, alpha=0.9)

    ax.axvline(0, color="black", linewidth=0.7, linestyle="--", alpha=0.4)
    ax.axvspan(-6, 6, alpha=0.07, color="grey")   # approx motif half-width

    ax.set_xlabel("Distance from motif centre (bp)", fontsize=10)
    ax.set_ylabel("Predicted accessibility (probability)", fontsize=10)
    ax.set_title(
        "neg_patterns.pattern_0 (suppressive)\n"
        "ChromBPNet footprint across NKT cell types (KLRC1_NK peak regions)",
        fontsize=10, fontweight="bold"
    )
    ax.legend(fontsize=8, frameon=False)
    ax.tick_params(labelsize=9)
    ax.spines[["top", "right"]].set_visible(False)

    os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUT_PDF, bbox_inches="tight", dpi=150)
    print(f"\nSaved: {OUT_PDF}")


if __name__ == "__main__":
    main()
