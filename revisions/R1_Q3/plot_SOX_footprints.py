#!/usr/bin/env python3
"""
plot_SOX_footprints.py  (R1_Q3)
Averages fold H5s per cell type and plots SOX/HMG motif footprints
comparing SOX9_high_P23 vs SOX9_low_P23.
"""

import argparse, glob, os
import h5py, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CELLTYPES = ["SOX9_high_P23", "SOX9_low_P23"]
COLORS    = {"SOX9_high_P23": "#E41A1C", "SOX9_low_P23": "#377EB8"}

MOTIF_PAIRS = [
    ("SOX8_HMG_p48",  "SOX8_HMG_p48_rc",  "SOX8_HMG (de novo high, pattern_48)"),
    ("Sox1_HMG_p59",  "Sox1_HMG_p59_rc",  "Sox1_HMG (de novo low, pattern_59)"),
    ("SOX9_MA0077",   "SOX9_MA0077_rc",   "SOX9 canonical (JASPAR MA0077.1)"),
]

OUT_DIR_DEFAULT = "/sc/arion/scratch/giottb01/meso_SOX9_P23_chromBPnet/footprints_SOX_motifs"
OUT_PDF = ("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/"
           "scATAC_PM/git_repo_claude/R1_Q3/Plots/SOX_motifs_footprints.pdf")


def average_folds(out_dir, celltype):
    files = sorted(glob.glob(os.path.join(out_dir, f"{celltype}_fold*_footprints.h5")))
    if not files:
        return {}
    with h5py.File(files[0], "r") as f:
        groups = list(f.keys())
    result = {}
    for grp in groups:
        arrays = []
        for fn in files:
            with h5py.File(fn, "r") as f:
                if grp in f and "i0" in f[grp]:
                    arrays.append(f[grp]["i0"][()])
        if arrays:
            result[grp] = np.mean(arrays, axis=0)
    avg_path = os.path.join(out_dir, f"{celltype}_average_footprints.h5")
    with h5py.File(avg_path, "w") as f:
        for grp, arr in result.items():
            g = f.create_group(grp)
            g.create_dataset("i0", data=arr)
    print(f"  [{celltype}] averaged {len(files)} folds")
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out_dir", default=OUT_DIR_DEFAULT)
    args = p.parse_args()

    data = {ct: average_folds(args.out_dir, ct) for ct in CELLTYPES}

    n = len(MOTIF_PAIRS)
    fig, axes = plt.subplots(n, 1, figsize=(8, n * 3.2), sharex=True)
    if n == 1: axes = [axes]

    x = np.arange(-100, 100)

    for ax, (fwd_key, rc_key, title) in zip(axes, MOTIF_PAIRS):
        for ct in CELLTYPES:
            if not data[ct]: continue
            color = COLORS[ct]
            center = None
            if fwd_key in data[ct]:
                fp = data[ct][fwd_key]
                center = fp.shape[0] // 2
                ax.plot(x, fp[center-100:center+100],
                        color=color, linewidth=1.6, linestyle="-", alpha=0.9,
                        label=ct if fwd_key == MOTIF_PAIRS[0][0] else "")
            if rc_key in data[ct] and center is not None:
                fp_rc = data[ct][rc_key]
                ax.plot(x, fp_rc[center-100:center+100],
                        color=color, linewidth=1.0, linestyle="--", alpha=0.5)

        ax.axvline(0, color="black", linewidth=0.7, linestyle=":", alpha=0.4)
        ax.axvspan(-5, 5, alpha=0.06, color="grey")
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.set_ylabel("Predicted accessibility", fontsize=8)
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    axes[-1].set_xlabel("Distance from motif centre (bp)", fontsize=9)

    ct_handles  = [Line2D([0],[0], color=COLORS[ct], linewidth=2, label=ct)
                   for ct in CELLTYPES]
    str_handles = [Line2D([0],[0], color="grey", linewidth=1.5, linestyle="-",  label="forward"),
                   Line2D([0],[0], color="grey", linewidth=1.0, linestyle="--", label="reverse complement")]
    fig.legend(handles=ct_handles + str_handles,
               loc="upper right", fontsize=8, frameon=False,
               bbox_to_anchor=(1.0, 1.0), ncol=1)

    fig.suptitle("ChromBPNet footprints — SOX/HMG motifs\nSOX9_high_P23 vs SOX9_low_P23",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 0.82, 0.96])

    os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches="tight", dpi=150)
    print(f"Saved: {OUT_PDF}")


if __name__ == "__main__":
    main()
