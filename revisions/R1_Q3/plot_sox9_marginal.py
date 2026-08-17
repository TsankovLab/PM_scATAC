#!/usr/bin/env python3
"""plot_sox9_marginal.py
Compare SOX9 marginal footprints between the SOX9pos and SOX9neg ChromBPNet
models (single fold). Overlays the predicted marginal footprint per SOX motif
and reports a footprint-amplitude metric (max-min of the central window) as a
quick read on whether/how strongly each model has learned SOX9 responsiveness.
"""
import argparse, os
import h5py, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(path):
    out = {}
    with h5py.File(path, "r") as f:
        for k in f.keys():
            if isinstance(f[k], h5py.Group) and "i0" in f[k]:
                out[k] = f[k]["i0"][()]
    return out


def amp(profile, half=50):
    c = len(profile) // 2
    w = profile[c - half:c + half]
    return float(w.max() - w.min())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pos", required=True)
    ap.add_argument("--neg", required=True)
    ap.add_argument("--out_pdf", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--fold", default="?")
    args = ap.parse_args()

    P, N = load(args.pos), load(args.neg)
    motifs = [m for m in P if not m.endswith("_rc") and m != "control"]
    motifs = [m for m in motifs if m in N]

    x = np.arange(-100, 100)
    rows = []
    fig, axes = plt.subplots(len(motifs), 1, figsize=(6, 2.6 * len(motifs)), squeeze=False)
    for ax, m in zip(axes[:, 0], motifs):
        cp, cn = len(P[m]) // 2, len(N[m]) // 2
        ax.plot(x, P[m][cp-100:cp+100], color="#E41A1C", lw=1.6, label="SOX9pos")
        ax.plot(x, N[m][cn-100:cn+100], color="#377EB8", lw=1.6, label="SOX9neg")
        ax.axvline(0, color="k", lw=0.6, ls=":", alpha=0.4)
        ap_, an_ = amp(P[m]), amp(N[m])
        ratio = ap_ / an_ if an_ > 0 else float("nan")
        rows.append((m, round(ap_, 5), round(an_, 5), round(ratio, 3)))
        ax.set_title(f"{m}   amp pos={ap_:.4f}  neg={an_:.4f}  (pos/neg={ratio:.2f})",
                     fontsize=9, fontweight="bold")
        ax.set_ylabel("pred. accessibility", fontsize=8)
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(fontsize=8, frameon=False)
    axes[-1, 0].set_xlabel("Distance from motif centre (bp)", fontsize=9)
    fig.suptitle(f"SOX9 marginal footprints — SOX9pos vs SOX9neg (fold {args.fold})\n"
                 "amplitude = max-min of central 100 bp of the predicted footprint",
                 fontsize=10, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs(os.path.dirname(args.out_pdf), exist_ok=True)
    fig.savefig(args.out_pdf, dpi=150)

    with open(args.out_csv, "w") as f:
        f.write("motif,amp_SOX9pos,amp_SOX9neg,ratio_pos_over_neg\n")
        for r in rows:
            f.write(",".join(map(str, r)) + "\n")
    print(f"Wrote {args.out_pdf} and {args.out_csv}")
    for r in rows:
        print(f"  {r[0]:20s} pos={r[1]:.4f} neg={r[2]:.4f} ratio={r[3]}")


if __name__ == "__main__":
    main()
