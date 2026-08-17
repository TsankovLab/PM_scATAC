#!/usr/bin/env python3
"""plot_sox9_marginal_allfolds.py
Aggregate SOX9 marginal-footprint amplitudes across ALL trained folds for the
SOX9pos vs SOX9neg ChromBPNet models.

For each fold we already computed per-model marginal footprints of the SOX motif
panel (chrombpnet footprints, *_footprints.h5). Here we:
  * read the predicted marginal footprint per motif per fold,
  * compute a footprint-amplitude metric (max-min of the central 100 bp),
  * summarise per motif as mean +/- SD across folds for SOX9pos and SOX9neg,
  * report the per-fold pos/neg ratio (mean +/- SD) and a paired t-test across
    folds (is the SOX9pos amplitude different from SOX9neg?).

Outputs:
  --out_pdf : one panel per motif, mean footprint across folds with +/-SD shading
  --out_csv : per-motif per-fold amplitudes + summary stats
"""
import argparse, os, glob, re
import h5py, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    from scipy.stats import ttest_rel
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False


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


def fold_of(path):
    m = re.search(r"fold(\d+)", os.path.basename(path))
    return int(m.group(1)) if m else -1


def collect(template):
    """template has a {ct} slot; returns {fold: {motif: profile}}"""
    res = {}
    for ct_path in sorted(glob.glob(template)):
        res[fold_of(ct_path)] = load(ct_path)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="sox9_marginal_test dir")
    ap.add_argument("--out_pdf", required=True)
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()

    pos = collect(os.path.join(args.dir, "SOX9pos_fold*_footprints.h5"))
    neg = collect(os.path.join(args.dir, "SOX9neg_fold*_footprints.h5"))
    folds = sorted(set(pos) & set(neg))
    if not folds:
        raise SystemExit("No folds present for BOTH cell types in " + args.dir)
    print(f"Folds with both models: {folds}")

    # motifs = real SOX motifs present in every fold's pos+neg h5 (exclude rc/control)
    common = None
    for f in folds:
        ms = {m for m in pos[f] if not m.endswith("_rc") and m != "control"}
        ms &= {m for m in neg[f] if not m.endswith("_rc") and m != "control"}
        common = ms if common is None else (common & ms)
    motifs = sorted(common)
    print(f"Motifs: {motifs}")

    # amplitudes: motif -> {fold -> (amp_pos, amp_neg)}
    rows = []
    summary = {}
    for m in motifs:
        ap_list, an_list, ratios = [], [], []
        for f in folds:
            a_p, a_n = amp(pos[f][m]), amp(neg[f][m])
            ap_list.append(a_p); an_list.append(a_n)
            ratios.append(a_p / a_n if a_n > 0 else np.nan)
            rows.append((m, f, round(a_p, 6), round(a_n, 6),
                         round(a_p / a_n, 4) if a_n > 0 else np.nan))
        ap_arr, an_arr = np.array(ap_list), np.array(an_list)
        pval = ttest_rel(ap_arr, an_arr).pvalue if (HAVE_SCIPY and len(folds) > 1) else np.nan
        summary[m] = dict(
            pos_mean=ap_arr.mean(), pos_sd=ap_arr.std(ddof=1) if len(folds) > 1 else 0.0,
            neg_mean=an_arr.mean(), neg_sd=an_arr.std(ddof=1) if len(folds) > 1 else 0.0,
            ratio_mean=np.nanmean(ratios), ratio_sd=np.nanstd(ratios, ddof=1) if len(folds) > 1 else 0.0,
            pval=pval)

    # ---- CSV ----
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
    with open(args.out_csv, "w") as fh:
        fh.write("motif,fold,amp_SOX9pos,amp_SOX9neg,ratio_pos_over_neg\n")
        for r in rows:
            fh.write(",".join(map(str, r)) + "\n")
        fh.write("\n# summary across folds (mean +/- SD)\n")
        fh.write("motif,pos_mean,pos_sd,neg_mean,neg_sd,ratio_mean,ratio_sd,paired_t_pval\n")
        for m in motifs:
            s = summary[m]
            fh.write(f"{m},{s['pos_mean']:.6f},{s['pos_sd']:.6f},"
                     f"{s['neg_mean']:.6f},{s['neg_sd']:.6f},"
                     f"{s['ratio_mean']:.4f},{s['ratio_sd']:.4f},{s['pval']:.4g}\n")

    # ---- PDF: mean footprint +/- SD across folds ----
    x = np.arange(-100, 100)
    fig, axes = plt.subplots(len(motifs), 1, figsize=(6.2, 2.7 * len(motifs)), squeeze=False)
    for ax, m in zip(axes[:, 0], motifs):
        for grp, store, col, lab in [("pos", pos, "#E41A1C", "SOX9pos"),
                                     ("neg", neg, "#377EB8", "SOX9neg")]:
            mat = []
            for f in folds:
                p = store[f][m]; c = len(p) // 2
                mat.append(p[c - 100:c + 100])
            mat = np.vstack(mat)
            mean, sd = mat.mean(0), mat.std(0, ddof=1) if len(folds) > 1 else (mat.mean(0), np.zeros(200))
            ax.plot(x, mean, color=col, lw=1.7, label=lab)
            ax.fill_between(x, mean - sd, mean + sd, color=col, alpha=0.18, lw=0)
        s = summary[m]
        ax.axvline(0, color="k", lw=0.6, ls=":", alpha=0.4)
        ax.set_title(f"{m}   pos={s['pos_mean']:.4f}±{s['pos_sd']:.4f}  "
                     f"neg={s['neg_mean']:.4f}±{s['neg_sd']:.4f}  "
                     f"ratio={s['ratio_mean']:.2f}  p={s['pval']:.2g}",
                     fontsize=8.5, fontweight="bold")
        ax.set_ylabel("pred. accessibility", fontsize=8)
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(fontsize=8, frameon=False)
    axes[-1, 0].set_xlabel("Distance from motif centre (bp)", fontsize=9)
    fig.suptitle(f"SOX9 marginal footprints across folds {folds} — SOX9pos vs SOX9neg\n"
                 "mean ± SD across folds; amplitude = max-min of central 100 bp",
                 fontsize=10, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(os.path.dirname(args.out_pdf), exist_ok=True)
    fig.savefig(args.out_pdf, dpi=150)

    print(f"Wrote {args.out_pdf} and {args.out_csv}")
    for m in motifs:
        s = summary[m]
        print(f"  {m:18s} pos={s['pos_mean']:.4f}±{s['pos_sd']:.4f} "
              f"neg={s['neg_mean']:.4f}±{s['neg_sd']:.4f} "
              f"ratio={s['ratio_mean']:.2f}±{s['ratio_sd']:.2f} p={s['pval']:.3g}")


if __name__ == "__main__":
    main()
