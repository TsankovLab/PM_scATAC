#!/usr/bin/env python3
"""
plot_HDMA_footprints.py  (R1_Q3)
Averages fold H5s per cell type and plots all HDMA harmonized motif footprints
comparing SOX9_high_P23 vs SOX9_low_P23.

Each motif gets one row: the CWM sequence logo + predicted TF (top TOMTOM match)
on the left, and the ChromBPNet footprint (high vs low) on the right.
Outputs a multi-page PDF.
"""

import argparse, glob, os, re, pickle, csv
import h5py, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_pdf import PdfPages

CELLTYPES = ["SOX9_high_P23", "SOX9_low_P23"]
COLORS    = {"SOX9_high_P23": "#E41A1C", "SOX9_low_P23": "#377EB8"}

OUT_DIR_DEFAULT = "/sc/arion/scratch/giottb01/meso_SOX9_P23_chromBPnet/footprints_HDMA_all"
OUT_PDF = ("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/"
           "git_repo_claude/R1_Q3/Plots/HDMA_all_footprints.pdf")
MOTIF_FILE = "/sc/arion/scratch/giottb01/meso_SOX9_P23_chromBPnet/HDMA_all_motifs_footprint.txt"

# HDMA harmonization output (CWM logos + TOMTOM TF annotations) per modisco type
HDMA_DIRS = {
    "counts":  "/sc/arion/scratch/giottb01/meso_SOX9_P23_chromBPnet/modisco_merged_counts/compiled_tomtom",
    "profile": "/sc/arion/scratch/giottb01/meso_SOX9_P23_chromBPnet/modisco_merged_profile/compiled_tomtom",
}

N_PER_PAGE = 5  # motifs (rows) per PDF page


def load_annotations():
    """Map (type, pattern_id) -> (top_TF_match, qval) from modisco_compiled.tsv."""
    ann = {}
    for typ, d in HDMA_DIRS.items():
        tsv = os.path.join(d, "modisco_compiled.tsv")
        if not os.path.exists(tsv):
            continue
        with open(tsv) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                ann[(typ, row["pattern"])] = (row.get("match0", "") or "NA",
                                              row.get("qval0", "") or "")
    return ann


def find_logo(typ, pattern_id):
    """Path to the forward CWM logo PNG for a merged pattern, or None."""
    hits = glob.glob(os.path.join(HDMA_DIRS[typ], "trimmed_logos",
                                  f"*{pattern_id}.cwm.fwd.png"))
    return hits[0] if hits else None


def key_to_pattern(key):
    """hdma_counts_<pid> / hdma_profile_<pid>  ->  (type, pattern_id)."""
    if key.startswith("hdma_counts_"):
        return "counts", key[len("hdma_counts_"):]
    if key.startswith("hdma_profile_"):
        return "profile", key[len("hdma_profile_"):]
    return None, key


def load_motif_pairs(motif_file):
    """Return list of (fwd_key, rc_key, display_name) from motif file."""
    keys = []
    with open(motif_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2: continue
            keys.append(parts[0])
    # pair fwd with rc
    keyset = set(keys)
    pairs = []
    seen = set()
    for k in keys:
        if k in seen: continue
        if k.endswith('_rc'):
            continue
        rc_key = k + '_rc'
        typ, pid = key_to_pattern(k)
        # clean display name: keep the merged_pattern index so sibling
        # sub-patterns of the same cluster (e.g. Average_65) stay distinguishable
        short = re.sub(r'__merged_pattern_(\d+)$', r' · mp\1', pid)
        name = f"[{typ}] {short}" if typ else k
        pairs.append((k, rc_key if rc_key in keyset else None, name, typ, pid))
        seen.add(k)
        seen.add(rc_key)
    return pairs


def _load_fold_i0(path):
    """Return {motif: i0_array} from one fold H5.

    chrombpnet footprints writes two layouts depending on the motif set:
      (a) one HDF5 group per motif with datasets 'i0' (footprint) and 'i1' (scalar)
      (b) a single 'data' dataset holding a pickled {motif: [i0_array, i1_scalar]}
          (seen with the large HDMA motif set; some names contain '.').
    Handle both so averaging works regardless of layout.
    """
    out = {}
    with h5py.File(path, "r") as f:
        keys = list(f.keys())
        if keys == ["data"]:
            obj = pickle.loads(f["data"][0].tobytes())
            for motif, val in obj.items():
                out[motif] = np.asarray(val[0])
        else:
            for grp in keys:
                if isinstance(f[grp], h5py.Group) and "i0" in f[grp]:
                    out[grp] = f[grp]["i0"][()]
    return out


def average_folds(out_dir, celltype):
    files = sorted(glob.glob(os.path.join(out_dir, f"{celltype}_fold*_footprints.h5")))
    if not files:
        return {}
    per_fold = [_load_fold_i0(fn) for fn in files]
    motifs = list(per_fold[0].keys())
    result = {}
    for grp in motifs:
        arrays = [d[grp] for d in per_fold if grp in d]
        if arrays:
            result[grp] = np.mean(arrays, axis=0)
    avg_path = os.path.join(out_dir, f"{celltype}_average_footprints.h5")
    with h5py.File(avg_path, "w") as f:
        for grp, arr in result.items():
            g = f.create_group(grp)
            g.create_dataset("i0", data=arr)
    print(f"  [{celltype}] averaged {len(files)} folds, {len(result)} motifs")
    return result


def _fmt_q(q):
    try:
        return f"q={float(q):.1e}"
    except (TypeError, ValueError):
        return ""


def _draw_logo(ax, logo_path):
    """Render a CWM logo PNG into ax (or a placeholder if missing)."""
    if logo_path and os.path.exists(logo_path):
        ax.imshow(mpimg.imread(logo_path), aspect="auto")
    else:
        ax.text(0.5, 0.5, "(logo n/a)", ha="center", va="center",
                fontsize=7, color="grey", transform=ax.transAxes)
    ax.axis("off")


def _draw_footprint(ax, data, fwd_key, rc_key, x):
    for ct in CELLTYPES:
        if not data[ct]:
            continue
        color = COLORS[ct]
        center = None
        if fwd_key in data[ct]:
            fp = data[ct][fwd_key]
            center = fp.shape[0] // 2
            ax.plot(x, fp[center-100:center+100],
                    color=color, linewidth=1.4, linestyle="-", alpha=0.9)
        if rc_key and rc_key in data[ct] and center is not None:
            fp_rc = data[ct][rc_key]
            ax.plot(x, fp_rc[center-100:center+100],
                    color=color, linewidth=0.9, linestyle="--", alpha=0.5)
    ax.axvline(0, color="black", linewidth=0.7, linestyle=":", alpha=0.4)
    ax.axvspan(-5, 5, alpha=0.06, color="grey")
    ax.set_ylabel("Accessibility", fontsize=7)
    ax.tick_params(labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out_dir",    default=OUT_DIR_DEFAULT)
    p.add_argument("--motif_file", default=MOTIF_FILE)
    args = p.parse_args()

    data = {ct: average_folds(args.out_dir, ct) for ct in CELLTYPES}
    ann = load_annotations()
    motif_pairs = load_motif_pairs(args.motif_file)
    print(f"Plotting {len(motif_pairs)} motifs (logo + TF + footprint)...")

    x = np.arange(-100, 100)
    ct_handles  = [Line2D([0],[0], color=COLORS[ct], linewidth=2, label=ct) for ct in CELLTYPES]
    str_handles = [Line2D([0],[0], color="grey", linewidth=1.5, linestyle="-",  label="forward"),
                   Line2D([0],[0], color="grey", linewidth=1.0, linestyle="--", label="reverse complement")]

    os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)
    pages = [motif_pairs[i:i+N_PER_PAGE] for i in range(0, len(motif_pairs), N_PER_PAGE)]

    with PdfPages(OUT_PDF) as pdf:
        for page_motifs in pages:
            n = len(page_motifs)
            h = n * 1.7 + 1.4                       # extra headroom for title + legend
            fig = plt.figure(figsize=(9, h))
            gs = GridSpec(n, 2, width_ratios=[1.5, 2.2], wspace=0.18, hspace=0.85,
                          left=0.04, right=0.99, top=1 - 1.05 / h, bottom=0.06)

            for i, (fwd_key, rc_key, name, typ, pid) in enumerate(page_motifs):
                tf, q = ann.get((typ, pid), ("NA", ""))
                qtxt = _fmt_q(q)

                # left: CWM logo + predicted TF label (pattern id below the logo)
                ax_logo = fig.add_subplot(gs[i, 0])
                _draw_logo(ax_logo, find_logo(typ, pid))
                ax_logo.set_title(f"{tf}" + (f"  ({qtxt})" if qtxt else ""),
                                  fontsize=8.5, fontweight="bold", color="#222222")
                ax_logo.text(0.0, -0.12, name, transform=ax_logo.transAxes,
                             fontsize=6.5, color="grey", va="top")

                # right: footprint high vs low
                ax_fp = fig.add_subplot(gs[i, 1])
                _draw_footprint(ax_fp, data, fwd_key, rc_key, x)
                if i == n - 1:
                    ax_fp.set_xlabel("Distance from motif centre (bp)", fontsize=8)

            fig.suptitle("ChromBPNet footprints — HDMA harmonized motifs   "
                         "(CWM logo + predicted TF + footprint)\n"
                         "SOX9_high_P23 vs SOX9_low_P23",
                         fontsize=9.5, fontweight="bold", y=0.995)
            fig.legend(handles=ct_handles + str_handles,
                       loc="upper center", fontsize=7, frameon=False,
                       bbox_to_anchor=(0.5, 1 - 0.62 / h), ncol=4)
            pdf.savefig(fig, dpi=150)
            plt.close(fig)

    print(f"Saved: {OUT_PDF}")


if __name__ == "__main__":
    main()
