import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repro_common import *
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

THREE = ["MSTO-211H", "NCI-H2052", "NCI-H2452"]
scr = pd.read_csv(f"{FIG}/screen_var5000_log2avg_by_KO.csv", index_col=0)
ncell = pd.read_csv(f"{FIG}/screen_ncells_by_KO.csv").set_index("Var1")["Freq"]
cp10k, meta = load_bulk()

zs = zscore_rows_of(scr); zs = zs.loc[zs.abs().sum(axis=1) > 0]
genes_all = zs.index[corr_order(zs, axis=0)]                       # gene order = screen clustering (all 5000)
genes = [g for g in genes_all if g in cp10k.index]                 # keep genes quantified in bulk
print("genes shown:", len(genes), "of", len(genes_all))

def panel(z_gene_by_ko, counts):
    ko = [k for k in z_gene_by_ko.columns if k != "NTC"]
    o = corr_order(z_gene_by_ko[ko].T, axis=0)
    order = ["NTC"] + [ko[i] for i in o]
    return z_gene_by_ko[order], counts.reindex(order)

panels = [("CROP-seq screen", LINE_COL["CROP-seq"], *panel(zs.loc[genes].T.T, ncell), "cells")]
for cl in THREE:
    g = group_log2avg(cp10k, meta, cl).loc[genes]
    z = g.sub(g.mean(axis=1), axis=0).div(g.std(axis=1, ddof=1), axis=0).fillna(0.0)
    panels.append((f"{cl} bulk (RUN1)", LINE_COL[cl], *panel(z, meta[meta.Cell_Line == cl]["Gene_Control"].value_counts()), "reps"))

widths = [p[2].shape[1] for p in panels]
gap = 2.2
fig = plt.figure(figsize=(13.5, 11), facecolor="white")
gs = GridSpec(1, 5, width_ratios=[widths[0], gap] + widths[1:], wspace=0.06, left=0.05, right=0.97, top=0.80, bottom=0.085)
vmax = 2.5
axes_idx = [0, 2, 3, 4]
for (title, colr, m, n, lab), gi in zip(panels, axes_idx):
    ax = fig.add_subplot(gs[0, gi])
    im = ax.imshow(m.values, aspect="auto", cmap=DIV, vmin=-vmax, vmax=vmax, interpolation="nearest")
    ax.set_xticks(range(m.shape[1]))
    ax.set_xticklabels([f"{k} ({int(n[k])})" for k in m.columns], rotation=90, fontsize=8, color=INK)
    ax.xaxis.tick_top(); ax.tick_params(axis="x", length=0, pad=2)
    ax.axvline(0.5, color="white", lw=3)
    ax.set_yticks([])
    for s in ax.spines.values(): s.set_edgecolor(MUTED); s.set_linewidth(0.6)
    ax.text(0.0, 1.125, title, transform=ax.transAxes, fontsize=11, fontweight="bold", color=colr, va="bottom")
    if gi == 0: ax.set_ylabel(f"{len(genes)} top variable genes (screen order)", fontsize=9, color=INK2)
# side labels
fig.text(0.05, 0.982, "KO pattern of the top variable genes: CROP-seq (left) vs bulk RNA-seq, three cell lines, RUN1 (right)", fontsize=14.5, fontweight="bold", color=INK)
fig.text(0.05, 0.960, "rows = the screen's 5000 variable genes present in bulk, same order in every panel (clustered on the CROP-seq screen); columns = KO groups, NTC first, others clustered per panel;", fontsize=8.5, color=INK2)
fig.text(0.05, 0.946, "values = log2(mean CP10K + 1) per KO group, z-scored per gene across KO groups; number in brackets = cells (screen) or replicates (bulk). NCI-H28 excluded; outlier samples removed (RUN1).", fontsize=8.5, color=INK2)
cax = fig.add_axes([0.05, 0.05, 0.2, 0.008])
cb = fig.colorbar(im, cax=cax, orientation="horizontal"); cb.set_label("gene z-score across KO groups", fontsize=8, color=INK2); cb.ax.tick_params(labelsize=7)
fig.savefig(f"{FIG}/fig1_side_by_side_CROPseq_vs_3lines_RUN1.png", dpi=220, facecolor="white")
fig.savefig(f"{FIG}/fig1_side_by_side_CROPseq_vs_3lines_RUN1.pdf", facecolor="white")
print("saved")
