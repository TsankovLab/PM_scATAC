import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repro_common import *
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

scr = pd.read_csv(f"{FIG}/screen_var5000_log2avg_by_KO.csv", index_col=0)      # genes x KO
ncell = pd.read_csv(f"{FIG}/screen_ncells_by_KO.csv").set_index("Var1")["Freq"]
cp10k, meta = load_bulk()
print("screen genes:", scr.shape, "| in bulk:", scr.index.isin(cp10k.index).sum())

# screen: scaled across KOs, drop zero-variance genes
zs = zscore_rows_of(scr); zs = zs.loc[zs.abs().sum(axis=1) > 0]
col_order = corr_order(zs, axis=0); genes = zs.index[col_order]                # gene order from the screen
print("screen genes kept after variance filter:", len(genes))

panels = []   # (title, color, matrix[KO x gene], n_annot, n_label)
def prep(z_ko_by_gene, ncounts):
    ko = [k for k in z_ko_by_gene.index if k != "NTC"]
    ro = corr_order(z_ko_by_gene.loc[ko], axis=0)
    order = (["NTC"] if "NTC" in z_ko_by_gene.index else []) + [ko[i] for i in ro]
    return z_ko_by_gene.loc[order], ncounts.reindex(order)
m, n = prep(zs.loc[genes].T, ncell)
panels.append(("CROP-seq screen", LINE_COL["CROP-seq"], m, n, "cells"))
for cl in LINES:
    g = group_log2avg(cp10k, meta, cl)
    g = g.reindex(genes)                                       # keep screen gene order; missing -> NaN
    miss = g.isna().all(axis=1).sum()
    z = zscore_rows_of(g.fillna(0).where(~g.isna().all(axis=1), np.nan)) if False else None
    gg = g.copy()
    z = gg.sub(gg.mean(axis=1), axis=0).div(gg.std(axis=1, ddof=1), axis=0).fillna(0.0)
    nrep = meta[meta["Cell_Line"] == cl]["Gene_Control"].value_counts()
    m, n = prep(z.T, nrep)
    panels.append((f"{cl} bulk (RUN1)  [{miss} genes absent]", LINE_COL[cl], m, n, "reps"))

heights = [p[2].shape[0] for p in panels]
fig = plt.figure(figsize=(14, 0.24 * sum(heights) + 4.0), facecolor="white")
gs = GridSpec(len(panels), 2, width_ratios=[40, 1.6], height_ratios=heights, hspace=0.55, wspace=0.03,
              left=0.10, right=0.95, top=0.925, bottom=0.06)
vmax = 2.5
for i, (title, colr, m, n, lab) in enumerate(panels):
    ax = fig.add_subplot(gs[i, 0]); axb = fig.add_subplot(gs[i, 1])
    im = ax.imshow(m.values, aspect="auto", cmap=DIV, vmin=-vmax, vmax=vmax, interpolation="nearest")
    ax.set_yticks(range(m.shape[0])); ax.set_yticklabels(m.index, fontsize=8.5, color=INK)
    ax.set_xticks([]); 
    if "NTC" in m.index: ax.axhline(0.5, color="white", linewidth=3)
    for s in ax.spines.values(): s.set_edgecolor(MUTED); s.set_linewidth(0.6)
    ax.set_title(title, loc="left", fontsize=10.5, fontweight="bold", color=colr, pad=4)
    axb.barh(range(m.shape[0]), n.values, color=colr, height=0.7); axb.set_ylim(m.shape[0] - 0.5, -0.5)
    axb.set_yticks([]); axb.tick_params(labelsize=7, colors=INK2)
    for s in ["top", "right", "left"]: axb.spines[s].set_visible(False)
    axb.set_xlabel(lab, fontsize=7, color=INK2, labelpad=1)
cax = fig.add_axes([0.10, 0.03, 0.25, 0.006])
cb = fig.colorbar(im, cax=cax, orientation="horizontal"); cb.set_label("gene z-score across KO groups", fontsize=8, color=INK2); cb.ax.tick_params(labelsize=7)
fig.text(0.10, 0.972, "KO pattern of the top-5000 variable genes: CROP-seq screen vs bulk RNA-seq", fontsize=15, fontweight="bold", color=INK)
fig.text(0.10, 0.953, "columns = the screen's 5000 variable genes (same order in every panel, clustered on the screen); rows = KO groups (NTC on top, others clustered per panel)",
         fontsize=8.5, color=INK2)
fig.text(0.10, 0.940, "values = log2(mean CP10K + 1) per KO group, z-scored per gene across KO groups (as in the original heatmap); bars = number of cells (screen) or replicates (bulk)",
         fontsize=8.5, color=INK2)
fig.savefig(f"{FIG}/fig1_KO_pattern_top5000_variable_genes.png", dpi=200, facecolor="white")
fig.savefig(f"{FIG}/fig1_KO_pattern_top5000_variable_genes.pdf", facecolor="white")
print("saved fig1")
