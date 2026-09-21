import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repro_common import *
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

THREE = ["MSTO-211H", "NCI-H2052", "NCI-H2452"]
ORDER = ["TCF3", "TEAD4", "TWIST1", "PITX1", "SOX9", "TEAD2", "BPTF", "HMGA1"]      # KOs in both, CROP-seq (fig1) order
EXTRA = ["MEF2A", "MEF2D"]                                                          # screen-only
cp10k, meta = load_bulk()

# ---------- KO - NTC profiles (same construction as fig3) ----------
scr_cm = -pd.read_csv(f"{PP}/screen_cm_diff_matrix.csv", index_col=0); cm_cols = list(scr_cm.columns)
sc = pd.read_csv(f"{PP}/cm_sample_scores_AddModuleScore_bulk_RUN1.csv").set_index("sample")[cm_cols]
mod = {"CROP-seq": scr_cm.T}
for cl in THREE:
    m = meta[meta.Cell_Line == cl]; ntc = sc.loc[m.index[m.Gene_Control == "NTC"]].median()
    mod[cl] = pd.DataFrame({g: sc.loc[i].mean() - ntc for g, i in m[m.Gene_Control != "NTC"].groupby("Gene_Control").groups.items()})
scr_g = pd.read_csv(f"{FIG}/screen_var5000_log2avg_by_KO.csv", index_col=0)
genes = [g for g in scr_g.index if g in cp10k.index]
gen = {"CROP-seq": scr_g.loc[genes].drop(columns="NTC").sub(scr_g.loc[genes, "NTC"], axis=0)}
for cl in THREE:
    g = group_log2avg(cp10k, meta, cl).loc[genes]; gen[cl] = g.drop(columns="NTC").sub(g["NTC"], axis=0)

def normalise(D, centre):
    out = {}
    for k, m in D.items():
        x = m.sub(m.mean(axis=1), axis=0) if centre else m
        out[k] = x.div(np.sqrt((x ** 2).mean(axis=1)).replace(0, np.nan), axis=0).fillna(0)
    return out

def cross(N):
    cols = ORDER + EXTRA; S = N["CROP-seq"][cols]
    rows, lab = [], []
    for cl in THREE:
        for ko in ORDER:
            if ko in N[cl].columns:
                rows.append([np.corrcoef(N[cl][ko].values, S[c].values)[0, 1] for c in cols]); lab.append((cl, ko))
    return np.array(rows), lab, cols

def draw(feat, title, fname):
    res = {c: cross(normalise(feat, c)) for c in [False, True]}
    vmax = np.ceil(np.percentile(np.abs(np.concatenate([res[c][0].ravel() for c in res])), 99) * 10) / 10
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 9.2), facecolor="white")
    stats = []
    for ax, centre in zip(axes, [False, True]):
        C, lab, cols = res[centre]
        im = ax.imshow(C, cmap=DIV, vmin=-vmax, vmax=vmax, aspect="auto", interpolation="nearest")
        ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=90, fontsize=9)
        ax.set_yticks(range(len(lab))); ax.set_yticklabels([k for _, k in lab], fontsize=8.5)
        ax.xaxis.tick_top()
        for i, (cl, ko) in enumerate(lab):
            for j, c in enumerate(cols):
                ax.text(j, i, f"{C[i, j]:+.2f}".replace("+0.", "+.").replace("-0.", "-."), ha="center", va="center", fontsize=6.3, color=INK2 if abs(C[i, j]) < 0.6 * vmax else "white")
            ax.add_patch(Rectangle((ORDER.index(ko) - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=INK, lw=1.5))
        ax.axvline(len(ORDER) - 0.5, color="white", lw=4)
        edges = [i for i in range(1, len(lab)) if lab[i][0] != lab[i - 1][0]]
        for e in edges: ax.axhline(e - 0.5, color="white", lw=4)
        bounds = [0] + edges + [len(lab)]
        for a, b, cl in zip(bounds[:-1], bounds[1:], THREE):
            ax.add_patch(Rectangle((-1.05, a - 0.5), 0.35, b - a, color=LINE_COL[cl], clip_on=False))
            ax.text(-1.25, (a + b - 1) / 2, cl, rotation=90, ha="right", va="center", fontsize=9.5, fontweight="bold", color=LINE_COL[cl])
        for s in ax.spines.values(): s.set_visible(False)
        match = np.array([C[i, ORDER.index(ko)] for i, (_, ko) in enumerate(lab)])
        mask = np.ones_like(C, bool)
        for i, (_, ko) in enumerate(lab): mask[i, ORDER.index(ko)] = False; mask[i, len(ORDER):] = False
        ttl = "raw (KO − NTC)" if not centre else "centred across KOs (shared NTC / global component removed)"
        ax.set_title(f"{ttl}\nmatched r = {match.mean():+.2f}  vs  other-KO r = {C[mask].mean():+.2f}", fontsize=10, color=INK, pad=48, loc="left")
        stats.append((ttl, match.mean(), C[mask].mean()))
    cax = fig.add_axes([0.35, 0.055, 0.3, 0.012]); cb = fig.colorbar(im, cax=cax, orientation="horizontal"); cb.set_label("Pearson r between bulk KO profile and CROP-seq KO profile", fontsize=9)
    fig.suptitle(title, x=0.02, ha="left", fontsize=14.5, fontweight="bold", y=0.985)
    fig.text(0.02, 0.945, "rows = bulk KO (RUN1) per cell line, columns = CROP-seq KO, same KO order (matched pair outlined: a diagonal = correspondence); MEF2A/MEF2D exist only in the screen.", fontsize=9, color=INK2)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.80, bottom=0.11, wspace=0.16)
    fig.savefig(f"{FIG}/{fname}.png", dpi=190, facecolor="white"); fig.savefig(f"{FIG}/{fname}.pdf", facecolor="white"); plt.close(fig)
    # per-line diagnostics
    C, lab, cols = res[False]
    rows = []
    for cl in THREE:
        idx = [i for i, (c, _) in enumerate(lab) if c == cl]
        ranks = [int((C[i, :len(ORDER)] > C[i, ORDER.index(lab[i][1])]).sum() + 1) for i in idx]
        rows.append((cl, len(idx), np.mean([C[i, ORDER.index(lab[i][1])] for i in idx]), np.median(ranks), sum(r == 1 for r in ranks)))
    print(title); print(pd.DataFrame(rows, columns=["line", "n", "mean matched r", "median rank (of 8)", "n rank1"]).round(3).to_string(index=False)); print(stats)

draw(mod, "Bulk KO vs CROP-seq KO similarity: 20 Cm modules", "fig4_bulk_vs_CROPseq_corr_modules_RUN1")
draw(gen, f"Bulk KO vs CROP-seq KO similarity: screen's top variable genes (n={len(genes)})", "fig4_bulk_vs_CROPseq_corr_variable_genes_RUN1")
