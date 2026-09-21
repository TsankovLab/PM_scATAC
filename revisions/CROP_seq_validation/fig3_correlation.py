import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repro_common import *
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

cp10k, meta = load_bulk()

# ------------------------ build KO profiles (KO - NTC) ------------------------
# (a) 20 Cm modules
scr_cm = -pd.read_csv(f"{PP}/screen_cm_diff_matrix.csv", index_col=0)                # stored as NTC-KO -> flip to KO-NTC
cm_cols = list(scr_cm.columns)
sc = pd.read_csv(f"{PP}/cm_sample_scores_AddModuleScore_bulk_RUN1.csv")
sc = sc.set_index("sample")[cm_cols]
mod = {"CROP-seq": scr_cm.T}                                                        # modules x KO
for cl in LINES:
    m = meta[meta.Cell_Line == cl]; ntc = sc.loc[m.index[m.Gene_Control == "NTC"]].median()
    mod[cl] = pd.DataFrame({g: sc.loc[idx].mean() - ntc for g, idx in m[m.Gene_Control != "NTC"].groupby("Gene_Control").groups.items()})

# (b) the screen's 5000 variable genes (log2(mean CP10K+1), KO - NTC)
scr_g = pd.read_csv(f"{FIG}/screen_var5000_log2avg_by_KO.csv", index_col=0)
genes = [g for g in scr_g.index if g in cp10k.index]
gen = {"CROP-seq": scr_g.loc[genes].drop(columns="NTC").sub(scr_g.loc[genes, "NTC"], axis=0)}
for cl in LINES:
    g = group_log2avg(cp10k, meta, cl).loc[genes]
    gen[cl] = g.drop(columns="NTC").sub(g["NTC"], axis=0)
print("modules:", {k: v.shape for k, v in mod.items()}); print("genes:", {k: v.shape for k, v in gen.items()})

def normalise(D, centre):
    """per-dataset: optionally centre each feature across KOs (removes shared NTC error), then scale features to equal RMS."""
    out = {}
    for k, m in D.items():
        x = m.sub(m.mean(axis=1), axis=0) if centre else m
        out[k] = x.div(np.sqrt((x ** 2).mean(axis=1)).replace(0, np.nan), axis=0).fillna(0)
    return out

def stack(D):
    cols, lab = [], []
    for ds in ["CROP-seq"] + LINES:
        for ko in sorted(D[ds].columns):
            cols.append(D[ds][ko].values); lab.append((ds, ko))
    return np.array(cols), lab

def draw(C, lab, title, sub, fname):
    n = len(lab); fig, ax = plt.subplots(figsize=(13.2, 12.2), facecolor="white")
    im = ax.imshow(C, cmap=DIV, vmin=-1, vmax=1, interpolation="nearest")
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels([k for _, k in lab], rotation=90, fontsize=6.5); ax.set_yticklabels([k for _, k in lab], fontsize=6.5)
    ds_arr = np.array([d for d, _ in lab]); edges = [0] + [i for i in range(1, n) if ds_arr[i] != ds_arr[i - 1]] + [n]
    for e in edges[1:-1]:
        ax.axhline(e - 0.5, color="white", lw=2.2); ax.axvline(e - 0.5, color="white", lw=2.2)
    for a, b in zip(edges[:-1], edges[1:]):
        ds = ds_arr[a]
        ax.add_patch(plt.Rectangle((-2.6, a - 0.5), 1.5, b - a, color=LINE_COL[ds], clip_on=False, transform=ax.transData))
        ax.add_patch(plt.Rectangle((a - 0.5, -2.6), b - a, 1.5, color=LINE_COL[ds], clip_on=False, transform=ax.transData))
        ax.text(-3.2, (a + b - 1) / 2, ds, rotation=90, ha="right", va="center", fontsize=9, fontweight="bold", color=LINE_COL[ds])
        ax.text((a + b - 1) / 2, -3.5, ds, ha="center", va="bottom", fontsize=9, fontweight="bold", color=LINE_COL[ds])
    for s in ax.spines.values(): s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.set_label("Pearson r between KO profiles", fontsize=9)
    fig.text(0.02, 0.985, title, fontsize=15, fontweight="bold", color=INK, va="top"); fig.text(0.02, 0.955, sub, fontsize=9, color=INK2, va="top")
    fig.subplots_adjust(top=0.86, left=0.09, right=0.93, bottom=0.06)
    fig.savefig(f"{FIG}/{fname}.png", dpi=180, facecolor="white"); fig.savefig(f"{FIG}/{fname}.pdf", facecolor="white"); plt.close(fig)

def matched_rank(C, lab, tag):
    ds_arr = np.array([d for d, _ in lab]); ko_arr = np.array([k for _, k in lab])
    si = np.where(ds_arr == "CROP-seq")[0]; rows = []
    for i in np.where(ds_arr != "CROP-seq")[0]:
        if ko_arr[i] not in ko_arr[si]: continue
        r = C[i, si]; same = r[list(ko_arr[si]).index(ko_arr[i])]
        rows.append([tag, ds_arr[i], ko_arr[i], same, int((r > same).sum() + 1), len(si), np.mean(np.delete(r, list(ko_arr[si]).index(ko_arr[i])))])
    return pd.DataFrame(rows, columns=["features", "cell_line", "KO", "r_same_KO", "rank_of_same_KO", "n_screen_KOs", "mean_r_other_KOs"])

summ = []
for feat, D in [("modules", mod), ("variable_genes", gen)]:
    for centre in [False, True]:
        N = normalise(D, centre); X, lab = stack(N)
        C = np.corrcoef(X); tag = f"{feat}_{'centred' if centre else 'raw'}"
        pd.DataFrame(C, index=[f"{a}|{b}" for a, b in lab], columns=[f"{a}|{b}" for a, b in lab]).to_csv(f"{FIG}/fig3_corr_{tag}.csv")
        ttl = {"modules": "KO–KO similarity from the 20 Cm modules", "variable_genes": f"KO–KO similarity from the screen's top variable genes (n={len(genes)})"}[feat]
        sub = ("profiles = KO − NTC per dataset; each feature scaled to equal RMS within dataset" +
               ("; features centred across KOs within each dataset first (removes the shared NTC/global component)" if centre else " (NTC noise is shared by every KO of a line)"))
        draw(C, lab, ttl + ("  [centred]" if centre else ""), sub, f"fig3_corr_{tag}")
        mr = matched_rank(C, lab, tag); summ.append(mr)
res = pd.concat(summ); res.to_csv(f"{FIG}/fig3_matched_KO_rank.csv", index=False)
g = res.groupby("features").agg(n=("KO", "size"), top1=("rank_of_same_KO", lambda r: int((r == 1).sum())), top3=("rank_of_same_KO", lambda r: int((r <= 3).sum())),
                                 median_rank=("rank_of_same_KO", "median"), mean_r_same=("r_same_KO", "mean"), mean_r_other=("mean_r_other_KOs", "mean"))
print(g.round(3).to_string())
print(res[res.features == "modules_centred"].groupby("cell_line").rank_of_same_KO.median())
