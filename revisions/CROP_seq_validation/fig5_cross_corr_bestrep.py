import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repro_common import *
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

THREE = ["MSTO-211H", "NCI-H2052", "NCI-H2452"]
ORDER = ["TCF3", "TEAD4", "TWIST1", "PITX1", "SOX9", "TEAD2", "BPTF", "HMGA1"]
EXTRA = ["MEF2A", "MEF2D"]
cp10k, meta = load_bulk()
lg = np.log2(cp10k + 1)

# ---------------- per-REPLICATE profiles (sample - NTC) ----------------
scr_cm = -pd.read_csv(f"{PP}/screen_cm_diff_matrix.csv", index_col=0); cm_cols = list(scr_cm.columns)
sc = pd.read_csv(f"{PP}/cm_sample_scores_AddModuleScore_bulk_RUN1.csv").set_index("sample")[cm_cols]
scr_g = pd.read_csv(f"{FIG}/screen_var5000_log2avg_by_KO.csv", index_col=0)
genes = [g for g in scr_g.index if g in cp10k.index]

feats = {}
# modules
rep = {}; avg = {}
for cl in THREE:
    m = meta[meta.Cell_Line == cl]; ntc = sc.loc[m.index[m.Gene_Control == "NTC"]].median()
    rep[cl] = {s: sc.loc[s] - ntc for s in m.index[m.Gene_Control != "NTC"]}
    avg[cl] = pd.DataFrame({g: sc.loc[i].mean() - ntc for g, i in m[m.Gene_Control != "NTC"].groupby("Gene_Control").groups.items()})
feats["modules"] = (scr_cm.T, rep, avg)
# variable genes
rep = {}; avg = {}
for cl in THREE:
    m = meta[meta.Cell_Line == cl]; ntc = np.log2(cp10k.loc[genes, m.index[m.Gene_Control == "NTC"]].mean(axis=1) + 1)
    rep[cl] = {s: lg.loc[genes, s] - ntc for s in m.index[m.Gene_Control != "NTC"]}
    g = group_log2avg(cp10k, meta, cl).loc[genes]; avg[cl] = g.drop(columns="NTC").sub(g["NTC"], axis=0)
feats["variable_genes"] = (scr_g.loc[genes].drop(columns="NTC").sub(scr_g.loc[genes, "NTC"], axis=0), rep, avg)

def rms_scale(m):   # per-feature RMS across KOs (dataset-wise)
    return np.sqrt((m ** 2).mean(axis=1)).replace(0, np.nan)

def build(name):
    S, rep, avg = feats[name]
    sS = S.div(rms_scale(S), axis=0).fillna(0)
    chosen = {"CROP-seq": S}; picks = {}; allrep_r = {}
    for cl in THREE:
        sc_ = rms_scale(avg[cl]); m = meta[meta.Cell_Line == cl]; cols = {}
        for ko in ORDER:
            samples = [s for s in m.index[m.Gene_Control == ko] if s in rep[cl]]
            if not samples: continue
            rs = {s: np.corrcoef((rep[cl][s] / sc_).fillna(0).values, sS[ko].values)[0, 1] for s in samples}
            best = max(rs, key=rs.get); cols[ko] = rep[cl][best]; picks[(cl, ko)] = (best, rs)
            allrep_r[(cl, ko)] = {s: (rep[cl][s] / sc_).fillna(0).values for s in samples}
        chosen[cl] = pd.DataFrame(cols)
    return chosen, picks, allrep_r, sS

def normalise(D, centre):
    out = {}
    for k, m in D.items():
        x = m.sub(m.mean(axis=1), axis=0) if centre else m
        out[k] = x.div(rms_scale(x), axis=0).fillna(0)
    return out

def cross(N, cols):
    S = N["CROP-seq"][cols]; rows, lab = [], []
    for cl in THREE:
        for ko in ORDER:
            if ko in N[cl].columns:
                rows.append([np.corrcoef(N[cl][ko].values, S[c].values)[0, 1] for c in cols]); lab.append((cl, ko))
    return np.array(rows), lab

def draw(name, title, fname):
    chosen, picks, allrep, sS = build(name); cols = ORDER + EXTRA
    res = {c: cross(normalise(chosen, c), cols) for c in [False, True]}
    # bias control: best replicate chosen separately for every cell (raw)
    ctl = []
    for (cl, ko), reps in allrep.items():
        row = [max(np.corrcoef(v, sS[c].values)[0, 1] for v in reps.values()) for c in cols]; ctl.append((ko, row))
    ctl_match = np.mean([r[ORDER.index(k)] for k, r in ctl]); ctl_other = np.mean([r[j] for k, r in ctl for j in range(len(ORDER)) if j != ORDER.index(k)])
    vmax = np.ceil(np.percentile(np.abs(np.concatenate([res[c][0].ravel() for c in res])), 99) * 10) / 10
    fig, axes = plt.subplots(1, 2, figsize=(17.5, 9.2), facecolor="white"); out = {}
    for ax, centre in zip(axes, [False, True]):
        C, lab = res[centre]
        im = ax.imshow(C, cmap=DIV, vmin=-vmax, vmax=vmax, aspect="auto", interpolation="nearest")
        ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=90, fontsize=9)
        ax.set_yticks(range(len(lab))); ax.set_yticklabels([f"{ko}  ({picks[(cl, ko)][0].split('_')[-1]})" for cl, ko in lab], fontsize=8.3)
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
            ax.add_patch(Rectangle((-2.75, a - 0.5), 0.35, b - a, color=LINE_COL[cl], clip_on=False))
            ax.text(-2.95, (a + b - 1) / 2, cl, rotation=90, ha="right", va="center", fontsize=9.5, fontweight="bold", color=LINE_COL[cl])
        for s in ax.spines.values(): s.set_visible(False)
        match = np.array([C[i, ORDER.index(ko)] for i, (_, ko) in enumerate(lab)]); mask = np.ones_like(C, bool)
        for i, (_, ko) in enumerate(lab): mask[i, ORDER.index(ko)] = False; mask[i, len(ORDER):] = False
        ranks = [int((C[i, :len(ORDER)] > C[i, ORDER.index(ko)]).sum() + 1) for i, (_, ko) in enumerate(lab)]
        ttl = "raw (best replicate − NTC)" if not centre else "centred across KOs (shared NTC / global component removed)"
        ax.set_title(f"{ttl}\nmatched r = {match.mean():+.2f} vs other-KO r = {C[mask].mean():+.2f};  matched ranks 1st in {sum(r == 1 for r in ranks)}/{len(ranks)}, median {np.median(ranks):.1f} of 8", fontsize=9.8, color=INK, pad=48, loc="left")
        out[centre] = (match.mean(), C[mask].mean(), sum(r == 1 for r in ranks), len(ranks), np.median(ranks))
    cax = fig.add_axes([0.36, 0.055, 0.3, 0.012]); cb = fig.colorbar(im, cax=cax, orientation="horizontal"); cb.set_label("Pearson r between bulk KO profile (best replicate) and CROP-seq KO profile", fontsize=9)
    fig.suptitle(title, x=0.02, ha="left", fontsize=14.5, fontweight="bold", y=0.985)
    fig.text(0.02, 0.945, "rows = bulk KO (RUN1), ONLY the replicate best correlated with the matching CROP-seq KO (replicate shown in brackets); columns = CROP-seq KO in the same order (matched pair outlined); MEF2A/MEF2D screen-only.", fontsize=8.8, color=INK2)
    fig.text(0.02, 0.925, f"selection-bias control: if the best replicate is picked separately for every cell, matched r = {ctl_match:+.2f} vs other-KO r = {ctl_other:+.2f} (the gap that selection alone produces).", fontsize=8.8, color=INK2)
    fig.subplots_adjust(left=0.10, right=0.99, top=0.80, bottom=0.11, wspace=0.42)
    fig.savefig(f"{FIG}/{fname}.png", dpi=190, facecolor="white"); fig.savefig(f"{FIG}/{fname}.pdf", facecolor="white"); plt.close(fig)
    pd.DataFrame([(cl, ko, picks[(cl, ko)][0], *[picks[(cl, ko)][1][s] for s in [picks[(cl, ko)][0]]]) for cl, ko in res[False][1]], columns=["line", "KO", "best_replicate", "r_to_matched_screen_KO"]).to_csv(f"{FIG}/{fname}_selected_replicates.csv", index=False)
    print(name, {("centred" if k else "raw"): tuple(round(float(x), 3) for x in v) for k, v in out.items()}, "| control matched/other:", round(ctl_match, 3), round(ctl_other, 3))

draw("modules", "Bulk KO (best replicate) vs CROP-seq KO similarity: 20 Cm modules", "fig5_bestreplicate_corr_modules_RUN1")
draw("variable_genes", f"Bulk KO (best replicate) vs CROP-seq KO similarity: screen's top variable genes (n={len(genes)})", "fig5_bestreplicate_corr_variable_genes_RUN1")
