import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repro_common import *
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr

scr = pd.read_csv(f"{FIG}/screen_cellcycle_per_cell.csv", index_col=0)
bulk = pd.read_csv(f"{FIG}/bulk_cellcycle_per_sample_RUN1.csv")

# ---- screen: Wilcoxon vs NTC, BH-FDR across KOs ----
ntc = scr.loc[scr.merged_call == "NTC", "cc"]
rows = []
for g, d in scr.groupby("merged_call"):
    if g == "NTC": continue
    rows.append([g, len(d), d.cc.mean(), d.cc.mean() - ntc.mean(), mannwhitneyu(d.cc, ntc, alternative="two-sided").pvalue])
st = pd.DataFrame(rows, columns=["KO", "n", "mean_cc", "delta_vs_NTC", "p"])
o = st.p.rank(method="first").astype(int); st["padj"] = np.minimum.accumulate((st.p * len(st) / o).sort_values(ascending=False).reindex(st.index)[::-1][::-1]) if False else None
p = st.p.values; idx = np.argsort(p); adj = np.empty(len(p)); prev = 1.0
for r, i in zip(range(len(p), 0, -1), idx[::-1]): prev = min(prev, p[i] * len(p) / r); adj[i] = prev
st["padj"] = adj
st["stars"] = np.where(st.padj < 0.001, "***", np.where(st.padj < 0.01, "**", np.where(st.padj < 0.05, "*", "")))
st.to_csv(f"{FIG}/fig2_screen_cc_stats.csv", index=False)
print(st.round(4).to_string(index=False))

# ---------------- Figure 2A: distribution panels ----------------
fig, axes = plt.subplots(1, 5, figsize=(21, 4.6), gridspec_kw={"width_ratios": [1.5, 1.1, 1.1, 1.1, 1.0]}, facecolor="white")
ax = axes[0]
order = scr.groupby("merged_call").cc.mean().sort_values(ascending=False).index.tolist()
data = [scr.loc[scr.merged_call == g, "cc"].values for g in order]
vp = ax.violinplot(data, positions=range(len(order)), widths=0.9, showextrema=False)
for b, g in zip(vp["bodies"], order):
    b.set_facecolor("#a0a0a0" if g == "NTC" else "#8b2e2e"); b.set_alpha(0.55); b.set_edgecolor("none")
ax.boxplot(data, positions=range(len(order)), widths=0.12, showfliers=False, medianprops=dict(color=INK, lw=1.2),
           boxprops=dict(color=INK2, lw=0.8), whiskerprops=dict(color=INK2, lw=0.8), capprops=dict(color=INK2, lw=0.8))
ytop = max(np.percentile(d, 99.5) for d in data)
for i, g in enumerate(order):
    if g == "NTC": continue
    s_ = st.loc[st.KO == g, "stars"].iloc[0]
    if s_: ax.text(i, ytop + 0.05, s_, ha="center", fontsize=10, color=INK)
ax.set_xticks(range(len(order))); ax.set_xticklabels([f"{g}\n(n={len(scr[scr.merged_call==g])})" for g in order], rotation=90, fontsize=7.5)
ax.set_ylabel("cell-cycle index (S.Score + G2M.Score)", fontsize=9.5)
ax.set_title("CROP-seq screen (single cells)", loc="left", fontsize=10.5, fontweight="bold", color=LINE_COL["CROP-seq"])
ax.axhline(ntc.mean(), color=MUTED, ls=(0, (4, 3)), lw=0.9)

bstats = []
for ax, cl in zip(axes[1:], ["MSTO-211H", "NCI-H2052", "NCI-H2452", "NCI-H28"]):
    d = bulk[bulk.Cell_Line == cl]
    order = d.groupby("Gene_Control").cc.mean().sort_values(ascending=False).index.tolist()
    nm = d.loc[d.Gene_Control == "NTC", "cc"].mean()
    for i, g in enumerate(order):
        v = d.loc[d.Gene_Control == g, "cc"].values
        ax.bar(i, v.mean(), width=0.7, color="#a0a0a0" if g == "NTC" else LINE_COL[cl], alpha=0.35, zorder=1)
        ax.scatter(np.full(len(v), i) + np.linspace(-0.12, 0.12, len(v)) * (len(v) > 1), v, s=22, color=INK, zorder=3, edgecolor="white", lw=0.5)
        bstats.append([cl, g, len(v), v.mean(), v.mean() - nm])
    ax.axhline(nm, color=MUTED, ls=(0, (4, 3)), lw=0.9)
    ax.set_xticks(range(len(order))); ax.set_xticklabels([f"{g}\n(n={len(d[d.Gene_Control==g])})" for g in order], rotation=90, fontsize=7.5)
    ax.set_title(f"{cl} bulk (RUN1)", loc="left", fontsize=10.5, fontweight="bold", color=LINE_COL[cl])
    ax.axhline(0, color=GRID, lw=0.8, zorder=0)
for a in axes:
    for s in ["top", "right"]: a.spines[s].set_visible(False)
    a.spines["left"].set_color(GRID); a.spines["bottom"].set_color(GRID); a.tick_params(colors=INK2)
fig.suptitle("Cell-cycle index by KO relative to NTC (sorted by mean; dashed = NTC mean)", x=0.01, ha="left", fontsize=14, fontweight="bold", y=1.04)
plt.tight_layout()
fig.savefig(f"{FIG}/fig2A_cell_cycle_index_by_KO.png", dpi=200, bbox_inches="tight", facecolor="white")
fig.savefig(f"{FIG}/fig2A_cell_cycle_index_by_KO.pdf", bbox_inches="tight", facecolor="white")

# ---------------- Figure 2B: bulk delta vs screen delta ----------------
bs = pd.DataFrame(bstats, columns=["cell_line", "KO", "n", "mean_cc", "delta_vs_NTC"]); bs = bs[bs.KO != "NTC"]
bs.to_csv(f"{FIG}/fig2_bulk_cc_deltas.csv", index=False)
m = bs.merge(st[["KO", "delta_vs_NTC", "padj"]].rename(columns={"delta_vs_NTC": "screen_delta", "padj": "screen_padj"}), on="KO")
fig, axes = plt.subplots(1, 5, figsize=(21, 4.4), facecolor="white")
summ = []
for ax, cl in zip(axes[:4], ["MSTO-211H", "NCI-H2052", "NCI-H2452", "NCI-H28"]):
    d = m[m.cell_line == cl]
    ax.axhline(0, color=GRID, lw=0.8); ax.axvline(0, color=GRID, lw=0.8)
    ax.scatter(d.screen_delta, d.delta_vs_NTC, s=60, color=LINE_COL[cl], edgecolor="white", zorder=3)
    for _, r in d.iterrows(): ax.annotate(r.KO, (r.screen_delta, r.delta_vs_NTC), xytext=(4, 4), textcoords="offset points", fontsize=8, color=INK)
    rho, pv = spearmanr(d.screen_delta, d.delta_vs_NTC); agree = (np.sign(d.screen_delta) == np.sign(d.delta_vs_NTC)).mean()
    summ.append([cl, len(d), rho, pv, agree])
    ax.set_title(f"{cl}\nSpearman ρ={rho:+.2f} (p={pv:.2f}) · sign agreement {agree:.0%}", loc="left", fontsize=9, fontweight="bold", color=LINE_COL[cl])
    ax.set_xlabel("screen Δ cell-cycle index (KO − NTC)", fontsize=9); ax.set_ylabel("bulk Δ cell-cycle index (KO − NTC)", fontsize=9)
ax = axes[4]
ax.axhline(0, color=GRID, lw=0.8); ax.axvline(0, color=GRID, lw=0.8)
for cl in ["MSTO-211H", "NCI-H2052", "NCI-H2452", "NCI-H28"]:
    d = m[m.cell_line == cl]; ax.scatter(d.screen_delta, d.delta_vs_NTC, s=40, color=LINE_COL[cl], edgecolor="white", label=cl, zorder=3)
rho, pv = spearmanr(m.screen_delta, m.delta_vs_NTC); summ.append(["ALL (pooled)", len(m), rho, pv, (np.sign(m.screen_delta) == np.sign(m.delta_vs_NTC)).mean()])
ax.set_title(f"all lines pooled\nρ={rho:+.2f} (p={pv:.2f})", loc="left", fontsize=9, fontweight="bold", color=INK)
ax.legend(frameon=False, fontsize=8); ax.set_xlabel("screen Δ cell-cycle index (KO − NTC)", fontsize=9)
for a in axes:
    for s in ["top", "right"]: a.spines[s].set_visible(False)
    a.tick_params(colors=INK2)
fig.suptitle("Does bulk reproduce the screen's proliferation shifts? (shared KOs only)", x=0.01, ha="left", fontsize=14, fontweight="bold", y=1.04)
plt.tight_layout()
fig.savefig(f"{FIG}/fig2B_cell_cycle_bulk_vs_screen.png", dpi=200, bbox_inches="tight", facecolor="white")
fig.savefig(f"{FIG}/fig2B_cell_cycle_bulk_vs_screen.pdf", bbox_inches="tight", facecolor="white")
pd.DataFrame(summ, columns=["cell_line", "n_KO", "spearman_rho", "p", "sign_agreement"]).to_csv(f"{FIG}/fig2B_summary.csv", index=False)
print(pd.DataFrame(summ, columns=["cell_line", "n_KO", "rho", "p", "sign_agree"]).round(3).to_string(index=False))
