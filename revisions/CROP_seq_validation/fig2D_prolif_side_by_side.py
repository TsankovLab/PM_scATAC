import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repro_common import *
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

GREY, BROWN = "#BEBEBE", "#A52A2A"
scr = pd.read_csv(f"{FIG}/screen_cellcycle_per_cell.csv", index_col=0)
st = pd.read_csv(f"{FIG}/fig2_screen_cc_stats.csv").set_index("KO")
bulk = pd.read_csv(f"{FIG}/bulk_cellcycle_per_sample_RUN1.csv")
THREE = ["MSTO-211H", "NCI-H2052", "NCI-H2452"]
shared = set(scr.merged_call) - {"NTC"}
bulk = bulk[bulk.Gene_Control.isin(shared | {"NTC"})]
print("bulk KOs kept:", sorted(set(bulk.Gene_Control) - {"NTC"}), "| screen-only:", sorted(shared - set(bulk.Gene_Control)))

def style_box(bp, colors):
    for patch, c in zip(bp["boxes"], colors): patch.set_facecolor(c); patch.set_alpha(0.6); patch.set_edgecolor(INK2); patch.set_linewidth(0.8)
    for k in ["whiskers", "caps"]:
        for l in bp[k]: l.set_color(INK2); l.set_linewidth(0.8)
    for l in bp["medians"]: l.set_color(INK); l.set_linewidth(1.4)

def order_of(means):                       # NTC first, KOs sorted by mean (high -> low)
    return ["NTC"] + [g for g in means.sort_values(ascending=False).index if g != "NTC"]

so = order_of(scr.groupby("merged_call").cc.mean())
bo = {cl: [g for g in so if g in set(bulk[bulk.Cell_Line == cl].Gene_Control)] for cl in THREE}   # same order as the screen
widths = [len(so), 2.2] + [len(bo[cl]) for cl in THREE]
fig = plt.figure(figsize=(19, 5.6), facecolor="white")
gs = GridSpec(1, 5, width_ratios=widths, wspace=0.16, left=0.045, right=0.99, top=0.80, bottom=0.25)

ax = fig.add_subplot(gs[0, 0])
data = [scr.loc[scr.merged_call == g, "cc"].values for g in so]; cols = [GREY if g == "NTC" else BROWN for g in so]
vp = ax.violinplot(data, positions=range(len(so)), widths=0.95, showextrema=False)
for b, c in zip(vp["bodies"], cols): b.set_facecolor(c); b.set_alpha(0.6); b.set_edgecolor("none")
style_box(ax.boxplot(data, positions=range(len(so)), widths=0.12, showfliers=False, patch_artist=True), ["white"] * len(so))
ytop = max(np.percentile(d, 99.7) for d in data)
for i, g in enumerate(so):
    if g != "NTC" and st.loc[g, "stars"] == st.loc[g, "stars"]: ax.text(i, ytop + 0.04, st.loc[g, "stars"], ha="center", fontsize=10)
ax.set_xticks(range(len(so))); ax.set_xticklabels([f"{g} ({len(d)})" for g, d in zip(so, data)], rotation=90, fontsize=8.5)
ax.axvline(0.5, color="white", lw=3); ax.set_ylabel("cell-cycle index (S.Score + G2M.Score)", fontsize=9.5)
ax.text(0, 1.06, "CROP-seq screen", transform=ax.transAxes, fontsize=11, fontweight="bold", color=LINE_COL["CROP-seq"])
axs = [ax]

rng = np.random.default_rng(1)
for gi, cl in zip([2, 3, 4], THREE):
    ax = fig.add_subplot(gs[0, gi]); d = bulk[bulk.Cell_Line == cl]; order = bo[cl]
    data = [d.loc[d.Gene_Control == g, "cc"].values for g in order]; cols = [GREY if g == "NTC" else BROWN for g in order]
    style_box(ax.boxplot(data, positions=range(len(order)), widths=0.55, showfliers=False, patch_artist=True), cols)
    for i, v in enumerate(data): ax.scatter(i + rng.uniform(-0.13, 0.13, len(v)), v, s=26, color=INK, zorder=4, edgecolor="white", linewidth=0.6)
    ax.axhline(d.loc[d.Gene_Control == "NTC", "cc"].mean(), color=MUTED, ls=(0, (4, 3)), lw=0.9, zorder=0)
    ax.set_xticks(range(len(order))); ax.set_xticklabels([f"{g} ({len(v)})" for g, v in zip(order, data)], rotation=90, fontsize=8.5)
    ax.axvline(0.5, color="white", lw=3)
    ax.text(0, 1.06, f"{cl} bulk (RUN1)", transform=ax.transAxes, fontsize=11, fontweight="bold", color=LINE_COL[cl]); axs.append(ax)
for a in axs:
    for s in ["top", "right"]: a.spines[s].set_visible(False)
    a.spines["left"].set_color(GRID); a.spines["bottom"].set_color(GRID); a.tick_params(colors=INK2)
    a.yaxis.grid(True, color=GRID, lw=0.7); a.set_axisbelow(True); a.set_xlim(-0.6, len(a.get_xticks()) - 0.4)
fig.text(0.045, 0.965, "Proliferation index by KO: CROP-seq (left) vs bulk RNA-seq, three cell lines, RUN1 (right)", fontsize=14, fontweight="bold", color=INK)
fig.text(0.045, 0.925, "cell-cycle index = S.Score + G2M.Score (same as CROP-seq); bulk restricted to KOs present in the screen; grey = NTC, brown = KO; NTC first, then in the CROP-seq order (sorted by screen mean);", fontsize=8.8, color=INK2)
fig.text(0.045, 0.902, "bulk: each dot = one replicate, dashed = NTC mean, y-axis standardised per line; screen: violin + box, stars = Wilcoxon vs NTC (FDR); bracket = cells (screen) or replicates (bulk).", fontsize=8.8, color=INK2)
fig.savefig(f"{FIG}/fig2D_proliferation_side_by_side_CROPseq_vs_3lines_RUN1.png", dpi=200, facecolor="white")
fig.savefig(f"{FIG}/fig2D_proliferation_side_by_side_CROPseq_vs_3lines_RUN1.pdf", facecolor="white")
print("saved")
