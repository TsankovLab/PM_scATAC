import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repro_common import *
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

GREY, BROWN = "#BEBEBE", "#A52A2A"          # R's 'grey' / 'brown', as in CRISPR_cropseq.R
scr = pd.read_csv(f"{FIG}/screen_cellcycle_per_cell.csv", index_col=0)
st = pd.read_csv(f"{FIG}/fig2_screen_cc_stats.csv").set_index("KO")
bulk = pd.read_csv(f"{FIG}/bulk_cellcycle_per_sample_RUN1.csv")
THREE = ["MSTO-211H", "NCI-H2052", "NCI-H2452"]

def style_box(bp, colors):
    for patch, c in zip(bp["boxes"], colors): patch.set_facecolor(c); patch.set_alpha(0.6); patch.set_edgecolor(INK2); patch.set_linewidth(0.8)
    for k in ["whiskers", "caps"]:
        for l in bp[k]: l.set_color(INK2); l.set_linewidth(0.8)
    for l in bp["medians"]: l.set_color(INK); l.set_linewidth(1.4)

fig, axes = plt.subplots(1, 4, figsize=(21, 5.2), gridspec_kw={"width_ratios": [1.45, 1.15, 1.15, 1.05]}, facecolor="white")

# ---- CROP-seq: violin + box, sorted by mean (as original) ----
ax = axes[0]
order = scr.groupby("merged_call").cc.mean().sort_values(ascending=False).index.tolist()
data = [scr.loc[scr.merged_call == g, "cc"].values for g in order]; cols = [GREY if g == "NTC" else BROWN for g in order]
vp = ax.violinplot(data, positions=range(len(order)), widths=0.95, showextrema=False)
for b, c in zip(vp["bodies"], cols): b.set_facecolor(c); b.set_alpha(0.6); b.set_edgecolor("none")
style_box(ax.boxplot(data, positions=range(len(order)), widths=0.12, showfliers=False, patch_artist=True), ["white"] * len(order))
ytop = max(np.percentile(d, 99.7) for d in data)
for i, g in enumerate(order):
    if g != "NTC" and st.loc[g, "stars"] == st.loc[g, "stars"]: ax.text(i, ytop + 0.04, st.loc[g, "stars"], ha="center", fontsize=10)
ax.set_xticks(range(len(order))); ax.set_xticklabels([f"{g}\n(n={len(d)})" for g, d in zip(order, data)], rotation=90, fontsize=8)
ax.set_title("CROP-seq screen (cells; Wilcoxon vs NTC, FDR)", loc="left", fontsize=10.5, fontweight="bold", color=INK)

# ---- bulk: box + every replicate as a point ----
rng = np.random.default_rng(1)
for ax, cl in zip(axes[1:], THREE):
    d = bulk[bulk.Cell_Line == cl]
    order = d.groupby("Gene_Control").cc.mean().sort_values(ascending=False).index.tolist()
    data = [d.loc[d.Gene_Control == g, "cc"].values for g in order]; cols = [GREY if g == "NTC" else BROWN for g in order]
    style_box(ax.boxplot(data, positions=range(len(order)), widths=0.55, showfliers=False, patch_artist=True), cols)
    for i, v in enumerate(data):
        ax.scatter(i + rng.uniform(-0.13, 0.13, len(v)), v, s=26, color=INK, zorder=4, edgecolor="white", linewidth=0.6)
    ax.axhline(d.loc[d.Gene_Control == "NTC", "cc"].mean(), color=MUTED, ls=(0, (4, 3)), lw=0.9, zorder=0)
    ax.set_xticks(range(len(order))); ax.set_xticklabels([f"{g}\n(n={len(v)})" for g, v in zip(order, data)], rotation=90, fontsize=8)
    ax.set_title(f"{cl} bulk (RUN1; each dot = one replicate)", loc="left", fontsize=10.5, fontweight="bold", color=INK)
for a in axes:
    for s in ["top", "right"]: a.spines[s].set_visible(False)
    a.spines["left"].set_color(GRID); a.spines["bottom"].set_color(GRID); a.tick_params(colors=INK2)
    a.set_ylabel("cell-cycle index (S.Score + G2M.Score)", fontsize=9.5)
    a.yaxis.grid(True, color=GRID, lw=0.7); a.set_axisbelow(True)
fig.suptitle("Proliferation index by KO: CROP-seq vs bulk RNA-seq (S + G2M cell-cycle marker score, sorted by mean; grey = NTC, brown = KO)", x=0.01, ha="left", fontsize=13.5, fontweight="bold", y=1.03)
plt.tight_layout()
fig.savefig(f"{FIG}/fig2C_proliferation_boxplots_CROPseq_vs_bulk_RUN1.png", dpi=200, bbox_inches="tight", facecolor="white")
fig.savefig(f"{FIG}/fig2C_proliferation_boxplots_CROPseq_vs_bulk_RUN1.pdf", bbox_inches="tight", facecolor="white")
print("saved")
