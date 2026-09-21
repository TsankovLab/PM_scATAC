import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from matplotlib.colors import LinearSegmentedColormap

PROC = "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CROPseq_validation_bulkRNA/processed"
PP   = f"{PROC}/de/perturbation_panel"
FIG  = f"{PP}/figs_repro"
LINES = ["NCI-H2052", "MSTO-211H", "NCI-H2452", "NCI-H28"]
LINE_COL = {"CROP-seq": "#4a3aa7", "NCI-H2052": "#eda100", "MSTO-211H": "#1baf7a", "NCI-H2452": "#e87ba4", "NCI-H28": "#008300"}
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
DIV = LinearSegmentedColormap.from_list("blue_orange", ["#2a78d6", "#f0efec", "#eb6834"])

def load_bulk():
    counts = pd.read_csv(f"{PROC}/counts/gene_counts_matrix_stranded_reverse.tsv", sep="\t", index_col=0)
    sym = pd.read_csv(f"{PROC}/counts/gene_id_to_symbol.tsv", sep="\t", header=None, names=["gene_id", "symbol"])
    sym = sym[sym["gene_id"].isin(counts.index) & sym["symbol"].notna() & (sym["symbol"] != "")]
    cs = counts.loc[sym["gene_id"]].groupby(sym["symbol"].values).sum()
    meta = pd.read_csv(f"{PROC}/sample_sheet/sample_metadata_RUN1.tsv", sep="\t")
    meta = meta[(meta["status"] == "OK") & meta["Cell_Line"].isin(LINES)].set_index("fastq_sample_id")
    meta = meta[[s in cs.columns for s in meta.index]]
    cp10k = cs.div(cs.sum(axis=0), axis=1) * 1e4
    return cp10k[meta.index], meta

def group_log2avg(cp10k, meta, cl):
    """log2(mean CP10K + 1) per KO group within a cell line (genes x KO)."""
    m = meta[meta["Cell_Line"] == cl]
    return pd.DataFrame({g: np.log2(cp10k[idx].mean(axis=1) + 1)
                         for g, idx in m.groupby("Gene_Control").groups.items()})

def zscore_rows_of(mat):
    """R's scale(t(x)): z-score each gene (column of t(x)) across KO groups; NaN->0."""
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1, ddof=1), axis=0)
    return z.fillna(0.0)

def corr_order(x, axis=0):
    """Hierarchical order (complete linkage, 1-Pearson) of rows of x (axis=0) or columns (axis=1)."""
    a = x.values if axis == 0 else x.values.T
    if a.shape[0] < 3:
        return np.arange(a.shape[0])
    return leaves_list(linkage(pdist(a, "correlation"), "complete"))
