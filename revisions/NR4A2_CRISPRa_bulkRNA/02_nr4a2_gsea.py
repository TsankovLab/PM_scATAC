"""02_nr4a2_gsea.py: preranked GSEA of the NR4A2 enhancer CRISPRa vs SCRAMBLE DESeq2 statistic (01_nr4a2_deseq2.R).
ranked list: DESeq2_results_full.csv (padj order) -> genes with a symbol -> sorted by |stat| -> first row per symbol (larger |stat|) -> sorted by stat,
decreasing (pandas default sorts; this reproduces the tie order of the original ranked_gene_list_stat.rnk, which affects the permutation results).
gseapy prerank (gseapy 1.3.1) with the Enrichr libraries MSigDB_Hallmark_2020 and GO_Biological_Process_2023 (downloaded by gseapy; needs internet),
1,000 permutations, seed 42, gene sets of 10-1,000 genes. Reconstruction of the original run; verified to reproduce the saved
gsea_MSigDB_Hallmark_2020.csv and gsea_GO_Biological_Process_2023.csv (ES, NES, p, FDR, leading edge).
outputs (DE_DIR/gsea): ranked_gene_list_stat.rnk, gsea_MSigDB_Hallmark_2020.csv, gsea_GO_Biological_Process_2023.csv
Usage: [DE_DIR=...] python3 02_nr4a2_gsea.py
"""
import os
import pandas as pd
import gseapy as gp

DE_DIR = os.environ.get("DE_DIR", "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CROPseq_validation_bulkRNA/processed/de/NR4A2_1_2_vs_SCRAMBLE")
os.makedirs(os.path.join(DE_DIR, "gsea"), exist_ok=True)
de = pd.read_csv(os.path.join(DE_DIR, "DESeq2_results_full.csv")).dropna(subset=["stat", "symbol"])
de = de[de.symbol != ""]
rnk = de.assign(a=de.stat.abs()).sort_values("a", ascending=False).drop_duplicates("symbol").set_index("symbol").stat.sort_values(ascending=False)
rnk.to_csv(os.path.join(DE_DIR, "gsea", "ranked_gene_list_stat.rnk"), sep="\t", header=False)
rnk = pd.read_csv(os.path.join(DE_DIR, "gsea", "ranked_gene_list_stat.rnk"), sep="\t", header=None, index_col=0).iloc[:, 0]   # GSEA on the list as written: the
print(f"ranked genes: {len(rnk)}")   # text round trip changes which values tie, and with it the permutation results (the original ran on the file)
for lib in ["MSigDB_Hallmark_2020", "GO_Biological_Process_2023"]:
    res = gp.prerank(rnk=rnk, gene_sets=lib, permutation_num=1000, min_size=10, max_size=1000, seed=42, threads=1, outdir=None, verbose=False).res2d
    res.to_csv(os.path.join(DE_DIR, "gsea", f"gsea_{lib}.csv"), index=False)
    print(f"{lib}: {len(res)} sets, FDR < 0.05: {(res['FDR q-val'].astype(float) < 0.05).sum()}")
