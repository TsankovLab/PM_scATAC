import sys
import os
import glob
import pickle
import pandas as pd
import seaborn as sns
from dask.distributed import Client
from dask.diagnostics import ProgressBar

from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase


def main():
    # -------------------------------
    # Command-line arguments
    # -------------------------------
    projdir_SC = sys.argv[1]
    motifs_tss = sys.argv[2]
    motifs_weights = sys.argv[3]
    TFs = sys.argv[4]
    expr_mat = sys.argv[5]
    vg = int(sys.argv[6])            # numeric variable
    motif_window = sys.argv[7]       # string variable
    scheduler_file = sys.argv[8]   

    # -------------------------------
    # Connect to Dask cluster
    # -------------------------------
    client = Client(scheduler_file=scheduler_file)
    print("[DASK] Connected:", client)

    # -------------------------------
    # File paths
    # -------------------------------
    DATA_FOLDER = projdir_SC
    RESOURCES_FOLDER = "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/SCENIC_motifs/"
    DATABASE_FOLDER = RESOURCES_FOLDER
    DATABASES_GLOB = os.path.join(DATABASE_FOLDER, motifs_tss)
    MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, motifs_weights)
    MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, TFs)
    SC_EXP_FNAME = os.path.join(DATA_FOLDER, expr_mat)
    REGULONS_FNAME = os.path.join(DATA_FOLDER, "regulons.p")
    MOTIFS_FNAME = os.path.join(DATA_FOLDER, "motifs.csv")

    # -------------------------------
    # Load expression matrix and TFs
    # -------------------------------
    ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0)
    if ex_matrix.shape[0] < ex_matrix.shape[1]:
        ex_matrix = ex_matrix.T

    tf_names = load_tf_names(MM_TFS_FNAME)
    tf_names = [tf for tf in tf_names if tf in ex_matrix.columns]
    print("Number of TFs to use:", len(tf_names))

    # -------------------------------
    # Load motif databases
    # -------------------------------
    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]

    db_fnames = glob.glob(DATABASES_GLOB)
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

    # -------------------------------
    # Run GRNBoost2
    # -------------------------------
    adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)

    # -------------------------------
    # Build modules & regulons
    # -------------------------------
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
    with ProgressBar():
        df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)
    regulons = df2regulons(df)

    # Save results
    df.to_csv(MOTIFS_FNAME)
    with open(REGULONS_FNAME, "wb") as f:
        pickle.dump(regulons, f)

    # -------------------------------
    # AUC calculation & clustermap
    # -------------------------------
    auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
    clustermap = sns.clustermap(auc_mtx, figsize=(8, 8))
    clustermap.savefig(os.path.join(DATA_FOLDER, "clustermap.png"), format="png")
    auc_mtx.to_csv(os.path.join(DATA_FOLDER, "auc_mtx.csv"), index=True, header=True)


if __name__ == "__main__":
    main()
