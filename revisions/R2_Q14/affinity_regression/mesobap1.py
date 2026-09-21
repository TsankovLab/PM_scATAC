"""Shared helpers for reproducing TCGA MESO supplementary tab 2B (BAP1-associated inferred
TF activity, Hmeljak et al. 2018) with affinity regression."""
import os, csv, numpy as np
from scipy import stats
from affreg import ar_train, ar_model2w, ar_predict

HERE = os.path.dirname(os.path.abspath(__file__))
XLSX = os.path.join(HERE, "..", "21598290cd180804-sup-205173_2_supp_5073240_pg14sl.xlsx")
VARIANTS = [(r, u) for r in ("meso_tcga", "meso_tcga_pan_can_atlas_2018") for u in ("top5000", "motifhit")]
FOCUS = ("IRF8", "EGR2", "YY1")

def paper_2B():
    f = os.path.join(HERE, "prepared", "paper_2B.csv")
    if not os.path.exists(f):
        import openpyxl; wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
        rows = list(wb["2B_bap1_inferred_TF_activity"].iter_rows(values_only=True))[1:]
        with open(f, "w", newline="") as h:
            w = csv.writer(h); w.writerow(["TF", "p_tf", "p_adj_tf", "estimate_tf"])
            for r in rows:
                if r and r[0]: w.writerow(r[:4])
    return {r["TF"]: (float(r["estimate_tf"]), float(r["p_tf"]), float(r["p_adj_tf"]))
            for r in csv.DictReader(open(f))}

def load_variant(rel, univ):
    z = np.load(os.path.join(HERE, "prepared", f"input_{rel}_{univ}.npz"), allow_pickle=True)
    return {k: z[k] for k in z.files}

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out

def fit_full(V, lam, rsL2, sA, sB):
    m = ar_train(V["D"], V["P"], V["Y"], lam, rsL2, sA, sB)
    W = ar_model2w(m)
    return m, W, W @ V["P"].T                      # TF activity: TFs x samples

def ttest(act, tfs, status, grp=("inactivated",), ref=("no_inactivation",)):
    g, r = np.isin(status, grp), np.isin(status, ref)
    t, p = stats.ttest_ind(act[:, g], act[:, r], axis=1, equal_var=False)
    return dict(TF=[str(x) for x in tfs], est=act[:, g].mean(1) - act[:, r].mean(1),
                t=t, p=p, padj=bh(p), n=(int(g.sum()), int(r.sum())))

def compare(res, paper, thr=0.01):
    tfs = [t for t in res["TF"] if t in paper]
    ix = [res["TF"].index(t) for t in tfs]
    oe, op, oa = res["est"][ix], res["p"][ix], res["padj"][ix]
    pe = np.array([paper[t][0] for t in tfs]); pp = np.array([paper[t][1] for t in tfs])
    pa = np.array([paper[t][2] for t in tfs])
    same = np.sign(oe) == np.sign(pe); ps, os_ = pa < thr, oa < thr
    out = dict(n_shared=len(tfs), rho_est=stats.spearmanr(oe, pe)[0],
               rho_signed_logp=stats.spearmanr(np.sign(oe) * -np.log10(op), np.sign(pe) * -np.log10(pp))[0],
               pct_same_sign=100 * same.mean(), paper_sig=int(ps.sum()), our_sig=int(os_.sum()),
               both_sig_same_sign=int((ps & os_ & same).sum()),
               pct_same_sign_paper_sig=100 * same[ps].mean() if ps.any() else np.nan)
    order = np.array(tfs)[np.argsort(op)]
    for t in FOCUS:
        k = tfs.index(t)
        out[f"{t}_est"], out[f"{t}_p"], out[f"{t}_rank"] = oe[k], op[k], int(np.where(order == t)[0][0]) + 1
    return out

def cv(V, lam, rsL2, sA, sB, k=10, seed=1):
    """10-fold CV as in Osmanbeyoglu et al.: train on k-1 folds (re-centred on the training
    mean), reconstruct held-out expression from held-out RPPA, score per-sample correlation."""
    Y, D, P = V["Y"], V["D"], V["P"]; n = Y.shape[1]
    perm = np.random.default_rng(seed).permutation(n)
    edges = np.round(np.linspace(0, n, k + 1)).astype(int)
    pr, sp = [], []
    for f in range(k):
        te = perm[edges[f]:edges[f + 1]]; tr = np.setdiff1d(np.arange(n), te)
        mu, pm = Y[:, tr].mean(1, keepdims=True), P[tr].mean(0, keepdims=True)
        Ytr, Yte, Ptr, Pte = Y[:, tr] - mu, Y[:, te] - mu, P[tr] - pm, P[te] - pm
        rec = ar_predict(D, Pte, Ytr, ar_train(D, Ptr, Ytr, lam, rsL2, sA, sB))["rec"]
        for j in range(len(te)):
            pr.append(np.corrcoef(Yte[:, j], rec[:, j])[0, 1]); sp.append(stats.spearmanr(Yte[:, j], rec[:, j])[0])
    return float(np.mean(pr)), float(np.mean(sp))
