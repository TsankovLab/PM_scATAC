"""Second CV grid: preprocessing choices the methods leave unstated (per-sample unit-norm Y,
column-normalised D, training on all 63 vs the 54 labelled tumours), crossed with a finer
lambda/spectrum grid.  Selection criterion is 10-fold CV reconstruction of held-out
expression, as in Osmanbeyoglu et al.; agreement with tab 2B is recorded, not used."""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[v] = "1"
import csv, itertools, time, numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy import stats
from affreg import ar_train, ar_predict, ar_model2w
from mesobap1 import HERE, load_variant, paper_2B, ttest, compare

def prepare(V, train, Dnorm, Ynorm):
    Y, D, P, st = V["Y"], V["D"].astype(float), V["P"], V["bap1_status"]
    if train == "labelled54":
        k = st != "NA"; Y, P, st = Y[:, k], P[k], st[k]
        Y = Y - Y.mean(1, keepdims=True); P = P - P.mean(0, keepdims=True)
    if Dnorm == "colnorm": D = D / np.sqrt((D ** 2).sum(0, keepdims=True))
    return Y, D, P, st

def ynorm(Y, how):
    return Y / np.sqrt((Y ** 2).sum(0, keepdims=True)) if how == "colnorm" else Y

def cv(Y, D, P, how, lam, sA, sB, k=10, seed=1):
    n = Y.shape[1]; perm = np.random.default_rng(seed).permutation(n)
    edges = np.round(np.linspace(0, n, k + 1)).astype(int); pr = []
    for f in range(k):
        te = perm[edges[f]:edges[f + 1]]; tr = np.setdiff1d(np.arange(n), te)
        mu, pm = Y[:, tr].mean(1, keepdims=True), P[tr].mean(0, keepdims=True)
        Ytr, Yte = ynorm(Y[:, tr] - mu, how), ynorm(Y[:, te] - mu, how)
        rec = ar_predict(D, P[te] - pm, Ytr, ar_train(D, P[tr] - pm, Ytr, lam, 0.0, sA, sB))["rec"]
        pr += [np.corrcoef(Yte[:, j], rec[:, j])[0, 1] for j in range(len(te))]
    return float(np.mean(pr))

def task(a):
    rel, univ, train, Dn, Yn, lam, sA, sB = a
    V = load_variant(rel, univ); paper = paper_2B(); t0 = time.time()
    base = dict(release=rel, universe=univ, train=train, Dnorm=Dn, Ynorm=Yn, lam=lam, specA=sA, specB=sB)
    try:
        Y, D, P, st = prepare(V, train, Dn, Yn)
        cvp = cv(Y, D, P, Yn, lam, sA, sB)
        m = ar_train(D, P, ynorm(Y, Yn), lam, 0.0, sA, sB)
        c = compare(ttest(ar_model2w(m) @ P.T, V["tfs"], st), paper)
        return dict(base, da=m["da"], db=m["db"], n_iter=m["n_iter"], cv_pearson=cvp, secs=round(time.time() - t0, 1), **c)
    except Exception as e:
        return dict(base, error=repr(e))

if __name__ == "__main__":
    grid = list(itertools.product(["meso_tcga_pan_can_atlas_2018", "meso_tcga"], ["top5000"],
                                  ["all63", "labelled54"], ["binary", "colnorm"], ["none", "colnorm"],
                                  [1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0],
                                  [(1.0, 1.0), (0.95, 0.95), (0.95, 0.9), (0.9, 0.9), (0.8, 0.8)]))
    grid = [g[:6] + g[6] for g in grid]
    out = os.path.join(HERE, "prepared", "cv_grid2_results.csv"); rows = []
    print(f"{len(grid)} grid points x 10 folds", flush=True)
    with ProcessPoolExecutor(max_workers=int(os.environ.get("LSB_DJOB_NUMPROC", 8))) as ex:
        for i, fut in enumerate(as_completed([ex.submit(task, g) for g in grid]), 1):
            rows.append(fut.result())
            if i % 20 == 0 or i == len(grid):
                keys = sorted({k for r in rows for k in r})
                with open(out, "w", newline="") as h:
                    w = csv.DictWriter(h, fieldnames=keys); w.writeheader(); w.writerows(rows)
                print(f"[{i}/{len(grid)}]", flush=True)
