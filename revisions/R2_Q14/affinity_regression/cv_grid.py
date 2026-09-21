"""Hyperparameter selection exactly as Osmanbeyoglu et al.: 10-fold CV on expression
reconstruction.  Agreement with tab 2B is recorded per grid point for transparency but is
NOT used to select."""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[v] = "1"
import csv, itertools, time, numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from mesobap1 import *

LAMBDAS = [1e-4, 1e-3, 1e-2, 1e-1, 1.0]
RSL2 = [0.0, 0.1]
SPECS = [(1.0, 1.0), (0.95, 0.95), (0.9, 0.9), (0.8, 0.8), (1.0, 0.9), (0.9, 1.0)]

def task(args):
    rel, univ, lam, rs, sA, sB = args
    V = load_variant(rel, univ); paper = paper_2B(); t0 = time.time()
    try:
        cvp, cvs = cv(V, lam, rs, sA, sB)
        m, W, act = fit_full(V, lam, rs, sA, sB)
        c = compare(ttest(act, V["tfs"], V["bap1_status"]), paper)
        return dict(release=rel, universe=univ, lam=lam, rsL2=rs, specA=sA, specB=sB,
                    da=m["da"], db=m["db"], n_iter=m["n_iter"], cv_pearson=cvp, cv_spearman=cvs,
                    secs=round(time.time() - t0, 1), **c)
    except Exception as e:
        return dict(release=rel, universe=univ, lam=lam, rsL2=rs, specA=sA, specB=sB, error=repr(e))

if __name__ == "__main__":
    grid = [(r, u, l, rs, a, b) for (r, u) in VARIANTS for l in LAMBDAS for rs in RSL2 for (a, b) in SPECS]
    out = os.path.join(HERE, "prepared", "cv_grid_results.csv"); rows = []
    print(f"{len(grid)} grid points x 10 folds", flush=True)
    with ProcessPoolExecutor(max_workers=int(os.environ.get("LSB_DJOB_NUMPROC", 8))) as ex:
        for i, fut in enumerate(as_completed([ex.submit(task, g) for g in grid]), 1):
            rows.append(fut.result())
            keys = sorted({k for r in rows for k in r})
            with open(out, "w", newline="") as h:
                w = csv.DictWriter(h, fieldnames=keys); w.writeheader(); w.writerows(rows)
            r = rows[-1]
            print(f"[{i}/{len(grid)}] {r.get('release')} {r.get('universe')} lam={r.get('lam')} rs={r.get('rsL2')} "
                  f"spec={r.get('specA')}/{r.get('specB')} cv_r={r.get('cv_pearson', float('nan')):.3f} "
                  f"rho_slogp={r.get('rho_signed_logp', float('nan')):+.3f} {r.get('error','')}", flush=True)
