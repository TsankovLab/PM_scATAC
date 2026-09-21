"""Sensitivity of the tab-2B reproduction to implementation details the methods do not
state.  One factor is varied at a time from a base configuration, at three lambdas so a
change of scale is not mistaken for a change of method.  Agreement with 2B is reported for
transparency; hyperparameters are selected separately by cross-validation (cv_grid.py)."""
import csv, itertools, numpy as np
from scipy import stats
from affreg import ar_train, ar_model2w
from mesobap1 import load_variant, paper_2B, ttest, compare

paper = paper_2B()
BASE = dict(release="meso_tcga_pan_can_atlas_2018", universe="top5000", train="all63", Dnorm="binary",
            keep_diag=False, Ynorm="none", test="welch", groups="inact_vs_no")
FACTORS = dict(train=["all63", "labelled54"], Dnorm=["binary", "colnorm"], keep_diag=[False, True],
               Ynorm=["none", "colnorm"], test=["welch", "student"], groups=["inact_vs_no", "inact+possibly_vs_no"],
               release=["meso_tcga_pan_can_atlas_2018", "meso_tcga"], universe=["top5000", "motifhit"])
LAMS, SPEC = [0.01, 0.1, 1.0], (0.95, 0.9)

def run(cfg, lam):
    V = load_variant(cfg["release"], cfg["universe"])
    Y, D, P, st, tfs = V["Y"], V["D"].astype(float), V["P"], V["bap1_status"], V["tfs"]
    if cfg["train"] == "labelled54":
        k = st != "NA"; Y, P, st = Y[:, k], P[k], st[k]
        Y = Y - Y.mean(1, keepdims=True); P = P - P.mean(0, keepdims=True)
    if cfg["Dnorm"] == "colnorm": D = D / np.sqrt((D ** 2).sum(0, keepdims=True))
    if cfg["Ynorm"] == "colnorm": Y = Y / np.sqrt((Y ** 2).sum(0, keepdims=True))
    m = ar_train(D, P, Y, lam, 0.0, SPEC[0], SPEC[1], old_version=cfg["keep_diag"])
    act = ar_model2w(m) @ P.T
    grp = ("inactivated", "possibly_inactivated") if cfg["groups"] != "inact_vs_no" else ("inactivated",)
    res = ttest(act, tfs, st, grp=grp)
    if cfg["test"] == "student":
        g, r = np.isin(st, grp), np.isin(st, ("no_inactivation",))
        t, p = stats.ttest_ind(act[:, g], act[:, r], axis=1, equal_var=True)
        from mesobap1 import bh
        res.update(t=t, p=p, padj=bh(p))
    return m, compare(res, paper)

rows = []
for factor, levels in FACTORS.items():
    for lev in levels:
        cfg = dict(BASE, **{factor: lev})
        for lam in LAMS:
            m, c = run(cfg, lam)
            rows.append(dict(factor=factor, level=str(lev), lam=lam, da=m["da"], db=m["db"], iters=m["n_iter"], **c))
            print(f"{factor:9s} {str(lev):22s} lam={lam:<5} da={m['da']:2d} db={m['db']:2d} | rho_est={c['rho_est']:+.3f} "
                  f"rho_slogp={c['rho_signed_logp']:+.3f} same_sign={c['pct_same_sign']:.0f}% "
                  f"sig both/paper={c['both_sig_same_sign']}/{c['paper_sig']} | "
                  + " ".join(f"{t}:{c[t+'_est']:+.2g}(#{c[t+'_rank']})" for t in ("IRF8", "EGR2", "YY1")), flush=True)
with open("prepared/sensitivity.csv", "w", newline="") as h:
    w = csv.DictWriter(h, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
