"""Is the agreement with tab 2B driven by the BAP1 labels?  Fit one configuration once, then
permute the BAP1 labels among the labelled tumours and recompute the t-test and agreement
statistics on the SAME activity matrix.  Writes results/<tag>/permutation_null.json."""
import argparse, json, os, numpy as np
from affreg import ar_train, ar_model2w
from mesobap1 import HERE, load_variant, paper_2B, ttest, compare

ap = argparse.ArgumentParser()
ap.add_argument("--release", default="meso_tcga_pan_can_atlas_2018"); ap.add_argument("--universe", default="top5000")
ap.add_argument("--train", default="all63"); ap.add_argument("--dnorm", default="binary"); ap.add_argument("--ynorm", default="colnorm")
ap.add_argument("--lam", type=float, default=0.1); ap.add_argument("--specA", type=float, default=0.95)
ap.add_argument("--specB", type=float, default=0.9); ap.add_argument("--nperm", type=int, default=2000)
ap.add_argument("--tag", default="provisional"); a = ap.parse_args()

paper = paper_2B(); rng = np.random.default_rng(2018)
V = load_variant(a.release, a.universe)
Y, D, P, st = V["Y"], V["D"].astype(float), V["P"], V["bap1_status"]
if a.train == "labelled54":
    k = st != "NA"; Y, P, st = Y[:, k], P[k], st[k]
    Y = Y - Y.mean(1, keepdims=True); P = P - P.mean(0, keepdims=True)
if a.dnorm == "colnorm": D = D / np.sqrt((D ** 2).sum(0, keepdims=True))
if a.ynorm == "colnorm": Y = Y / np.sqrt((Y ** 2).sum(0, keepdims=True))
act = ar_model2w(ar_train(D, P, Y, a.lam, 0.0, a.specA, a.specB)) @ P.T
obs = compare(ttest(act, V["tfs"], st), paper)
lab = np.where(np.isin(st, ["inactivated", "no_inactivation"]))[0]
KEYS = ("rho_est", "rho_signed_logp", "both_sig_same_sign", "pct_same_sign")
null = {k: [] for k in KEYS}
for _ in range(a.nperm):
    s = st.copy(); s[lab] = rng.permutation(st[lab])
    c = compare(ttest(act, V["tfs"], s), paper)
    for k in KEYS: null[k].append(c[k])
res = {}
print(f"[{a.tag}] {a.nperm} label permutations")
for k in KEYS:
    v = np.array(null[k], float)
    res[k] = dict(observed=float(obs[k]), null_mean=float(np.nanmean(v)), null_sd=float(np.nanstd(v)),
                  null_95=float(np.nanpercentile(v, 95)), null_max=float(np.nanmax(v)),
                  p=float((1 + np.sum(v >= obs[k])) / (a.nperm + 1)))
    r = res[k]; print(f"  {k:20s} observed {r['observed']:8.3f} | null mean {r['null_mean']:7.3f} sd {r['null_sd']:6.3f} "
                      f"95th {r['null_95']:7.3f} max {r['null_max']:7.3f} | p = {r['p']:.4g}")
os.makedirs(os.path.join(HERE, "results", a.tag), exist_ok=True)
json.dump(dict(config=vars(a), **res), open(os.path.join(HERE, "results", a.tag, "permutation_null.json"), "w"), indent=2)
