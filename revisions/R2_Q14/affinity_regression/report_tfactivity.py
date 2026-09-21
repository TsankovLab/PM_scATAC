"""Fit one affinity-regression configuration on TCGA MESO and write everything needed to
compare with Hmeljak et al. 2018 tab 2B / Fig 2D-G:
  results/<tag>/tf_activity.csv        TF x tumour inferred activity (W P')
  results/<tag>/bap1_tf_ttest.csv      our 2B-equivalent, merged with the paper's values
  results/<tag>/summary.json           configuration, model size, agreement statistics
Run with the CV-selected configuration for the final result."""
import argparse, json, os, csv, numpy as np
from affreg import ar_train, ar_model2w
from mesobap1 import HERE, load_variant, paper_2B, ttest, compare

ap = argparse.ArgumentParser()
ap.add_argument("--release", default="meso_tcga_pan_can_atlas_2018"); ap.add_argument("--universe", default="top5000")
ap.add_argument("--train", default="all63"); ap.add_argument("--dnorm", default="binary"); ap.add_argument("--ynorm", default="colnorm")
ap.add_argument("--lam", type=float, default=0.1); ap.add_argument("--specA", type=float, default=0.95)
ap.add_argument("--specB", type=float, default=0.9); ap.add_argument("--tag", default="provisional")
a = ap.parse_args()

V = load_variant(a.release, a.universe)
Y, D, P, st, samples = V["Y"], V["D"].astype(float), V["P"], V["bap1_status"], V["samples"]
if a.train == "labelled54":
    k = st != "NA"; Y, P, st, samples = Y[:, k], P[k], st[k], samples[k]
    Y = Y - Y.mean(1, keepdims=True); P = P - P.mean(0, keepdims=True)
if a.dnorm == "colnorm": D = D / np.sqrt((D ** 2).sum(0, keepdims=True))
if a.ynorm == "colnorm": Y = Y / np.sqrt((Y ** 2).sum(0, keepdims=True))
m = ar_train(D, P, Y, a.lam, 0.0, a.specA, a.specB)
act = ar_model2w(m) @ P.T
tfs = [str(t) for t in V["tfs"]]
res = ttest(act, tfs, st); paper = paper_2B(); comp = compare(res, paper)

out = os.path.join(HERE, "results", a.tag); os.makedirs(out, exist_ok=True)
with open(os.path.join(out, "tf_activity.csv"), "w", newline="") as h:
    w = csv.writer(h); w.writerow(["TF"] + [str(s) for s in samples])
    for i, t in enumerate(tfs): w.writerow([t] + [f"{v:.6g}" for v in act[i]])
with open(os.path.join(out, "samples.csv"), "w", newline="") as h:
    w = csv.writer(h); w.writerow(["sample", "bap1_status"]); w.writerows(zip(samples, st))
with open(os.path.join(out, "bap1_tf_ttest.csv"), "w", newline="") as h:
    w = csv.writer(h); w.writerow(["TF", "estimate", "t", "p", "padj", "paper_estimate", "paper_p", "paper_padj"])
    for i, t in enumerate(tfs):
        pv = paper.get(t, (np.nan, np.nan, np.nan))
        w.writerow([t, res["est"][i], res["t"][i], res["p"][i], res["padj"][i], pv[0], pv[1], pv[2]])
json.dump(dict(config=vars(a), da=m["da"], db=m["db"], n_iter=m["n_iter"], n_groups=res["n"],
               agreement={k: (float(v) if isinstance(v, (float, np.floating)) else v) for k, v in comp.items()}),
          open(os.path.join(out, "summary.json"), "w"), indent=2, default=float)
print(json.dumps(dict(tag=a.tag, n=res["n"], **{k: comp[k] for k in ("rho_est", "rho_signed_logp", "pct_same_sign",
      "paper_sig", "our_sig", "both_sig_same_sign", "IRF8_rank", "EGR2_rank", "YY1_rank")}), default=float))
