"""Quick lambda landscape under per-sample unit-norm Y (diagnostic, not selection)."""
import itertools, numpy as np
from affreg import ar_train, ar_model2w
from mesobap1 import load_variant, paper_2B, ttest, compare
paper = paper_2B()
V = load_variant("meso_tcga_pan_can_atlas_2018", "top5000")
for train, Dn in itertools.product(["all63", "labelled54"], ["binary", "colnorm"]):
    Y, D, P, st = V["Y"], V["D"].astype(float), V["P"], V["bap1_status"]
    if train == "labelled54":
        k = st != "NA"; Y, P, st = Y[:, k], P[k], st[k]
        Y = Y - Y.mean(1, keepdims=True); P = P - P.mean(0, keepdims=True)
    if Dn == "colnorm": D = D / np.sqrt((D ** 2).sum(0, keepdims=True))
    Y = Y / np.sqrt((Y ** 2).sum(0, keepdims=True))
    for sA, sB in [(1.0, 1.0), (0.95, 0.9), (0.9, 0.9)]:
        for lam in [0.001, 0.003, 0.01, 0.03, 0.1, 0.2, 0.3, 0.5]:
            m = ar_train(D, P, Y, lam, 0.0, sA, sB)
            c = compare(ttest(ar_model2w(m) @ P.T, V["tfs"], st), paper)
            nz = int((np.abs(m["beta"]) > 0).sum())
            print(f"{train:10s} D={Dn:7s} spec={sA}/{sB} lam={lam:<5} nz_beta={nz:4d}/{m['beta'].size:4d} | "
                  f"rho_est={c['rho_est']:+.3f} rho_slogp={c['rho_signed_logp']:+.3f} same={c['pct_same_sign']:.0f}% "
                  f"sig both/ours/paper={c['both_sig_same_sign']}/{c['our_sig']}/{c['paper_sig']} | "
                  + " ".join(f"{t}#{c[t+'_rank']}{'+' if c[t+'_est']>0 else '-'}" for t in ("IRF8", "EGR2", "YY1")), flush=True)
