import time, numpy as np
from mesobap1 import *
paper = paper_2B()
print(f"paper tab 2B: {len(paper)} TFs, padj<0.01: {sum(v[2] < 0.01 for v in paper.values())}")
PARAMS = [(0.01, 0.0, 1.0, 1.0), (0.1, 0.0, 0.95, 0.9), (0.001, 0.0, 1.0, 1.0)]
best = None
for rel, univ in VARIANTS:
    V = load_variant(rel, univ)
    for lam, rs, sA, sB in PARAMS:
        t0 = time.time()
        m, W, act = fit_full(V, lam, rs, sA, sB)
        res = ttest(act, V["tfs"], V["bap1_status"])
        c = compare(res, paper)
        print(f"{rel:28s} {univ:8s} lam={lam:<6} spec={sA}/{sB} da={m['da']:2d} db={m['db']:2d} it={m['n_iter']:5d} {time.time()-t0:5.1f}s | "
              f"n={res['n']} rho_est={c['rho_est']:+.3f} rho_slogp={c['rho_signed_logp']:+.3f} same_sign={c['pct_same_sign']:.0f}% "
              f"sig(ours/paper/both)={c['our_sig']}/{c['paper_sig']}/{c['both_sig_same_sign']} | "
              + " ".join(f"{t}:{c[t+'_est']:+.3g}(#{c[t+'_rank']})" for t in FOCUS))
        if best is None or c["rho_signed_logp"] > best[0]:
            best = (c["rho_signed_logp"], rel, univ, lam, rs, sA, sB, res)
_, rel, univ, lam, rs, sA, sB, res = best
print(f"\nbest quick setting by agreement (for inspection only, not selection): {rel} {univ} lam={lam} spec={sA}/{sB}")
o = np.argsort(res["p"])
print(f"{'TF':8s} {'our_est':>9s} {'our_p':>9s} {'our_padj':>9s} | {'paper_est':>9s} {'paper_padj':>10s}")
for i in o[:20]:
    t = res["TF"][i]; pv = paper.get(t, (np.nan, np.nan, np.nan))
    print(f"{t:8s} {res['est'][i]:+9.3g} {res['p'][i]:9.2g} {res['padj'][i]:9.2g} | {pv[0]:+9.3g} {pv[2]:10.2g}")
print("\nsame best setting, grouping possibly_inactivated WITH inactivated:")
V = load_variant(rel, univ); _, _, act = fit_full(V, lam, rs, sA, sB)
r2 = ttest(act, V["tfs"], V["bap1_status"], grp=("inactivated", "possibly_inactivated"))
c2 = compare(r2, paper)
print(f"  n={r2['n']} rho_est={c2['rho_est']:+.3f} rho_slogp={c2['rho_signed_logp']:+.3f} sig both={c2['both_sig_same_sign']}/{c2['paper_sig']}")
