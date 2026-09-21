"""Pre-specified selection of the final affinity-regression configuration.

Rule (fixed before the complete grids were inspected):
  PRIMARY  the configuration with the highest mean 10-fold CV reconstruction correlation
           over both grids -- the criterion Osmanbeyoglu et al. used.  Where several grid
           points tie to 1e-4 (lambda is inert once the spectrum truncation dominates), the
           first in grid order is taken; they give the same TF activities.
  LITERAL  the best-CV configuration restricted to the most literal reading of the
           published methods: binary D, no Y/D normalisation, top-5000 genes, all tumours
           with RPPA -- reported alongside so the effect of the unstated choices is visible.
  SPREAD   agreement with tab 2B across the 20 best-CV configurations, to show how much the
           conclusion depends on the exact pick.
Agreement with tab 2B is never used to choose.
Writes prepared/final_selection.json."""
import csv, json, statistics, sys

def fl(x):
    try: return float(x)
    except (TypeError, ValueError): return float("nan")

rows = []
for f, g in [("prepared/cv_grid_results.csv", 1), ("prepared/cv_grid2_results.csv", 2)]:
    try:
        for r in csv.DictReader(open(f)):
            if r.get("error") or not r.get("cv_pearson"): continue
            r.setdefault("train", "all63"); r.setdefault("Dnorm", "binary"); r.setdefault("Ynorm", "none")
            if g == 1 and fl(r.get("rsL2", 0)) != 0: continue          # rsL2>0 not carried to reporting
            r["grid"] = g; rows.append(r)
    except FileNotFoundError:
        pass
partial = "--partial" in sys.argv
print(f"{len(rows)} completed grid points{' (PARTIAL grids)' if partial else ''}")
key = lambda r: -fl(r["cv_pearson"])
rows.sort(key=key)
def cfg(r):
    return dict(release=r["release"], universe=r["universe"], train=r["train"], dnorm=r["Dnorm"], ynorm=r["Ynorm"],
                lam=fl(r["lam"]), specA=fl(r["specA"]), specB=fl(r["specB"]), cv_pearson=fl(r["cv_pearson"]),
                rho_est=fl(r["rho_est"]), rho_signed_logp=fl(r["rho_signed_logp"]),
                both_sig_same_sign=int(fl(r["both_sig_same_sign"])), IRF8_rank=int(fl(r["IRF8_rank"])),
                EGR2_rank=int(fl(r["EGR2_rank"])), YY1_rank=int(fl(r["YY1_rank"])))
primary = cfg(rows[0])
lit = [r for r in rows if r["Dnorm"] == "binary" and r["Ynorm"] == "none" and r["universe"] == "top5000" and r["train"] == "all63"]
literal = cfg(lit[0]) if lit else None
top = [cfg(r) for r in rows[:20]]
spread = {k: (min(t[k] for t in top), statistics.median(t[k] for t in top), max(t[k] for t in top))
          for k in ("cv_pearson", "rho_est", "rho_signed_logp", "both_sig_same_sign", "IRF8_rank", "EGR2_rank", "YY1_rank")}
out = dict(primary=primary, literal=literal, top20_spread_min_median_max=spread, n_grid_points=len(rows), partial=partial)
json.dump(out, open("prepared/final_selection.json" if not partial else "prepared/final_selection_PARTIAL.json", "w"), indent=2)
for name, c in [("PRIMARY", primary), ("LITERAL", literal)]:
    if c: print(f"{name:8s}: {c['release']} {c['universe']} {c['train']} D={c['dnorm']} Y={c['ynorm']} lam={c['lam']} "
                f"spec={c['specA']}/{c['specB']} | cv_r={c['cv_pearson']:.4f} | rho_est={c['rho_est']:+.3f} "
                f"sig both={c['both_sig_same_sign']}/28 IRF8#{c['IRF8_rank']} EGR2#{c['EGR2_rank']} YY1#{c['YY1_rank']}")
print("top-20 CV spread (min / median / max):")
for k, v in spread.items(): print(f"   {k:20s} {v[0]:.4g} / {v[1]:.4g} / {v[2]:.4g}")
