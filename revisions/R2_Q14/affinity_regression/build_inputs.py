"""
Build the affinity-regression inputs for the TCGA MESO BAP1 reproduction.

D (genes x TFs): MSigDB v6.1 c3.tft (TRANSFAC) target sets -- the exact release that
  reproduces the paper's published YY1 (427/427) and IRF8 (248/248) target lists.
  One set per TF; for TFs with several motifs the LARGEST set is used (the rule that
  reproduces both published lists).  The 141 TFs are the paper's (tab 2B).  121 map by
  HGNC symbol/alias from the TRANSFAC factor name; the rest are assigned by hand from
  TRANSFAC nomenclature (MANUAL below) -- ambiguous ones are logged in
  prepared/tf_motif_crosswalk.csv with a 'confidence' column.
Y: log10(RSEM+1), mean-centred per gene, two gene universes:
  top5000   the 5,000 most variable genes (Osmanbeyoglu et al. 2017 methods)
  motifhit  every expressed gene hit by >=1 of the 141 TF sets
P: RPPA, mean-centred per antibody (from prepare step).
"""
import csv, re, collections, numpy as np, openpyxl

H = list(csv.DictReader(open("hgnc_complete_set.txt"), delimiter="\t"))
sym = {r["symbol"].upper(): r["symbol"] for r in H}
alias = collections.defaultdict(set)
for r in H:
    for k in ("alias_symbol", "prev_symbol"):
        for a in (r.get(k) or "").split("|"):
            a = a.strip().strip('"').upper()
            if a: alias[a].add(r["symbol"]); alias[a.replace("-", "")].add(r["symbol"])
sets = {}
for line in open("c3.tft.v6.1.symbols.gmt"):
    p = line.rstrip("\n").split("\t"); sets[p[0]] = set(p[2:])
IUPAC = set("ACGTRYKMSWBDHVN")
def factor(n):
    t = n.split("_")
    if len(t) > 1 and len(t[0]) >= 5 and set(t[0]) <= IUPAC: t = t[1:]
    while len(t) > 1 and re.fullmatch(r"(\d+|Q\d+|B\d*|C|DR\d+|UNKNOWN)", t[-1]): t = t[:-1]
    return "_".join(t)
def cands(f):
    F = f.upper(); out = set()
    for v in {F, F.replace("GAMMA","G").replace("ALPHA","A").replace("BETA","B").replace("DELTA","D"), F.replace("_","")}:
        if v in sym: out.add(sym[v])
        out |= alias.get(v, set())
    return out
auto = collections.defaultdict(list)
for n in sets:
    for s in cands(factor(n)): auto[s].append(n)

wb = openpyxl.load_workbook("../21598290cd180804-sup-205173_2_supp_5073240_pg14sl.xlsx", read_only=True, data_only=True)
TFS = [str(r[0]) for r in list(wb["2B_bap1_inferred_TF_activity"].iter_rows(values_only=True))[1:] if r and r[0]]
yy1_pub  = {str(r[1]) for r in list(wb["2C_YY1 target genes"].iter_rows(values_only=True))[1:] if r and r[1]}
irf8_pub = {str(r[1]) for r in list(wb["2D_IRF8 target genes"].iter_rows(values_only=True))[1:] if r and r[1]}

rx = lambda p: [n for n in sets if re.search(p, n)]
MANUAL = {   # TF: (candidate sets, confidence, note)
 "RELA":   (["NFKAPPAB65_01"], "high", "p65"),
 "TAL1":   (rx(r"TAL1"), "high", ""),
 "IKZF2":  (["IK2_01"], "medium", "TRANSFAC IK2 (Ikaros isoform); name-based"),
 "CREB1":  ([n for n in rx(r"(^|_)CREB_") if "CREBP1" not in n and "TAX" not in n], "high", ""),
 "USF1":   ([n for n in rx(r"USF") if "USF2" not in n], "high", ""),
 "CREBBP": (["P300_01"], "low", "only P300 motif; shared with EP300"),
 "NFYA":   (rx(r"NFY"), "medium", "NFY motifs; split between NFYA/NFYB"),
 "NFYB":   (rx(r"NFY"), "low", "NFY motifs; split between NFYA/NFYB"),
 "IRF9":   (["ISRE_01"], "medium", "ISGF3 element"),
 "LMO2":   (rx(r"LMO2"), "high", ""),
 "NFATC2": (rx(r"NFAT"), "medium", "NFAT motifs; split between NFATC2/NFATC4"),
 "NFATC4": (rx(r"NFAT"), "low", "NFAT motifs; split between NFATC2/NFATC4"),
 "NKX3.1": (["NKX3A_01"], "high", ""),
 "NR1H3":  (rx(r"LXR"), "high", ""),
 "SPI1":   (rx(r"PU1"), "high", ""),
 "RUNX2":  ([n for n in rx(r"AML") if "AML1" not in n], "low", "generic AML motif; AML1 sets go to RUNX1"),
 "SMAD1":  (["SMAD_Q6"], "medium", "generic SMAD motif"),
 "MEF2D":  ([n for n in rx(r"MEF2") if "RSRFC4" not in n], "low", "MEF2 motifs other than RSRFC4 (->MEF2A)"),
 "RARA":   (rx(r"(^|_)RAR|RARA|DR5|RXR"), "low", "RAR/RXR motif if present"),
 "ADD1":   (rx(r"SREBP"), "medium", "ADD1 is an SREBF1 alias; second SREBP set"),
}
used, cross = set(), []
order = [t for t in TFS if t not in MANUAL] + [t for t in TFS if t in MANUAL]   # auto first
for t in order:
    if t in MANUAL:
        cand, conf, note = MANUAL[t]
    else:
        cand, conf, note = auto.get(t, []), "high", "HGNC symbol/alias"
    cand = sorted(set(cand), key=lambda n: -len(sets[n]))
    free = [n for n in cand if n not in used]
    pick = (free or cand or [None])[0]
    if pick: used.add(pick)
    cross.append(dict(TF=t, motif_set=pick or "", n_targets=len(sets[pick]) if pick else 0,
                      n_candidates=len(cand), confidence=conf if pick else "unmapped",
                      shared_set=bool(pick and not free), note=note))
with open("prepared/tf_motif_crosswalk.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(cross[0])); w.writeheader(); w.writerows(cross)
C = {c["TF"]: c for c in cross}
print(f"TFs: {len(TFS)} | mapped {sum(bool(c['motif_set']) for c in cross)} | unmapped {[c['TF'] for c in cross if not c['motif_set']]}")
print("confidence:", collections.Counter(c["confidence"] for c in cross))
print("sets shared by >1 TF:", [(c['TF'], c['motif_set']) for c in cross if c['shared_set']])
print("low/medium manual:", [(c['TF'], c['motif_set']) for c in cross if c['confidence'] in ('low','medium')])
print(f"check published lists: YY1 -> {C['YY1']['motif_set']} ({len(sets[C['YY1']['motif_set']] & yy1_pub)}/427) | "
      f"IRF8 -> {C['IRF8']['motif_set']} ({len(sets[C['IRF8']['motif_set']] & irf8_pub)}/248)")

tfs = [c["TF"] for c in cross if c["motif_set"]]
tfs = [t for t in TFS if t in tfs]                     # paper order
for rel in ["meso_tcga", "meso_tcga_pan_can_atlas_2018"]:
    z = np.load(f"prepared/{rel}.npz", allow_pickle=True)
    genes = list(z["genes"]); gi = {g: i for i, g in enumerate(genes)}
    Yl, fe = z["Y_log10"], z["frac_expressed"]
    Dfull = np.zeros((len(genes), len(tfs)))
    for j, t in enumerate(tfs):
        for g in sets[C[t]["motif_set"]]:
            if g in gi: Dfull[gi[g], j] = 1
    tf_expr = {t: (fe[gi[t]] if t in gi else np.nan) for t in tfs}
    n_fail = sum(1 for v in tf_expr.values() if not (v >= 0.4))
    sd = Yl.std(axis=1)
    univ = {"top5000": np.argsort(-sd)[:5000],
            "motifhit": np.where((Dfull.sum(axis=1) > 0) & (fe >= 0.4))[0]}
    for u, rows in univ.items():
        Y = Yl[rows]; Y = Y - Y.mean(axis=1, keepdims=True)
        D = Dfull[rows]
        keep_tf = D.sum(axis=0) > 0
        np.savez_compressed(f"prepared/input_{rel}_{u}.npz", Y=Y, D=D, genes=np.array(genes)[rows],
                            tfs=np.array(tfs), P=z["P"], antibodies=z["antibodies"],
                            samples=z["samples"], bap1_status=z["bap1_status"])
        print(f"{rel:30s} {u:9s}: Y {Y.shape[0]} genes x {Y.shape[1]} | D {D.shape} | "
              f"TFs with >=1 target {int(keep_tf.sum())} | targets/TF median {int(np.median(D.sum(0)))} | "
              f"YY1 targets in universe {int(D[:, tfs.index('YY1')].sum())} | IRF8 {int(D[:, tfs.index('IRF8')].sum())}")
    print(f"   TFs below the 40%-expressed filter in {rel}: {n_fail} "
          f"({[t for t, v in tf_expr.items() if not (v >= 0.4)][:12]})")
