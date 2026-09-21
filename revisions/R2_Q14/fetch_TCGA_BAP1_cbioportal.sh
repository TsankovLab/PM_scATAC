#!/bin/bash
###############################################################################
# Fetch TCGA MESO (PanCancer Atlas) BAP1 genetic data from the public cBioPortal
# API into tcga_genetic/, and flatten it to TSV.  Study meso_tcga_pan_can_atlas_2018,
# 87 samples -- the same 87 tumours as bulkRNA_meso/bulk_RNA_studies.rds$tcga.
#   mutations       MAF-level somatic calls        (sample list: sequenced, 86)
#   gistic          discrete copy number -2..2     (sample list: cna, 87)
#   log2CNA         continuous copy number         (sample list: log2CNA, 87)
#   structural var. fusions / rearrangements       (profiled samples: sv, 86)
#   rppa            BAP1 PROTEIN, reverse-phase array (sample list: rppa) --
#                   the closest TCGA readout to clinical IHC
# BAP1 Entrez ID 8314.  Only public, de-identified TCGA data is queried.
###############################################################################
set -euo pipefail
cd "$(dirname "$0")"; O=tcga_genetic; mkdir -p $O
API=https://www.cbioportal.org/api; S=meso_tcga_pan_can_atlas_2018; G=8314
post(){ curl -sf -X POST -H "Content-Type: application/json" -d "$2" "$API/$1"; }
post "molecular-profiles/${S}_mutations/mutations/fetch?projection=DETAILED" \
     "{\"entrezGeneIds\":[$G],\"sampleListId\":\"${S}_sequenced\"}" > $O/BAP1_mutations.json
post "molecular-profiles/${S}_gistic/discrete-copy-number/fetch?discreteCopyNumberEventType=ALL&projection=SUMMARY" \
     "{\"entrezGeneIds\":[$G],\"sampleListId\":\"${S}_cna\"}" > $O/BAP1_gistic.json
post "molecular-profiles/${S}_log2CNA/molecular-data/fetch?projection=SUMMARY" \
     "{\"entrezGeneIds\":[$G],\"sampleListId\":\"${S}_log2CNA\"}" > $O/BAP1_log2CNA.json
post "structural-variant/fetch" \
     "{\"entrezGeneIds\":[$G],\"molecularProfileIds\":[\"${S}_structural_variants\"]}" > $O/BAP1_structural_variants.json
post "molecular-profiles/${S}_rppa/molecular-data/fetch?projection=SUMMARY" \
     "{\"entrezGeneIds\":[$G],\"sampleListId\":\"${S}_rppa\"}" > $O/BAP1_rppa.json
for L in sequenced cna log2CNA sv rppa; do
  curl -sf "$API/sample-lists/${S}_${L}/sample-ids" > $O/samplelist_${L}.json
done
python3 - <<PY
import json, csv
O="$O"
def w(name, rows, cols):
    with open(f"{O}/{name}", "w", newline="") as f:
        x = csv.writer(f, delimiter="\t"); x.writerow(cols)
        for r in rows: x.writerow([r.get(c, "") for c in cols])
m = json.load(open(f"{O}/BAP1_mutations.json"))
for r in m:
    a, b = r.get("tumorAltCount", -1), r.get("tumorRefCount", -1)
    r["vaf"] = round(a/(a+b), 3) if a not in (None,-1) and b not in (None,-1) and a+b > 0 else ""
w("BAP1_mutations.tsv", m, ["sampleId","mutationType","proteinChange","variantType","vaf","startPosition","referenceAllele","variantAllele"])
w("BAP1_gistic.tsv",    json.load(open(f"{O}/BAP1_gistic.json")),  ["sampleId","alteration"])
w("BAP1_log2CNA.tsv",   json.load(open(f"{O}/BAP1_log2CNA.json")), ["sampleId","value"])
w("BAP1_rppa.tsv",      json.load(open(f"{O}/BAP1_rppa.json")),    ["sampleId","value"])
sv = json.load(open(f"{O}/BAP1_structural_variants.json"))
w("BAP1_structural_variants.tsv", sv, ["sampleId","site1HugoSymbol","site2HugoSymbol","eventInfo","variantClass"])
for L in ["sequenced","cna","log2CNA","sv","rppa"]:
    with open(f"{O}/samplelist_{L}.txt","w") as f: f.write("\n".join(json.load(open(f"{O}/samplelist_{L}.json")))+"\n")
print("flattened to TSV")
PY
