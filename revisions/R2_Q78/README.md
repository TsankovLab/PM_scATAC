# R2_Q78 — browser tracks of cHub2711, cHub3722 and cHub3728

Reviewer 2, Q78 asks for track plots of three cHubs: 2711 (near HAND2) and 3722 / 3728 (near MEF2C).

- `R2_Q78_cHub_tracks.R` — loads the main ArchR project (`main/scatac_ArchR`) and the hub set of the main analysis
  (`hubs_obj_cor_0.3_md_12500_dgs_0_min_peaks_5`, 6,719 hubs), checks where each hub lies relative to the gene it is
  expected to be near, and draws the regions side by side with `plotTracksSideBySide()` (`git_repo/utils/scATAC_functions.R`).
- `run_step.sh` — runs a script of this folder in the `meso_scatac` environment.

## Hub positions (hg38)

| hub | region | annotated gene | distance to expected gene |
|---|---|---|---|
| HUB2711 | chr4:173,411,073-173,419,134 | SCRG1 | HAND2, 106 kb |
| HUB3722 | chr5:74,250,424-74,258,705 | LINC01335 | MEF2C, 14.5 Mb (not near MEF2C) |
| HUB3727 | chr5:89,267,692-89,279,891 | MEF2C | MEF2C, 363 kb |
| HUB3728 | chr5:89,290,877-89,303,441 | MEF2C | MEF2C, 387 kb |

HUB3722 is not near MEF2C in this hub set; the second MEF2C hub next to HUB3728 is HUB3727. Both versions are plotted.

## Outputs (`Plots/`)

- `R2_Q78_cHubs_2711_3722_3728.pdf` — the three ids as requested
- `R2_Q78_cHubs_2711_3727_3728.pdf` — with the two MEF2C hubs
- `R2_Q78_cHubs_gene_context.pdf` — wide windows showing each hub together with HAND2 / MEF2C
- `R2_Q78_hub_positions.csv` — coordinates, number of peaks, nearest gene, distance to HAND2 / MEF2C
