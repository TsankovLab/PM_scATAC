# general_improvements — side-by-side cHub browser tracks

`hub_browsertracks_side_by_side.R` draws ArchR browser tracks for several cHub regions as
panels in a single row, so loci can be compared across cell types at a glance.

## Inputs

- **Regions** — `main/scatac_ArchR/hubs_obj_cor_0.3_md_12500_dgs_0_min_peaks_5/`
  - `hub_regions.bed` — chr, start, end, HUB id, strand, overlapping gene symbols (6,719 hubs)
  - `global_hubs_obj.rds` — `$hubsMerged[[i]]` gives the peaks constituting hub *i*
- **Signal** — `main/scatac_ArchR`, `groupBy = celltype_lv1`, 49,849 cells, existing group
  coverages (no recomputation needed)
- **Cell-type order** (top to bottom, as supplied): Malignant, Mesothelium, Alveolar,
  Fibroblasts, SmoothMuscle, Endothelial, Myeloid, T_cells, NK, B_cells, Plasma, pDCs
- **Colours** — `palette_celltype_lv1` from `git_repo/utils/palettes.R`, passed to
  `plotBrowserTrack(pal = ...)`. That script renames `SmoothMuscle` to `Smooth Muscle`, so the
  name is mapped back before use and every cell type is asserted to have a colour.

## Outputs

| file | contents |
|---|---|
| `Plots/hub_browsertracks_set1_row.pdf` | the explicit 12-hub list, one row (11 panels, 21.3 x 7.8 in) |
| `Plots/hub_browsertracks_set2_row.pdf` | the full list, one row (22 panels, 40 x 7.8 in) |
| `Plots/hub_browsertracks_set2_wrapped.pdf` | same 22, wrapped 6 per row |
| `hub_regions_plotted.csv` | region, hub width, genes, n hub peaks, per-panel y range |
| `hub_browsertrack_panels_v2.rds` | cached rendered panels — delete to force a re-render |

Each panel carries a bulk coverage track per cell type, a `cHub peaks` feature track showing
that hub's own constituent peaks, and a gene track. Regions are padded by 5% either side so
edge peaks are not clipped.

## Gene track

ArchR hard-codes the gene track's strand palette (`red`/`dodgerblue2`) and label placement
inside its internal `.geneTracks()`, and `plotBrowserTrack()` exposes neither. The script
rewrites that function from its own source, so the track now shows:

- **forward strand in brown, reverse strand in black**
- **the gene symbol centred above its own gene body** (ArchR anchors each label at the gene
  start/end and pushes it sideways, which detaches the name from the feature)
- `geom_text_repel` rather than boxed labels, and no stray "strand" key under the page
- labels sitting close to the bodies (`nudge_y = 0.18`) and a compact track (`sizes` 1.9)

**The genomic axis is dropped from every track.** Repeated coordinate axes only cost height,
so each panel's coordinates are printed in its header instead. The header carries the cHub id
and the region only; hub ids are stored as `HUBnnn` in the hub object and displayed as
`cHubnnn`, and the per-panel coverage range is recorded in `hub_regions_plotted.csv` rather
than on the figure.

## Links track

`global_hubs_obj.rds$peakLinks$main` stores each pairwise link anchor-to-anchor with a `value`
column, which is exactly ArchR's loop format, so it feeds a `loopTrack` directly. Arcs are drawn
between co-accessible peaks of the hub and coloured by correlation with `palette_enrichment`
(pale to dark purple). `plotBrowserTrack()` never passes `pal` down to `.loopTracks()`, so that
palette is installed as the function's formal default. The colour key is kept on the right-most
panel only. Note these links are the global co-accessibility set used to build the hubs, not a
per-cell-type calculation: the arcs are identical whichever cell type carries the signal.

## Coverage profiles

Profiles carry a thin black outline. `colour` is mapped to group in ArchR's global aes, so an
explicit `colour` on the `geom_area` layer overrides it while keeping the group fill; the width
is 0.05, since anything heavier swamps the fill at this panel width. Both the labelled and
unlabelled bulk-track variants are derived from a pristine copy of `.bulkTracks`, because
re-patching an already-patched function fails to match its own source patterns.

## Coverage track

ArchR labels each coverage facet with a strip on the right of every panel. Side by side that is
12 repeated boxes per hub, so `.bulkTracks()` is patched the same way as the gene track to draw
the cell-type name **inside each facet at the top left**, and the strip column is dropped from
every panel. To label only the left-most panel instead, pass `i > 1` to the strip blanking in
`prep()` and restrict the `geom_text` layer to the first region.

Whitespace rows below the coverage block (the gaps around the peak and gene tracks) are
collapsed to 0.04 cm by `compact_rows()`, which is what closes the space between the peak track
and the gene annotation.

The patch is applied with `assignInNamespace`, so it affects only the running session and
leaves the installed ArchR untouched. It asserts on each source pattern it expects and stops
with a clear error if a future ArchR version changes them.

## Two points to note when using the figure

**HUB6929 does not exist in this hub object.** Ids run HUB1–HUB6719, so it is reported and
skipped; the 12-hub set therefore has 11 panels. If it came from the `cor_0.2` object or from
the permutation-FDR hub set, say so and it can be pulled from there.

**The coverage scale is chosen per region**, as ArchR does by default, so bar heights are
comparable between cell types within a panel but not between panels. ArchR states that range
in the y-axis title, which would otherwise be dropped when the repeated axes are removed, so
it is parsed out and printed in each panel header (`y 0-0.17`) and recorded in the CSV.
Ranges across these hubs run 0.15–0.33. If a single common scale is wanted, re-render with an
explicit `ylim` in `plotBrowserTrack`.

## Layout note

Each ArchR panel is a 92 x 14 gtable in which the y-axis title, the feature-row label and the
cell-type strip each occupy a fixed ~2 cm column. Repeated 22 times these consume the whole
page, so they are blanked on all but one panel (strip on the right-most, feature label on the
left-most) and their column widths zeroed. Panels are then combined with `cbind(size = "max")`
so that row heights stay aligned despite gene tracks of differing height.

## Redrawing from the cache

The cached panels store ggplot/ggrepel grobs whose draw methods live in those namespaces, so
any session that re-draws them must attach `ggplot2` and `ggrepel` as well as `grid` — without
them every panel renders blank with no error.

## Colour clash to be aware of

The cell-type palette gives Endothelial `brown` and Alveolar `black`, which are also the
forward- and reverse-strand gene colours. They never share a track, so the figure stays
readable, but if it matters for the legend either the strand colours or those two cell types
should be changed.
