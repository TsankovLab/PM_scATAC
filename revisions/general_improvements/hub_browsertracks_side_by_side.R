###############################################################################
# general_improvements -- browser tracks for several cHub regions displayed side by side.
#
# Regions come from the cHub call in
#   main/scatac_ArchR/hubs_obj_cor_0.3_md_12500_dgs_0_min_peaks_5/
#     hub_regions.bed      chr, start, end, HUB id, strand, overlapping gene symbols
#     global_hubs_obj.rds  $hubsMerged[[i]] = the peaks constituting hub i (data.frame)
#                          $hubs_id         = HUB ids, aligned to hubsMerged
#
# Tracks: ArchR plotBrowserTrack on the main project, groupBy celltype_lv1, cell types in the
# order supplied by the user. Each hub is one panel; panels are then placed in a single row
# (and, for the long list, also wrapped over several rows so the figure stays readable).
# The featureTrack under each panel shows that hub's own constituent peaks, so the reader sees
# which peaks define the hub rather than just the interval.
#
# NB coverage is normalised per group (ReadsInTSS) but the y-scale is chosen PER PANEL by
# ArchR, as is standard for browser tracks; heights are therefore comparable within a locus,
# not between loci.
#
# NB HUB6929 was requested but does not exist in this hub object (ids run HUB1..HUB6719);
# it is reported and skipped rather than silently dropped.
#
# Output: Plots/hub_browsertracks_set1_row.pdf        12 hubs, single row
#         Plots/hub_browsertracks_set2_row.pdf        full list, single row
#         Plots/hub_browsertracks_set2_wrapped.pdf    full list, 6 per row
#         hub_regions_plotted.csv
###############################################################################
suppressPackageStartupMessages({ library(ArchR); library(GenomicRanges); library(grid)
                                 library(gridExtra); library(ggplot2); library(ggrepel)
                                 library(RColorBrewer); library(paletteer); library(circlize) })
## ggplot2 + ggrepel must be attached even when only re-drawing cached panels: the stored
## gtables contain ggplot/ggrepel grobs whose draw methods live in those namespaces, and
## without them every panel silently renders blank.
addArchRThreads(4); addArchRGenome("hg38")
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUT  <- file.path(ROOT, "git_repo_claude", "general_improvements")
setwd(OUT); dir.create("Plots", showWarnings = FALSE)
HUBDIR <- file.path(ROOT, "main", "scatac_ArchR", "hubs_obj_cor_0.3_md_12500_dgs_0_min_peaks_5")

celltype_order <- c('Malignant','Mesothelium','Alveolar','Fibroblasts','SmoothMuscle','Endothelial',
                    'Myeloid','T_cells','NK','B_cells','Plasma','pDCs')
## project cell-type palette; palettes.R renames SmoothMuscle -> "Smooth Muscle", so map back
source(file.path(ROOT, "git_repo", "utils", "palettes.R"))
CT_PAL <- palette_celltype_lv1
names(CT_PAL)[names(CT_PAL) == "Smooth Muscle"] <- "SmoothMuscle"
stopifnot(all(celltype_order %in% names(CT_PAL)))
CT_PAL <- CT_PAL[celltype_order]
## SET1 is the requested set (the one to use); SET2 is the full list seen earlier, kept so the
## other hubs stay available from the same cache (deduplicated, first occurrence wins).
SET1 <- c('HUB2543','HUB5226','HUB1063','HUB15','HUB2387','HUB25','HUB85','HUB9','HUB32',
          'HUB34','HUB4599','HUB6622')
SET2 <- unique(c('HUB522','HUB5713','HUB277','HUB32','HUB68','HUB1459','HUB26','HUB16','HUB48','HUB77',
                 'HUB6929','HUB2371','HUB2543','HUB5226','HUB1063','HUB15','HUB25','HUB2387','HUB85',
                 'HUB9','HUB32','HUB34','HUB4599','HUB6622'))
PAD <- 0.05      # fraction of hub width added either side so edge peaks are not clipped
PER_ROW <- 6

## ---- gene track: strand colours and centred labels ---------------------------------------
## ArchR hard-codes the gene track's strand palette and label placement inside the internal
## .geneTracks(); plotBrowserTrack() exposes neither. The function is therefore re-written with
## its own source: forward strand brown, reverse strand black, and the gene symbol centred
## above its gene body (ArchR anchors labels at the gene start/end and pushes them sideways).
GENE_PLUS  <- "brown"      # forward strand
GENE_MINUS <- "black"      # reverse strand
patch_gene_track <- function(colorPlus = GENE_PLUS, colorMinus = GENE_MINUS,
                             labelSize = 3, nudge_y = 0.18) {
  ## deparse() wraps at width.cutoff, so the source is rejoined with newlines and every
  ## pattern tolerates a line break inside the expression (\\s+ between tokens).
  txt <- paste(deparse(ArchR:::.geneTracks), collapse = "\n")
  subs <- list(
    c("x\\s*=\\s*start,\\s*y\\s*=\\s*cluster,\\s*label\\s*=\\s*symbol",
      "x = (start + end)/2, y = cluster, label = symbol"),
    c("x\\s*=\\s*end,\\s*y\\s*=\\s*cluster,\\s*label\\s*=\\s*symbol",
      "x = (start + end)/2, y = cluster, label = symbol"),
    c("nudge_x\\s*=\\s*-0\\.01\\s*\\*\\s*\\(end\\(region\\)\\s*-\\s*start\\(region\\)\\)", "nudge_x = 0"),
    c("nudge_x\\s*=\\s*\\+0\\.01\\s*\\*\\s*\\(end\\(region\\)\\s*-\\s*start\\(region\\)\\)", "nudge_x = 0"),
    c("nudge_y\\s*=\\s*-0\\.25", paste0("nudge_y = ", nudge_y)),
    c("nudge_y\\s*=\\s*0\\.25",  paste0("nudge_y = ", nudge_y)),
    c("ggrepel::geom_label_repel", "ggrepel::geom_text_repel"))
  for (u in subs) {
    if (!grepl(u[1], txt)) stop("gene-track patch: pattern not found -> ", u[1])
    txt <- gsub(u[1], u[2], txt)
  }
  f <- eval(parse(text = txt))
  formals(f)$colorPlus <- colorPlus; formals(f)$colorMinus <- colorMinus
  formals(f)$labelSize <- labelSize
  environment(f) <- asNamespace("ArchR")     # keep ArchR internals in scope
  assignInNamespace(".geneTracks", f, ns = "ArchR")
  invisible(TRUE)
}
## ArchR labels the coverage facets with a strip on the right of every panel. Side by side that
## is 12 repeated boxes per hub, so the name is drawn inside each facet at the top left instead
## and the strip column is dropped altogether (see prep()).
## theme_ArchR sets panel.border = element_rect(...) for every track, and plotBrowserTrack's
## borderWidth argument does not switch it off, so the theme is wrapped to blank it. This drops
## the frame around each coverage row, the peak track and the gene track.
THEME_ORIG <- ArchR::theme_ArchR
assignInNamespace("theme_ArchR", function(...) {
  THEME_ORIG(...) + ggplot2::theme(panel.border = ggplot2::element_blank())
}, ns = "ArchR")

## plotBrowserTrack() never passes `pal` down to .loopTracks(), so its formal default is what
## gets used: set it to the project purple ramp (palette_enrichment, pale -> dark purple, so
## stronger co-accessibility draws darker).
LOOP_FN <- ArchR:::.loopTracks
formals(LOOP_FN)$pal <- palette_enrichment
assignInNamespace(".loopTracks", LOOP_FN, ns = "ArchR")

BULK_ORIG <- ArchR:::.bulkTracks        # pristine copy: used for the panels that carry no labels
patch_bulk_track <- function(labelSize = 3.4, labelled = TRUE) {
  txt <- paste(deparse(BULK_ORIG), collapse = "\n")   # always patch the pristine source
  ## black outline around each coverage profile: `colour` is mapped to group in the global
  ## aes, so an explicit colour on the layer overrides it while the group fill is kept
  ga <- 'geom_area\\(stat\\s*=\\s*"identity"\\)'
  if (!grepl(ga, txt)) stop("bulk-track patch: geom_area pattern not found")
  txt <- sub(ga, 'geom_area(stat = "identity", colour = "black", linewidth = 0.05)', txt)
  pat <- 'facet_wrap\\(facets\\s*=\\s*~group,\\s*strip\\.position\\s*=\\s*"right",\\s*ncol\\s*=\\s*1\\)'
  if (!grepl(pat, txt)) stop("bulk-track patch: facet_wrap pattern not found")
  if (labelled) txt <- sub(pat, paste0('facet_wrap(facets = ~group, strip.position = "right", ncol = 1) + ',
    'geom_text(data = data.frame(group = unique(df$group)), ',
    'mapping = aes(x = -Inf, y = Inf, label = group), hjust = -0.07, vjust = 1.45, ',
    'size = ', labelSize, ', inherit.aes = FALSE, colour = "black")'), txt)
  f <- eval(parse(text = txt))
  environment(f) <- asNamespace("ArchR")
  assignInNamespace(".bulkTracks", f, ns = "ArchR")
  invisible(TRUE)
}
## the plain variant differs only by the absence of the in-facet names: it still needs the
## black outline, so it is the same patch with the label layer switched off
unpatch_bulk_track <- function() patch_bulk_track(labelled = FALSE)
patch_gene_track()
message("gene track: + strand ", GENE_PLUS, ", - strand ", GENE_MINUS, ", labels centred above gene bodies")

## ---- hub regions ------------------------------------------------------------------------
bed <- read.table(file.path(HUBDIR, "hub_regions.bed"), sep = "\t", stringsAsFactors = FALSE)
colnames(bed) <- c("chr", "start", "end", "hub", "strand", "genes")
hobj <- readRDS(file.path(HUBDIR, "global_hubs_obj.rds"))

resolve <- function(ids) {
  miss <- setdiff(ids, bed$hub)
  if (length(miss)) message("NOT FOUND in this hub object, skipped: ", paste(miss, collapse = ", "))
  ids[ids %in% bed$hub]
}
mk_gr <- function(ids) {
  b <- bed[match(ids, bed$hub), ]
  w <- b$end - b$start
  gr <- GRanges(b$chr, IRanges(pmax(1, round(b$start - PAD * w)), round(b$end + PAD * w)))
  mcols(gr)$hub <- b$hub; mcols(gr)$genes <- b$genes; mcols(gr)$kb <- round(w / 1e3, 1)
  gr
}
## each hub's own peaks, as one GRanges for the feature track
hub_peaks <- function(ids) {
  i <- match(ids, hobj$hubs_id)
  do.call(c, lapply(i, function(k) {
    d <- hobj$hubsMerged[[k]]
    GRanges(as.character(d$seqnames), IRanges(d$start, d$end))
  }))
}

## pairwise co-accessibility links between hub peaks: already stored anchor-to-anchor with a
## `value` column, i.e. exactly the format ArchR's loopTrack expects
HUB_LINKS <- hobj$peakLinks[["main"]]

ids2 <- resolve(SET2); ids1 <- resolve(SET1)
gr2 <- mk_gr(ids2)
info <- data.frame(hub = mcols(gr2)$hub, chr = as.character(seqnames(gr2)),
                   plot_start = start(gr2), plot_end = end(gr2),
                   hub_kb = mcols(gr2)$kb, genes = mcols(gr2)$genes,
                   n_hub_peaks = sapply(match(ids2, hobj$hubs_id), function(k) nrow(hobj$hubsMerged[[k]])),
                   in_set1 = mcols(gr2)$hub %in% ids1)
write.csv(info, "hub_regions_plotted.csv", row.names = FALSE)
cat("=== hub regions plotted ===\n"); print(info, row.names = FALSE)

## ---- tracks -----------------------------------------------------------------------------
proj <- loadArchRProject(file.path(ROOT, "main", "scatac_ArchR"), showLogo = FALSE)
proj <- proj[getCellColData(proj)$celltype_lv1 %in% celltype_order, ]

tracks <- function(gr, ids, labelled) {
  if (labelled) patch_bulk_track() else unpatch_bulk_track()
  feats <- GRangesList(`cHub peaks` = hub_peaks(ids))
  lnk <- subsetByOverlaps(HUB_LINKS, gr, type = "within")   # keep the figure light
  loops <- GRangesList(`cHub links` = lnk)
  message("hub links in view: ", length(lnk))
  p <- plotBrowserTrack(proj, region = gr, groupBy = "celltype_lv1", useGroups = celltype_order,
                        features = feats, loops = loops,
                        plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
                        sizes = c(10, 0.45, 1.0, 1.9), tileSize = 200, baseSize = 11,
                        facetbaseSize = 11, pal = CT_PAL,
                        borderWidth = 0, tickWidth = 0,   # no panel/axis frames around the tracks
                        title = "")
  if (!is.list(p)) p <- list(p)
  p
}
## ---- layout -----------------------------------------------------------------------------
## Each ArchR panel is a 92x14 gtable. Side by side, columns 5 (y-axis title), 6 (feature-row
## label) and 8 (cell-type strip) would repeat 22 times and leave almost no width for data, so
## they are blanked on all but one panel and their column widths zeroed.
## NB ArchR states the coverage range in the y-axis TITLE rather than as tick labels, and the
## range is chosen per region. Blanking that title would silently drop it, so the range is
## parsed out and written into each panel's own header instead.
COL <- c(ylab = 5, axis_l = 6, strip = 8, title = 7)

labels_of <- function(x) {                                  # recursive grob label hunt
  out <- character(0)
  if (!is.null(x$label)) out <- c(out, as.character(x$label))
  for (ch in c(x$children, x$grobs)) out <- c(out, labels_of(ch))
  out
}
yrange <- function(g) {
  txt <- unlist(lapply(g$grobs[which(g$layout$l == COL[["ylab"]])], labels_of))
  m <- regmatches(txt, regexpr("Range \\([^)]*\\)", txt))
  if (length(m)) sub("Range \\((.*)\\)", "\\1", m[1]) else NA_character_
}
blank_cols <- function(g, cols) {
  if (!length(cols)) return(g)
  for (i in which(g$layout$l %in% cols)) g$grobs[[i]] <- nullGrob()
  g$widths[cols] <- unit(0, "cm")
  g
}
retitle <- function(g, txt) {
  i <- which(g$layout$name == "title")
  if (length(i)) {
    g$grobs[[i[1]]] <- textGrob(txt, gp = gpar(fontsize = 9.2, lineheight = 1.2))
    g$heights[g$layout$t[i[1]]] <- unit(2.2, "lines")
  }
  g
}
## one panel, trimmed: y-axis title never kept (range moves to the header); the feature-row
## label ("cHub peaks") only on the first panel; the cell-type strip only on the last.
## the genomic axis is dropped from every track; each panel's coordinates are printed in its
## header, so 22 repeated coordinate axes add nothing but height
drop_x_axis <- function(g) {
  i <- grep("^(axis-b|xlab-b)", g$layout$name)
  if (length(i)) { for (k in i) g$grobs[[k]] <- nullGrob(); g$heights[g$layout$t[i]] <- unit(0, "cm") }
  g
}
## the loop track carries its own colour key; side by side that is one legend per hub, so it is
## kept on the right-most panel only
drop_right_guide <- function(g) {
  i <- grep("^guide-box-right$", g$layout$name)
  if (length(i)) { for (k in i) g$grobs[[k]] <- nullGrob()
                   g$widths[unique(g$layout$l[i])] <- unit(0, "cm") }
  g
}
drop_bottom_guide <- function(g) {            # ArchR leaves a stray "strand" key under the page
  i <- grep("^guide-box-bottom$", g$layout$name)
  if (length(i)) { for (k in i) g$grobs[[k]] <- nullGrob(); g$heights[g$layout$t[i]] <- unit(0, "cm") }
  g
}
## rows occupied only by "background" are pure whitespace between the stacked tracks; below the
## coverage block they are the gaps around the peak and gene tracks, so they are collapsed.
compact_rows <- function(g, keep = unit(0.02, "cm")) {
  last_bulk <- max(g$layout$b[grepl("^panel-1-", g$layout$name)])
  for (r in seq_len(nrow(g))) {
    if (r <= last_bulk) next
    occ <- unique(g$layout$name[g$layout$t <= r & g$layout$b >= r])
    if (!length(setdiff(occ, "background"))) g$heights[r] <- keep
  }
  g
}
prep <- function(g, i, n, hdr) {
  g <- compact_rows(drop_x_axis(drop_bottom_guide(retitle(g, hdr))))
  if (i < n) g <- drop_right_guide(g)          # legend only on the right-most panel
  ## strip dropped on ALL panels now that cell types are named inside the facets
  g <- blank_cols(g, c(COL[["ylab"]], COL[["strip"]], if (i > 1) COL[["axis_l"]]))
  ## with the panel frames removed, white space is the only thing separating neighbouring
  ## hubs, so the outer gutter is widened rather than left at a hairline
  g$widths[c(1, ncol(g))] <- unit(0.13, "cm")
  g
}
## ids are stored as HUBnnn in the hub object but displayed as cHubnnn. The header carries the
## id and the genomic region only; per-panel coverage ranges stay in hub_regions_plotted.csv.
header <- function(gr, i, yr) sprintf("%s\n%s:%s-%s",
  sub("^HUB", "cHub", mcols(gr)$hub[i]), as.character(seqnames(gr))[i],
  format(start(gr)[i], big.mark = ","), format(end(gr)[i], big.mark = ","))

## assemble a row of panels, aligning row heights across panels (gene tracks differ in height)
row_of <- function(panels, gr, idx) {
  yr <- sapply(PP$plain, yrange)
  n <- length(idx)
  tr <- lapply(seq_along(idx), function(k) {
    ## cell-type names appear only on the left-most panel of the row
    g <- if (k == 1) PP$labelled[[idx[k]]] else PP$plain[[idx[k]]]
    prep(g, k, n, header(gr, idx[k], yr))                   # yr stays globally indexed, as idx is
  })
  do.call(cbind, c(tr, list(size = "max")))
}

## Two renders of every region: one carrying the in-facet cell-type names (used for whichever
## panel sits left-most in a given layout) and one without (every other panel).
CACHE <- "hub_browsertrack_panels_v11.rds"  # v11 = larger text again
if (file.exists(CACHE)) {
  message("reusing cached panels (", CACHE, "); delete it to re-render")
  PP <- readRDS(CACHE)
} else {
  message("plotting ", length(gr2), " regions x ", length(celltype_order), " cell types (x2 variants) ...")
  PP <- list(labelled = tracks(gr2, ids2, TRUE), plain = tracks(gr2, ids2, FALSE))
  saveRDS(PP, CACHE)
}
P2 <- PP$plain
YR <- sapply(PP$plain, yrange)
info$y_range <- YR
write.csv(info, "hub_regions_plotted.csv", row.names = FALSE)
cat("\nper-panel coverage range (ArchR auto-scales each region):\n")
print(data.frame(hub = info$hub, y_range = info$y_range), row.names = FALSE)

PW <- 2.1            # inches of drawing width per panel
LAB <- 2.6           # inches for the shared y label + cell-type strip

CAPTION <- sprintf(paste("Genes: %s = forward strand, %s = reverse strand.",
  "cHub peaks = constituent peaks of that hub; arcs = pairwise co-accessibility between them (value = correlation).",
  "Coverage is normalised by ReadsInTSS and scaled per panel; per-panel ranges are listed in hub_regions_plotted.csv."),
  GENE_PLUS, GENE_MINUS)
add_caption <- function() grid.text(CAPTION, x = unit(0.5, "npc"), y = unit(3, "mm"),
  gp = gpar(fontsize = 8.5, col = "grey25"))

draw_row <- function(file, idx) {
  g <- row_of(P2, gr2, idx)
  pdf(file, width = PW * length(idx) + LAB, height = 8.1)
  grid.newpage(); grid.draw(g); add_caption(); dev.off()
  cat(sprintf("%-46s %2d panels  %5.1f x 8.1 in  %d KB\n", basename(file), length(idx),
              PW * length(idx) + LAB, round(file.size(file) / 1024)))
}

## full list, one row
draw_row("Plots/hub_browsertracks_set2_row.pdf", seq_along(gr2))
## full list, wrapped: each row is its own aligned block, stacked on one page
idx_chunks <- split(seq_along(gr2), ceiling(seq_along(gr2) / PER_ROW))
pdf("Plots/hub_browsertracks_set2_wrapped.pdf", width = PW * PER_ROW + LAB, height = 7.8 * length(idx_chunks) + 0.3)
grid.newpage()
pushViewport(viewport(layout = grid.layout(length(idx_chunks), 1)))
for (k in seq_along(idx_chunks)) {
  pushViewport(viewport(layout.pos.row = k, layout.pos.col = 1))
  grid.draw(row_of(P2, gr2, idx_chunks[[k]])); popViewport()
}
popViewport(); add_caption(); dev.off()
cat(sprintf("%-46s %2d rows of %d\n", "hub_browsertracks_set2_wrapped.pdf", length(idx_chunks), PER_ROW))
## the explicit 12-hub set, in its own order
draw_row("Plots/hub_browsertracks_set1_row.pdf", match(ids1, mcols(gr2)$hub))

cat("DONE\n")
