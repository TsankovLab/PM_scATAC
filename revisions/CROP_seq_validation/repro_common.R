# Shared helpers for the CROP-seq validation figures (bulk RNA-seq RUN1 vs CROP-seq screen).
# Every script does: source(file.path(<script dir>, "repro_common.R"))
suppressPackageStartupMessages(library(matrixStats))
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))   # figures use unicode glyphs (delta, rho, minus)

PROC <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CROPseq_validation_bulkRNA/processed"
PP   <- file.path(PROC, "de", "perturbation_panel")
OUT  <- Sys.getenv("CROP_OUT", unset = file.path(PP, "figs_repro"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

LINES <- c("NCI-H2052", "MSTO-211H", "NCI-H2452", "NCI-H28")
THREE <- c("MSTO-211H", "NCI-H2052", "NCI-H2452")
LINE_COL <- c("CROP-seq" = "#4a3aa7", "NCI-H2052" = "#eda100", "MSTO-211H" = "#1baf7a",
              "NCI-H2452" = "#e87ba4", "NCI-H28" = "#008300")
KO_COL <- c(NTC = "#7a7a7a", TCF3 = "#2a78d6", TEAD4 = "#eb6834", TWIST1 = "#1baf7a", PITX1 = "#eda100", SOX9 = "#e87ba4", TEAD2 = "#4a3aa7", BPTF = "#008300", HMGA1 = "#b3400f",
            MEF2A = "#8fa88a", MEF2D = "#5a7d9a", RUNX1 = "#8c564b", RUNX2 = "#17becf", SNAI2 = "#bcbd22")       # one colour per KO, shared by the guide PCA and the cross-line heatmaps
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"
GREY <- "#BEBEBE"; BROWN <- "#A52A2A"                       # colours of the original CRISPR_cropseq.R plots
DIV <- colorRampPalette(c("#2a78d6", "#f0efec", "#eb6834"))(101)

# H28_NTC_R1/R2/R3: real, fully-processed samples (STAR ReadsPerGene.out.tab present in counts/) that are absent from sample_metadata_RUN1.tsv entirely.
# Index sequences recovered from sample_sheet/SampleSheet.csv and mapped to i5/i7_index_num via the existing rows' seq -> index_num correspondence.
# GENOTYPE REASSIGNMENT (registered spec25_NTC_reassignment_table.csv): PCA of the full dataset (fig34/fig35 in spec25_bulk_PCA_prefilter.R) showed
# these 3 samples cluster tightly with NCI-H2052 (not NCI-H28, despite the "H28_" name) -- specifically with H2052_NTC_R3, which sits with them at
# PC1 ~= -40 to -44 while every other NCI-H28 sample (incl. the KO samples) sits at PC1 ~= +65. Cell_Line is therefore corrected to NCI-H2052 here;
# Cell_Line_original stays NCI-H28 (the as-sequenced/nominal label) -- mirroring how the sheet itself already handles H2052_NTC_R1/R2 (sequenced as
# NCI-H2052, genotyped and corrected to NCI-H28: the reciprocal half of what looks like a wholesale NTC-batch mislabelling between these 2 lines).
H28_NTC_NATIVE <- data.frame(
  Sample_ID = c("H28_NTC_R1", "H28_NTC_R2", "H28_NTC_R3"), Cell_Line = "NCI-H2052", Gene_Control = "NTC",
  i5_index_num = c(18L, 17L, 17L), i7_index_num = c(9L, 5L, 23L), i5_seq = c("GCGCCTGT", "ATATTCAC", "ATATTCAC"), i7_seq = c("GATCAG", "ACAGTG", "GAGTGG"),
  status = "OK", fastq_sample_id = c("H28_NTC_R1", "H28_NTC_R2", "H28_NTC_R3"), orig_sample_ids_in_fastq = c("H28_NTC_R1", "H28_NTC_R2", "H28_NTC_R3"),
  Cell_Line_original = "NCI-H28", label_corrected = "True", stringsAsFactors = FALSE)

# Samples visually flagged as outliers in the per-line, cell-cycle-regressed PCA of the full audited universe (spec25_bulk_PCA_prefilter.R,
# fig37/38: load_full_audited_bulk() + top 5000 variable genes per line, cc regressed). For NCI-H2052: every sample with PC2 > 20 (the 3
# reassigned H28_NTC_R1/2/3 plus H2052_TEAD4_R3). For the other 3 lines: specific samples named directly by inspection of that same PCA.
# Used to highlight (not remove) these samples in spec25's fig39 and spec23_QC_RUN1.R's fig27_QC_RUN1_screenKOs_byKO.
VISUAL_OUTLIERS <- c("H2052_TEAD4_R3", "H28_NTC_R1", "H28_NTC_R2", "H28_NTC_R3",   # NCI-H2052 panel, PC2 > 20
                     "MSTO_TWIST1_R2", "MSTO_SOX9_R3",                             # MSTO-211H
                     "H2452_BPTF_R1", "H2452_TWIST1_R2",                           # NCI-H2452
                     "H28_TCF3_R1", "H28_HMGA1_R2")                                # NCI-H28

## ---------------------------------------------------------------- data
# bulk RUN1 counts collapsed to gene symbol -> CP10K (genes x samples) + sample metadata
load_bulk <- function() {
  counts <- read.delim(file.path(PROC, "counts/gene_counts_matrix_stranded_reverse.tsv"), row.names = 1, check.names = FALSE)
  sym <- read.delim(file.path(PROC, "counts/gene_id_to_symbol.tsv"), header = FALSE, col.names = c("gene_id", "symbol"))
  sym <- sym[sym$gene_id %in% rownames(counts) & !is.na(sym$symbol) & sym$symbol != "", ]
  cs <- rowsum(as.matrix(counts[sym$gene_id, ]), group = sym$symbol)
  meta <- read.delim(file.path(PROC, "sample_sheet/sample_metadata_RUN1.tsv"))
  # NOTE: H28_NTC_R1/2/3 (see H28_NTC_NATIVE below) are deliberately NOT added here -- several precomputed per-sample inputs used elsewhere in this
  # pipeline (e.g. cm_sample_scores_AddModuleScore_bulk_RUN1.csv) predate them and don't have rows for these samples, which breaks scripts that join on
  # sample ID (fig3 et al.). H28_NTC_NATIVE is added only in the one script that needs it and is self-contained (spec25_bulk_PCA_prefilter.R).
  meta <- meta[meta$status == "OK" & meta$Cell_Line %in% LINES & meta$fastq_sample_id %in% colnames(cs), ]
  rownames(meta) <- meta$fastq_sample_id
  cp10k <- sweep(cs, 2, colSums(cs), "/") * 1e4
  list(cp10k = cp10k[, rownames(meta), drop = FALSE], meta = meta, counts = cs[, rownames(meta), drop = FALSE])
}

# the FULL, AUDITED bulk sample universe: sample list built directly from counts/ (every *_ReadsPerGene.out.tab file), NOT just the curated sheet --
# catches samples entirely missing from the sheet (H2452_TWIST1_R2, H28_HMGA1_R2, H28_NTC_R1/2/3; see H28_NTC_NATIVE above). Every KO except
# SNAI2/RUNX1/RUNX2 (not in the CROP-seq screen), the 4 real lines, ALL statuses kept (meta$removed_run1 flags status != "OK"). This is the single
# source of truth for "how many samples/replicates exist per line" -- used by spec25_bulk_PCA_prefilter.R and spec23_QC_RUN1.R alike so they never
# silently disagree on the sample set. meta$Cell_Line is the corrected/genotyped line (includes the registered NTC reassignment via H28_NTC_NATIVE
# and the sheet's own label_corrected rows); meta$Cell_Line_original is the as-named/as-sequenced label; meta$flag records how each row was
# determined (clean / swapped (genotype-corrected) / mixed (pending deconv. w/ Lily) / undocumented (missing from sheet)).
load_full_audited_bulk <- function() {
  counts <- read.delim(file.path(PROC, "counts/gene_counts_matrix_stranded_reverse.tsv"), row.names = 1, check.names = FALSE)
  sym <- read.delim(file.path(PROC, "counts/gene_id_to_symbol.tsv"), header = FALSE, col.names = c("gene_id", "symbol"))
  sym <- sym[sym$gene_id %in% rownames(counts) & !is.na(sym$symbol) & sym$symbol != "", ]
  cs <- rowsum(as.matrix(counts[sym$gene_id, ]), group = sym$symbol)

  all_ids <- sub("_ReadsPerGene\\.out\\.tab$", "", list.files(file.path(PROC, "counts"), pattern = "_ReadsPerGene\\.out\\.tab$"))
  sheet_full <- read.delim(file.path(PROC, "sample_sheet/sample_metadata_RUN1.tsv"))
  undoc <- setdiff(all_ids, sheet_full$fastq_sample_id)
  PFX <- c(H2052 = "NCI-H2052", MSTO = "MSTO-211H", H2452 = "NCI-H2452", H28 = "NCI-H28")
  parse_orphan <- function(id) {
    m <- regmatches(id, regexec("^(H2052|MSTO|H2452|H28)_([A-Za-z0-9]+)_R([0-9]+)$", id))[[1]]
    if (length(m) != 4) return(NULL)
    data.frame(Sample_ID = id, Cell_Line = PFX[[m[2]]], Cell_Line_original = PFX[[m[2]]], Gene_Control = m[3], status = "OK", fastq_sample_id = id, stringsAsFactors = FALSE)
  }
  orphan_rows <- do.call(rbind, lapply(undoc, parse_orphan))

  meta <- sheet_full
  meta$flag <- ifelse(meta$status == "MIXED_PENDING_DECONV", "mixed (pending deconv. w/ Lily)", ifelse(meta$Cell_Line != meta$Cell_Line_original, "swapped (genotype-corrected)", "clean"))
  meta <- meta[, c("Sample_ID", "Cell_Line", "Cell_Line_original", "Gene_Control", "status", "fastq_sample_id", "flag")]
  orphan_rows$flag <- "undocumented (missing from sheet)"
  meta <- rbind(meta, orphan_rows[, colnames(meta)])
  i <- match(H28_NTC_NATIVE$Sample_ID, meta$Sample_ID); meta$Cell_Line[i] <- H28_NTC_NATIVE$Cell_Line   # registered NTC reassignment
  meta <- meta[meta$Cell_Line_original %in% LINES & meta$fastq_sample_id %in% colnames(cs) & !meta$Gene_Control %in% c("SNAI2", "RUNX1", "RUNX2"), ]
  rownames(meta) <- meta$fastq_sample_id
  meta$rep <- sub(".*_(R[0-9]+)$", "\\1", meta$Sample_ID)
  meta$swapped <- meta$Cell_Line != meta$Cell_Line_original
  meta$removed_run1 <- meta$status != "OK"
  cp10k <- sweep(cs, 2, colSums(cs), "/") * 1e4; cp10k <- cp10k[, rownames(meta), drop = FALSE]
  list(cs = cs[, rownames(meta), drop = FALSE], cp10k = cp10k, lg = log2(cp10k + 1), meta = meta)
}

# log2(mean CP10K + 1) per KO group inside one cell line (genes x KO groups)
group_log2avg <- function(cp10k, meta, cl) {
  m <- meta[meta$Cell_Line == cl, ]
  gs <- sort(unique(m$Gene_Control), method = "radix")
  sapply(setNames(gs, gs), function(g) log2(rowMeans(cp10k[, rownames(m)[m$Gene_Control == g], drop = FALSE]) + 1))
}

reindex_rows <- function(m, rows) { out <- m[match(rows, rownames(m)), , drop = FALSE]; rownames(out) <- rows; out }   # missing rows -> NA

zscore_by_gene <- function(mat) {           # R's scale(t(x)): z-score each gene across KO groups, NA -> 0
  z <- (mat - rowMeans(mat)) / rowSds(mat)
  z[!is.finite(z)] <- 0
  z
}

# hierarchical order (complete linkage, 1 - Pearson) of the ROWS of x
corr_order <- function(x) {
  if (nrow(x) < 3) return(seq_len(nrow(x)))
  hclust(as.dist(1 - cor(t(x))), "complete")$order
}

rms_scale <- function(m) { r <- sqrt(rowMeans(m^2)); r[r < 1e-12] <- NA; r }   # constant features (float noise after centring) -> dropped to 0

# per-dataset feature normalisation: optional centring across KOs, then equal-RMS scaling
normalise_profiles <- function(D, centre) {
  lapply(D, function(m) {
    x <- if (centre) m - rowMeans(m) else m
    x <- x / rms_scale(x); x[!is.finite(x)] <- 0; x
  })
}

## ---------------------------------------------------------------- plotting
save_fig <- function(name, w, h, fun, res = 200) {
  png(file.path(OUT, paste0(name, ".png")), width = w, height = h, units = "in", res = res, type = "cairo"); fun(); dev.off()
  cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h); fun(); dev.off()
  message("saved ", name)
}

heat_z <- function(M) t(M)[, nrow(M):1, drop = FALSE]        # matrix (rows top->bottom) -> image() layout

draw_colorbar <- function(x, zlim, label, cex = 0.6) {        # horizontal colour bar in the bottom margin; x = c(left, right) in figure fraction
  hin <- par("din")[2]
  op <- par(fig = c(x, 0, 0.75 / hin), new = TRUE, oma = c(0, 0, 0, 0), mar = c(2.0, 0, 1.0, 0), cex = 1); on.exit(par(op))
  image(seq(zlim[1], zlim[2], length.out = 101), 1, matrix(seq(zlim[1], zlim[2], length.out = 101), ncol = 1),
        col = DIV, zlim = zlim, axes = FALSE, xlab = "", ylab = "")
  axis(1, cex.axis = cex, lwd = 0.5, mgp = c(1, 0.3, 0)); mtext(label, side = 1, line = 1.2, cex = cex, col = INK2); box(lwd = 0.5)
}

violin_at <- function(v, at, col, half = 0.45) {              # width-scaled violin (trimmed at the data range), like the original
  d <- density(v, from = min(v), to = max(v), n = 256)
  w <- d$y / max(d$y) * half
  polygon(c(at - w, rev(at + w)), c(d$x, rev(d$x)), col = adjustcolor(col, 0.6), border = NA)
}
