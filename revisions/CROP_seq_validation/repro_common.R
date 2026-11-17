# Shared helpers for the CROP-seq validation figures (bulk RNA-seq KO cohort vs CROP-seq screen).
# Every script does: source(file.path(<script dir>, "repro_common.R"))
source(local({a <- grep("^--file=", commandArgs(FALSE), value = TRUE); for (d in c(if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))), getwd())) { while (!file.exists(file.path(d, "config.R")) && d != dirname(d)) d <- dirname(d); if (file.exists(file.path(d, "config.R"))) break }; Sys.getenv("MESO_CONFIG", file.path(d, "config.R"))}))   # data / results / reference roots: see config.R at the repo root
suppressPackageStartupMessages(library(matrixStats))
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))   # figures use unicode glyphs (delta, rho, minus)

PROC <- Sys.getenv("PROC", paste0(DATA_DIR, "/CROPseq_validation_bulkRNA/04_processed_counts"))   # star_counts/ (STAR GeneCounts), sample_sheet/
COUNTS_DIR <- Sys.getenv("COUNTS_DIR", file.path(PROC, "star_counts"))
OUT  <- Sys.getenv("CROP_OUT", unset = file.path(RESULTS_DIR, "CROP_seq_validation"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

LINES <- c("MSTO-211H", "NCI-H2452", "NCI-H28")
LINE_COL <- c("CROP-seq" = "#4a3aa7", "MSTO-211H" = "#1baf7a", "NCI-H2452" = "#e87ba4", "NCI-H28" = "#008300")
KO_COL <- c(NTC = "#7a7a7a", TEAD4 = "#eb6834", TWIST1 = "#1baf7a", PITX1 = "#eda100", SOX9 = "#e87ba4", TEAD2 = "#4a3aa7", BPTF = "#008300")   # one colour per KO
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"
GREY <- "#BEBEBE"; BROWN <- "#A52A2A"                       # colours of the original CRISPR_cropseq.R plots
DIV <- colorRampPalette(c("#2a78d6", "#f0efec", "#eb6834"))(101)

## ---------------------------------------------------------------- cohort
# The bulk RNA-seq KO cohort (43 samples, as in ../KO35_bulkRNA): 35 KO samples (6 KOs: TEAD4, TWIST1, PITX1, SOX9, TEAD2, BPTF x MSTO-211H, NCI-H2452, NCI-H28)
# + 8 non-targeting controls.
# library = sample name (FASTQ / STAR files); line = cell line; KO = knocked-out gene (NTC = non-targeting control).
COHORT <- local({
  ko <- c("MSTO_TEAD4_R1", "MSTO_TEAD4_R2", "MSTO_TWIST1_R2", "MSTO_TWIST1_R3", "MSTO_PITX1_R1", "MSTO_PITX1_R3", "MSTO_SOX9_R2", "MSTO_SOX9_R3", "MSTO_TEAD2_R2", "MSTO_TEAD2_R3",
          "MSTO_BPTF_R1", "MSTO_BPTF_R3",
          "H2452_TEAD4_R1", "H2452_TEAD4_R3", "H2452_TWIST1_R1", "H2452_PITX1_R1", "H2452_PITX1_R3", "H2452_SOX9_R2", "H2452_SOX9_R3", "H2452_TEAD2_R1", "H2452_TEAD2_R2",
          "H2452_BPTF_R1", "H2452_BPTF_R2",
          "H28_TEAD4_R1", "H28_TEAD4_R3", "H28_TWIST1_R1", "H28_TWIST1_R3", "H28_PITX1_R1", "H28_PITX1_R3", "H28_SOX9_R1", "H28_SOX9_R2", "H28_TEAD2_R1", "H28_TEAD2_R2",
          "H28_BPTF_R1", "H28_BPTF_R3")
  ntc <- c("MSTO_NTC_R1", "MSTO_NTC_R2", "MSTO_NTC_R3", "H2452_NTC_R1", "H2452_NTC_R2", "H2452_NTC_R3", "H28_NTC_R1", "H28_NTC_R2")
  PFX <- c(MSTO = "MSTO-211H", H2452 = "NCI-H2452", H28 = "NCI-H28"); lib <- c(ko, ntc)
  d <- data.frame(library = lib, line = unname(PFX[sub("_.*", "", lib)]), KO = sub("^[^_]+_([^_]+)_.*", "\\1", lib), stringsAsFactors = FALSE)
  rownames(d) <- d$library; stopifnot(nrow(d) == 43); d })

# STAR gene counts (reverse-stranded column) of the given libraries, collapsed to gene symbols (genes x samples)
load_counts <- function(samples = COHORT$library) {
  counts <- read.delim(file.path(COUNTS_DIR, "gene_counts_matrix_stranded_reverse.tsv"), row.names = 1, check.names = FALSE)
  sym <- read.delim(file.path(COUNTS_DIR, "gene_id_to_symbol.tsv"), header = FALSE, col.names = c("gene_id", "symbol"))
  sym <- sym[sym$gene_id %in% rownames(counts) & !is.na(sym$symbol) & sym$symbol != "", ]
  rowsum(as.matrix(counts[sym$gene_id, samples]), group = sym$symbol)
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
