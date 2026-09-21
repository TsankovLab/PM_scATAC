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
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"
GREY <- "#BEBEBE"; BROWN <- "#A52A2A"                       # colours of the original CRISPR_cropseq.R plots
DIV <- colorRampPalette(c("#2a78d6", "#f0efec", "#eb6834"))(101)

## ---------------------------------------------------------------- data
# bulk RUN1 counts collapsed to gene symbol -> CP10K (genes x samples) + sample metadata
load_bulk <- function() {
  counts <- read.delim(file.path(PROC, "counts/gene_counts_matrix_stranded_reverse.tsv"), row.names = 1, check.names = FALSE)
  sym <- read.delim(file.path(PROC, "counts/gene_id_to_symbol.tsv"), header = FALSE, col.names = c("gene_id", "symbol"))
  sym <- sym[sym$gene_id %in% rownames(counts) & !is.na(sym$symbol) & sym$symbol != "", ]
  cs <- rowsum(as.matrix(counts[sym$gene_id, ]), group = sym$symbol)
  meta <- read.delim(file.path(PROC, "sample_sheet/sample_metadata_RUN1.tsv"))
  meta <- meta[meta$status == "OK" & meta$Cell_Line %in% LINES & meta$fastq_sample_id %in% colnames(cs), ]
  rownames(meta) <- meta$fastq_sample_id
  cp10k <- sweep(cs, 2, colSums(cs), "/") * 1e4
  list(cp10k = cp10k[, rownames(meta), drop = FALSE], meta = meta)
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
