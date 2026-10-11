# Spec 53: screen-vs-bulk agreement for specific KO/module pairs: TEAD4 in Cm20, TCF3 in Cm10, SOX9 in Cm1 (env PAIRS="TEAD4:Cm20,TCF3:Cm10,SOX9:Cm1",
# env EXCL_LINES default NCI-H28 removed). Same Cm score (AddModuleScore, top-50 genes) and NTC - KO direction as spec45/51/52.
# Per pair: NTC-KO diff in the screen and in each retained bulk line (H2052 bulk shown as same-line reference, 1 NTC), sign agreement, and how the KO ranks
# among all KOs for that module (rank 1 = largest NTC-KO diff) in the screen and in each line -- i.e. is the effect specific to that KO.
# needs: git_repo/files/scRNA_PM_modules.rds, spec_common.R, repro_common.R
# outputs: spec53_module_KO_pairs_agreement.csv, spec53_module_KO_pairs_rank.csv, fig3K_module_KO_pairs_agreement.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat))
PAIRS <- do.call(rbind, lapply(strsplit(strsplit(Sys.getenv("PAIRS", "TEAD4:Cm20,TCF3:Cm10,SOX9:Cm1"), ",")[[1]], ":"), function(x) data.frame(KO = x[1], MOD = x[2])))
EXCL <- setdiff(strsplit(Sys.getenv("EXCL_LINES", "NCI-H28"), ",")[[1]], ""); USE <- setdiff(LINES, c(EXCL, "NCI-H2052"))   # H2052 bulk (1 NTC) is shown separately as reference
cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds")
top50 <- lapply(cm_genes, head, 50); CM <- names(top50); options(Seurat.object.assay.version = "v3")
srt <- suppressWarnings(AddModuleScore(load_screen(), features = top50, name = "Cm_", seed = 1))
sc <- as.matrix(srt@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(sc) <- CM; sgrp <- as.character(srt$merged_call)
b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1)
sb <- CreateSeuratObject(counts = as.matrix(b$counts)); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(lg)); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1))
bs <- as.matrix(sb@meta.data[rownames(meta), paste0("Cm_", seq_along(CM))]); colnames(bs) <- CM
KOS <- intersect(ORDER, unique(sgrp)); SRC <- c("CROP-seq (H2052)", "NCI-H2052", USE)
getx <- function(s, mod) if (s == "CROP-seq (H2052)") data.frame(g = sgrp, x = sc[, mod]) else { i <- meta$Cell_Line == s; data.frame(g = meta$Gene_Control[i], x = bs[i, mod]) }
diffs <- function(s, mod) { d <- getx(s, mod); n <- d$x[d$g == "NTC"]; sapply(KOS, function(k) { v <- d$x[d$g == k]; if (!length(v) || !length(n)) NA else mean(n) - mean(v) }) }

rows <- list(); rk <- list()
for (i in seq_len(nrow(PAIRS))) { k <- PAIRS$KO[i]; m <- PAIRS$MOD[i]
  for (s in SRC) { d <- getx(s, m); n <- d$x[d$g == "NTC"]; v <- d$x[d$g == k]
    p <- if (s == "CROP-seq (H2052)") wilcox.test(n, v)$p.value else if (length(n) > 1 && length(v) > 1) tryCatch(t.test(n, v)$p.value, error = function(e) NA) else NA
    df <- diffs(s, m); rows[[length(rows) + 1]] <- data.frame(KO = k, module = m, source = s, n_NTC = length(n), n_KO = length(v), mean_NTC = mean(n), mean_KO = mean(v), NTC_minus_KO = mean(n) - mean(v), p = p,
                                                              rank_among_KOs = rank(-df, na.last = "keep")[k], n_KOs = sum(!is.na(df))) } }
res <- do.call(rbind, rows); write.csv(res, file.path(OUT, "spec53_module_KO_pairs_agreement.csv"), row.names = FALSE)
options(width = 220); cat(sprintf("\nlines used: %s (removed: %s); H2052 bulk = same-line reference (1 NTC)\n", paste(USE, collapse = ", "), paste(EXCL, collapse = ", ")))
for (i in seq_len(nrow(PAIRS))) { r <- res[res$KO == PAIRS$KO[i] & res$module == PAIRS$MOD[i], ]; cat(sprintf("\n=== %s in %s (NTC - KO; positive = higher in NTC) ===\n", PAIRS$KO[i], PAIRS$MOD[i]))
  print(format(r[, c("source", "n_NTC", "n_KO", "mean_NTC", "mean_KO", "NTC_minus_KO", "p", "rank_among_KOs", "n_KOs")], digits = 3), row.names = FALSE)
  bl <- r[r$source %in% USE, ]; scr <- r[r$source == "CROP-seq (H2052)", ]
  cat(sprintf("  sign agreement with screen: %d/%d retained lines (%s); H2052 bulk: %s\n", sum(sign(bl$NTC_minus_KO) == sign(scr$NTC_minus_KO)), nrow(bl), paste(bl$source[sign(bl$NTC_minus_KO) == sign(scr$NTC_minus_KO)], collapse = ", "),
              if (sign(r$NTC_minus_KO[r$source == "NCI-H2052"]) == sign(scr$NTC_minus_KO)) "agrees" else "disagrees"))
  d <- do.call(rbind, lapply(USE, function(s) { x <- getx(s, PAIRS$MOD[i]); x <- x[x$g %in% c("NTC", PAIRS$KO[i]), ]; x$line <- s; x })); d$g <- relevel(factor(d$g), ref = PAIRS$KO[i])
  if (length(USE) > 1) { m <- summary(lm(x ~ line + g, d))$coefficients["gNTC", ]; cat(sprintf("  bulk pooled over %s (line-adjusted lm): NTC - KO = %+.3f (p=%.3g)\n", paste(USE, collapse = "+"), m[1], m[4])) } }

save_fig("fig3K_module_KO_pairs_agreement", 12.5, 3.2 * nrow(PAIRS) + 1.2, function() {
  layout(matrix(seq_len(nrow(PAIRS) * length(SRC)), nrow = nrow(PAIRS), byrow = TRUE)); par(oma = c(0, 0, 3.2, 0), family = "sans")
  for (i in seq_len(nrow(PAIRS))) for (s in SRC) { k <- PAIRS$KO[i]; m <- PAIRS$MOD[i]; d <- getx(s, m); d <- d[d$g %in% c("NTC", k), ]; d$g <- factor(d$g, levels = c("NTC", k)); par(mar = c(3.6, 4.2, 2.8, 0.5))
    boxplot(x ~ g, d, outline = FALSE, col = c("grey80", adjustcolor(KO_COL[[k]], 0.55)), border = INK2, xlab = "", ylab = paste(m, "score"), ylim = range(d$x))
    set.seed(1); big <- s == "CROP-seq (H2052)"; stripchart(x ~ g, d, vertical = TRUE, add = TRUE, method = "jitter", jitter = 0.12, pch = 16, cex = if (big) 0.15 else 1, col = if (big) adjustcolor(INK, 0.3) else INK)
    r <- res[res$KO == k & res$module == m & res$source == s, ]; scr <- res[res$KO == k & res$module == m & res$source == "CROP-seq (H2052)", ]
    agree <- s != "CROP-seq (H2052)" && sign(r$NTC_minus_KO) == sign(scr$NTC_minus_KO)
    mtext(sprintf("%s in %s: %s%s", k, m, s, if (s == "NCI-H2052") " (1 NTC, ref.)" else ""), side = 3, adj = 0, line = 1.3, font = 2, cex = 0.72)
    mtext(sprintf("NTC-KO %+.3f%s | rank %d/%d%s", r$NTC_minus_KO, if (is.finite(r$p)) sprintf(", p=%.2g", r$p) else "", r$rank_among_KOs, r$n_KOs, if (s == "CROP-seq (H2052)") "" else if (agree) " | agrees" else " | disagrees"),
          side = 3, adj = 0, line = 0.2, cex = 0.62, col = if (s == "CROP-seq (H2052)") INK2 else if (agree) "#1a9850" else "#d62728") }
  mtext(sprintf("Screen vs bulk agreement for selected KO/module pairs (NTC-KO; %s removed; rank = position of this KO among all KOs by NTC-KO diff, 1 = largest)", paste(EXCL, collapse = ", ")), outer = TRUE, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.9)
})
