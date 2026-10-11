# Spec 94: CROP-seq vs bulk NCI-H2052 without NTC subtraction. Profiles per KO:
#   NTC - KO (reference)   : screen pseudobulk logCPM NTC - KO; bulk limma coef NTC - KO                                 (uncentred)
#   raw KO                 : screen pseudobulk logCPM of the KO; bulk mean TMM logCPM of the KO samples                     (uncentred, no NTC)
#   KO - mean(other KOs)   : the same raw KO profiles centred over the 8 KOs (= NTC - KO centred: NTC cancels exactly)
#   t vs other KOs         : screen Welch t of KO cells vs the cells of the other 7 KOs; bulk limma moderated t of mean(other KOs) - KO (t_rest)  [noise-weighted, no NTC]
# genes: screen top-5000 HVGs expressed in both | top 1000 / 3000 bulk variable genes | all shared genes
# metric: KO x KO Spearman, mean own - mean other, exact permutation p of the bulk KO labels (8!), own-KO rank per TF
# needs: spec88 (pseudobulk), spec72 (S70), spec_common.R (load_screen())
# outputs: spec94_H2052_no_NTC.csv, spec94_H2052_no_NTC_per_KO.csv, fig9P_H2052_no_NTC_subtraction.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages({ library(Seurat); library(ComplexHeatmap); library(circlize); library(grid) })
CL <- "NCI-H2052"; pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); L <- as.matrix(read.csv(file.path(OUT, "spec88_screen_pseudobulk_logCPM.csv"), row.names = 1, check.names = FALSE))
fab <- load_full_audited_bulk(); meta <- fab$meta; s <- pr$samples[meta[pr$samples, "Cell_Line"] == CL]; g <- meta[s, "Gene_Control"]; yb <- normalise(fab$cs[pr$genes, s], "tmm")$y
srt <- load_screen(); Y <- GetAssayData(srt, layer = "data"); grp <- as.character(srt$merged_call)
lc_all <- { X <- GetAssayData(srt, layer = "counts"); PB <- sapply(c("NTC", KOS), function(k) Matrix::rowSums(X[, grp == k, drop = FALSE])); edgeR::cpm(edgeR::DGEList(PB), log = TRUE, prior.count = 2) }
Gall <- Reduce(intersect, list(pr$genes, rownames(lc_all), rownames(Y)))
SPACES <- list(hvg = intersect(rownames(L), Gall), top1000 = Gall[order(pr$gene_var[Gall, CL], decreasing = TRUE)][1:1000], top3000 = Gall[order(pr$gene_var[Gall, CL], decreasing = TRUE)][1:3000], all = Gall)
SPL <- c(hvg = "screen top-5000 HVGs", top1000 = "top 1000 H2052 variable", top3000 = "top 3000 H2052 variable", all = "all shared genes")
welch_rest <- function(f) { Yf <- Y[f, ]; cells <- grp %in% KOS; sapply(KOS, function(k) { a <- which(grp == k); b <- which(cells & grp != k)
  ma <- Matrix::rowMeans(Yf[, a]); mb <- Matrix::rowMeans(Yf[, b]); va <- (Matrix::rowSums(Yf[, a]^2) - length(a) * ma^2) / (length(a) - 1); vb <- (Matrix::rowSums(Yf[, b]^2) - length(b) * mb^2) / (length(b) - 1)
  t <- (mb - ma) / sqrt(pmax(va, 1e-8) / length(a) + pmax(vb, 1e-8) / length(b)); t[!is.finite(t)] <- 0; t }) }        # mean(other KOs) - KO, like bulk t_rest
bulk_raw <- sapply(KOS, function(k) rowMeans(yb[, g == k, drop = FALSE]))
cen <- function(P) P - rowMeans(P)
perms <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(1:8)) })
VIEWS <- c(ntc = "NTC - KO (uncentred)", raw = "raw KO expression (uncentred)", rawc = "KO - mean(other KOs)", trest = "t: KO vs other KOs")
res <- list(); pk <- list()
for (sp in names(SPACES)) { f <- SPACES[[sp]]; tr <- welch_rest(f)
  for (v in names(VIEWS)) { P <- switch(v, ntc = list(S = SGN * (lc_all[f, "NTC"] - lc_all[f, KOS]), B = pr$P$tmm$none$coef[f, , CL]), raw = list(S = lc_all[f, KOS], B = bulk_raw[f, ]),
                                         rawc = list(S = cen(lc_all[f, KOS]), B = cen(bulk_raw[f, ])), trest = list(S = tr, B = pr$P$tmm$none$t_rest[f, , CL]))
    C <- cor(P$S, P$B, method = "spearman"); dm <- mean(diag(C)); nul <- apply(perms, 1, function(p) mean(C[cbind(1:8, p)]))
    res[[length(res) + 1]] <- data.frame(space = sp, view = v, diag = dm, off = mean(C[row(C) != col(C)]), perm_p = mean(nul >= dm - 1e-12), top1 = sum(apply(C, 1, which.max) == 1:8))
    pk[[length(pk) + 1]] <- data.frame(space = sp, view = v, KO = KOS, own = diag(C), median_other = sapply(1:8, function(i) median(C[i, -i])), rank = sapply(1:8, function(i) rank(-C[i, ])[i])) } }
res <- do.call(rbind, res); res$dmo <- res$diag - res$off; pk <- do.call(rbind, pk)
write.csv(res, file.path(OUT, "spec94_H2052_no_NTC.csv"), row.names = FALSE); write.csv(pk, file.path(OUT, "spec94_H2052_no_NTC_per_KO.csv"), row.names = FALSE)
options(width = 200); print(res[, c("space", "view", "diag", "off", "dmo", "perm_p", "top1")], digits = 3, row.names = FALSE)
for (v in names(VIEWS)) { x <- pk[pk$view == v, ]; cat("\n--", VIEWS[[v]], ": own rho (rank)\n"); print(reshape(transform(x[, c("KO", "space")], o = sprintf("%.2f (%d)", x$own, x$rank)), idvar = "KO", timevar = "space", direction = "wide"), row.names = FALSE) }

VC <- c(ntc = "grey55", raw = "#eda100", rawc = "#2a78d6", trest = "#1baf7a")
save_fig("fig9P_H2052_no_NTC_subtraction", 15, 8.6, function() {
  grid.newpage(); pushViewport(viewport(layout = grid.layout(3, 4, heights = unit(c(1.2, 1, 1.1), c("cm", "null", "null")))))
  grid.text("CROP-seq vs bulk NCI-H2052: KO profiles with and without NTC subtraction", x = unit(4, "mm"), y = unit(0.65, "npc"), just = "left", gp = gpar(fontsize = 11, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4))
  grid.text("KO x KO Spearman over genes; top: own KO vs other KOs (bars = mean own and mean other), exact permutation p of bulk KO labels; bottom: own-KO rho per TF, # = rank of own KO among 8 (screen HVGs)", x = unit(4, "mm"), y = unit(0.2, "npc"), just = "left", gp = gpar(fontsize = 7, col = INK2), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4))
  for (j in seq_along(SPACES)) { sp <- names(SPACES)[j]; d <- res[res$space == sp, ]; d <- d[match(names(VIEWS), d$view), ]
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = j)); pushViewport(plotViewport(c(4.8, 3.6, 2.4, 0.8))); yr <- range(c(0, d$diag, d$off)) + c(-0.03, 0.08)
    pushViewport(dataViewport(xscale = c(0.4, 4.6), yscale = yr)); grid.yaxis(gp = gpar(fontsize = 6.5)); grid.lines(unit(c(0.4, 4.6), "native"), unit(c(0, 0), "native"), gp = gpar(col = INK2))
    grid.rect(unit(1:4 - 0.17, "native"), unit(0, "native"), unit(0.3, "native"), unit(d$diag, "native"), just = c("centre", "bottom"), gp = gpar(fill = VC[d$view], col = NA))
    grid.rect(unit(1:4 + 0.17, "native"), unit(0, "native"), unit(0.3, "native"), unit(d$off, "native"), just = c("centre", "bottom"), gp = gpar(fill = adjustcolor(VC[d$view], 0.3), col = VC[d$view]))
    grid.text(sprintf("p %.2g\n%d/8", d$perm_p, d$top1), unit(1:4, "native"), unit(pmax(d$diag, d$off, 0), "native") + unit(4, "mm"), gp = gpar(fontsize = 6, fontface = ifelse(d$perm_p < 0.05, "bold", "plain")))
    grid.text(c("NTC - KO", "raw KO", "KO - other\nKOs", "t vs other\nKOs"), unit(1:4, "native"), unit(-1.6, "lines"), gp = gpar(fontsize = 6.2))
    grid.text(SPL[[sp]], y = unit(1, "npc") + unit(0.9, "lines"), gp = gpar(fontsize = 8, fontface = "bold")); if (j == 1) grid.text("Spearman (solid = own KO, light = other KOs)", x = unit(-2.7, "lines"), rot = 90, gp = gpar(fontsize = 6.5)); popViewport(3) }
  for (j in seq_along(VIEWS)) { v <- names(VIEWS)[j]; x <- pk[pk$view == v & pk$space == "hvg", ]; x <- x[match(KOS, x$KO), ]
    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = j)); pushViewport(plotViewport(c(3.6, 3.4, 2.4, 0.8))); yr <- c(min(-0.3, x$own), max(0.6, x$own) + 0.1)
    pushViewport(dataViewport(xscale = c(0.4, 8.6), yscale = yr)); grid.yaxis(gp = gpar(fontsize = 6.5)); grid.lines(unit(c(0.4, 8.6), "native"), unit(c(0, 0), "native"), gp = gpar(col = INK2))
    grid.rect(unit(1:8, "native"), unit(0, "native"), unit(0.7, "native"), unit(x$own, "native"), just = c("centre", "bottom"), gp = gpar(fill = KO_COL[KOS], col = NA))
    grid.points(unit(1:8, "native"), unit(x$median_other, "native"), pch = "-", size = unit(5, "mm"), gp = gpar(col = INK))
    grid.text(sprintf("#%d", x$rank), unit(1:8, "native"), unit(pmax(x$own, 0), "native") + unit(2.4, "mm"), gp = gpar(fontsize = 6, fontface = ifelse(x$rank == 1, "bold", "plain")))
    grid.text(KOS, unit(1:8, "native"), unit(-1, "lines"), rot = 45, just = "right", gp = gpar(fontsize = 6.2, col = KO_COL[KOS], fontface = "bold"))
    grid.text(VIEWS[[v]], y = unit(1, "npc") + unit(0.9, "lines"), gp = gpar(fontsize = 7.5, fontface = "bold", col = VC[[v]])); popViewport(3) }
})
