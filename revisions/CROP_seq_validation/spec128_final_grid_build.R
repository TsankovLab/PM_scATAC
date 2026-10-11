# Spec 128 (step 1 of 3): inputs of the final combination grid (see spec128_final_grid_run.R and spec129_final_grid_figures.R).
# sample set: S70 without H2052_NTC_R3 (carries BPTF guides), with the clean H2052 controls H28_NTC_R1-R3, plus the RUNX2 and SNAI2 KO samples; RUNX1 excluded (no RUNX1 guides in its samples)
# built here, per cell line (all samples of the line, controls included):
#   normalisations  tmm (TMM logCPM) | cp10k (log2 CP10K + 1) | tmmq (TMM logCPM, quantile-normalised over all samples of the set) | vst (DESeq2, per line)
#   corrections     none | cc | tech (log depth, genes detected, % mito, liver-gene signal) | decont (carry-over within the i5 group removed from the counts, RHO 0.7%) | ruv2 (replicate-noise directions)
#                   | cc+tech | decont+cc+ruv2            (i5-block centring is applied later, at the profile level)
#   module scores   AddModuleScore of the 20 Cm (top-50 and top-25 genes) on every normalisation x correction matrix
# per sample: editing evidence (spec118 counts: best-guide indel %, call), own-guide reads per million unmapped (spec119), i5 group; per KO x line: best-guide indel % of the KO samples
# screen (cells of NTC + the 8 shared KOs): group means of log-normalised expression as is and with S.Score + G2M.Score regressed out over those cells; vst ranking of the genes; Cm scores per cell
#   (top-50 and top-25, as is and cell cycle regressed)
# needs: spec72 (S70 lists), spec118, spec119, sample sheet, agree_common.R::load_bulk_S91, spec_common.R
# output: spec128_inputs.rds
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages({ library(Seurat); library(limma) }); options(Seurat.object.assay.version = "v3", width = 220); RHO <- 0.007; KOS10 <- c(KOS, "RUNX2", "SNAI2")
pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); keep <- pr$genes; B <- load_bulk_S91(pr$samples, drop_bad_ntc = TRUE, drop_kos = "RUNX1"); m <- B$meta; ids <- B$ids
ss <- readLines(file.path(PROC, "sample_sheet/SampleSheet.csv")); ss <- ss[(grep("^\\[BCLConvert_Data\\]", ss) + 2):length(ss)]; ss <- do.call(rbind, strsplit(ss[nzchar(ss)], ",")); i5 <- setNames(ss[, 4], ss[, 2])
cat("samples:", length(ids), "\n"); print(table(factor(m$Cell_Line, LINES), factor(m$Gene_Control, c("NTC", KOS10))))

## ---------------------------------------------------------------- sample annotation: technical covariates, editing evidence, guide reads
ct <- read.delim(file.path(PROC, "counts/gene_counts_matrix_stranded_reverse.tsv"), row.names = 1, check.names = FALSE); sy <- read.delim(file.path(PROC, "counts/gene_id_to_symbol.tsv"), header = FALSE, col.names = c("id", "symbol")); sy <- sy[sy$id %in% rownames(ct) & !is.na(sy$symbol) & sy$symbol != "", ]
csA <- rowsum(as.matrix(ct[sy$id, ]), group = sy$symbol)[rownames(B$cs), ]; dep <- colSums(csA); LIVG <- intersect(c("ALB", "APOA1", "APOA2", "APOB", "FGB", "FGG", "TF", "TTR", "SERPINA1", "AFP"), rownames(csA))
SA <- data.frame(sample = ids, line = m$Cell_Line, KO = m$Gene_Control, i5 = i5[ids], cc = B$cc[ids], depth = log10(dep[ids]), genes = colSums(csA[, ids] > 0), mito = 100 * colSums(csA[grepl("^MT-", rownames(csA)), ids]) / dep[ids], liver = log10(1e6 * colSums(csA[LIVG, ids]) / dep[ids] + 1), stringsAsFactors = FALSE); rownames(SA) <- ids
C <- read.csv(file.path(OUT, "spec118_site_counts.csv"), stringsAsFactors = FALSE); St <- read.csv(file.path(OUT, "spec118_gRNA_sites.csv")); St <- St[!is.na(St$chr), ]; C <- C[paste(C$guide, C$site) %in% paste(St$guide, sprintf("%s:%d", St$chr, St$cut)), ]; C$gene <- sub("_g.*", "", C$guide)
ev <- t(sapply(ids, function(s) { x <- C[C$sample == s & C$gene == SA[s, "KO"], ]; if (!nrow(x)) return(c(pairs = NA, best = NA, ng = NA)); c(pairs = sum(x$n_cover), best = max(ifelse(x$n_cover >= 5, 100 * x$n_indel / x$n_cover, 0)), ng = sum(x$n_indel >= 2 & x$n_indel / pmax(1, x$n_cover) >= 0.02)) }))
SA$indel_best <- ev[, "best"]; SA$cut_call <- ifelse(SA$KO == "NTC", "control", ifelse(is.na(ev[, "ng"]), "not assessable", ifelse(ev[, "ng"] >= 1, "cut", ifelse(ev[, "pairs"] >= 40, "no evidence", "not assessable"))))
gr <- read.csv(file.path(OUT, "spec119_guide_reads.csv"), stringsAsFactors = FALSE); gr$gene <- sub("_g.*", "", gr$guide); SA$guide_rpm <- sapply(ids, function(s) { x <- gr[gr$sample == s & gr$gene == SA[s, "KO"], ]; if (!nrow(x)) NA else 1e6 * sum(x$reads) / x$unmapped_reads[1] })
gsum <- read.csv(file.path(OUT, "spec118_guide_summary.csv"), stringsAsFactors = FALSE); gsum <- gsum[gsum$KO_pairs >= 10, ]; IND <- matrix(NA_real_, 10, 4, dimnames = list(KOS10, LINES)); for (i in seq_len(nrow(gsum))) { r <- gsum$target[i]; cl <- gsum$line[i]; if (r %in% KOS10 && (is.na(IND[r, cl]) || gsum$KO_indel_pct[i] > IND[r, cl])) IND[r, cl] <- gsum$KO_indel_pct[i] }
cat("\ncut call:\n"); print(table(SA$cut_call, SA$line)); cat("KO samples with own-guide reads < 200 per million:", SA$sample[SA$KO != "NTC" & !is.na(SA$guide_rpm) & SA$guide_rpm < 200], "\nbest-guide indel % per KO x line:\n"); print(round(IND, 1))

## ---------------------------------------------------------------- bulk matrices: normalisation x correction, per line
allS <- colnames(csA); gt <- sapply(unique(i5[allS]), function(g) rowSums(csA[, i5[allS] == g, drop = FALSE])); gd <- colSums(gt); dec <- sapply(ids, function(s) { g <- i5[[s]]; pmax(csA[, s] - RHO * dep[[s]] * (gt[, g] - csA[, s]) / (gd[[g]] - dep[[s]]), 0) }); rownames(dec) <- rownames(csA)
reg <- function(y, X) { X <- scale(as.matrix(X), scale = FALSE); X <- X[, apply(X, 2, sd) > 0, drop = FALSE]; y - y %*% X %*% solve(crossprod(X)) %*% t(X) }
ruv <- function(y, gl, k) { R <- y; for (kk in unique(gl)) { w <- gl == kk; R[, w] <- y[, w, drop = FALSE] - rowMeans(y[, w, drop = FALSE]) }; R <- R[, gl %in% names(which(table(gl) > 1)), drop = FALSE]; U <- svd(R, nu = k, nv = 0)$u; y - U %*% crossprod(U, y - rowMeans(y)) }
ln <- function(cnt, s, nm) switch(nm, tmm = normalise(cnt[keep, s], "tmm")$y, cp10k = normalise(cnt[keep, s], "cp10k")$y, vst = normalise(cnt[keep, s], "vst")$y)
qn <- function(L) { A <- do.call(cbind, L); Q <- limma::normalizeQuantiles(A); dimnames(Q) <- dimnames(A); lapply(L, function(x) Q[, colnames(x), drop = FALSE]) }
SL <- lapply(setNames(LINES, LINES), function(cl) ids[m$Cell_Line == cl]); BASE <- list(); for (nm in c("tmm", "cp10k", "vst")) for (src in c("raw", "dec")) BASE[[nm]][[src]] <- lapply(SL, function(s) ln(if (src == "raw") B$cs else dec, s, nm))
BASE$tmmq <- list(raw = qn(BASE$tmm$raw), dec = qn(BASE$tmm$dec)); NORMS <- c("tmm", "cp10k", "tmmq", "vst"); CORRS <- c("none", "cc", "tech", "decont", "ruv2", "cc+tech", "decont+cc+ruv2"); MAT <- list()
for (nm in NORMS) for (cl in LINES) { s <- SL[[cl]]; gl <- SA[s, "KO"]; y <- BASE[[nm]]$raw[[cl]]; yd <- BASE[[nm]]$dec[[cl]]; cc <- SA[s, "cc", drop = FALSE]; tech <- SA[s, c("depth", "genes", "mito", "liver")]
  MAT[[nm]][[cl]] <- list(none = y, cc = reg(y, cc), tech = reg(y, tech), decont = yd, ruv2 = ruv(y, gl, 2), `cc+tech` = reg(y, cbind(cc, tech)), `decont+cc+ruv2` = ruv(reg(yd, cc), gl, 2)) }
cat("\nbulk matrices built:", length(NORMS), "normalisations x", length(CORRS), "corrections x 4 lines;", nrow(MAT$tmm[[1]]$none), "genes\n")

## ---------------------------------------------------------------- module scores on every matrix
top50 <- cm_top50(); CM <- names(top50); top25 <- lapply(top50, head, 25)
ms <- function(y, sets) { y <- y - min(y); sb <- SetAssayData(CreateSeuratObject(counts = round(pmax(y, 0))), slot = "data", new.data = y); sb <- suppressWarnings(AddModuleScore(sb, features = sets, name = "Cm_", seed = 1)); x <- t(as.matrix(sb@meta.data[colnames(y), paste0("Cm_", seq_along(sets))])); rownames(x) <- names(sets); x }
MOD <- list(); for (nm in NORMS) for (cl in LINES) for (cr in CORRS) { y <- MAT[[nm]][[cl]][[cr]]; MOD[[nm]][[cl]][[cr]] <- list(top50 = ms(y, top50), top25 = ms(y, top25)) }; cat("module scores done\n")

## ---------------------------------------------------------------- screen
srt <- load_screen(); grp <- as.character(srt$merged_call); GR <- c("NTC", KOS); ck <- which(grp %in% GR); Y <- GetAssayData(srt, layer = "data")[, ck]; g <- grp[ck]; md <- srt@meta.data[ck, ]; Y <- Y[Matrix::rowSums(Y > 0) >= 0.01 * ncol(Y), ]; Gs <- intersect(keep, rownames(Y)); Y <- Y[Gs, ]
Xc <- scale(cbind(md$S.Score, md$G2M.Score), scale = FALSE); Bc <- as.matrix(Y %*% Xc) %*% solve(crossprod(Xc)); gm <- sapply(GR, function(k) Matrix::rowMeans(Y[, g == k, drop = FALSE])); gx <- sapply(GR, function(k) colMeans(Xc[g == k, , drop = FALSE]))
hv <- intersect(VariableFeatures(FindVariableFeatures(srt[, ck], selection.method = "vst", nfeatures = nrow(srt), verbose = FALSE)), Gs)
s50 <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds"))$score[colnames(srt)[ck], CM]; o25 <- AddModuleScore(srt[, ck], features = top25, name = "Cm_", seed = 1); s25 <- as.matrix(o25@meta.data[, paste0("Cm_", 1:20)]); colnames(s25) <- CM
rs <- function(S) apply(S, 2, function(v) resid(lm(v ~ Xc)) + mean(v)); SCS <- list(top50 = list(raw = s50, cc = rs(s50)), top25 = list(raw = s25, cc = rs(s25)))
cat("screen: cells", ncol(Y), "| genes shared with bulk", length(Gs), "\n"); print(table(g)[GR])
saveRDS(list(SA = SA, IND = IND, KOS10 = KOS10, NORMS = NORMS, CORRS = CORRS, MAT = MAT, MOD = MOD, genes = keep, Gs = Gs, hv = hv, screen = list(raw = gm, cc = gm - Bc %*% gx), screen_groups = g, screen_scores = SCS, CM = CM), file.path(OUT, "spec128_inputs.rds"))
cat("saved spec128_inputs.rds\n")
