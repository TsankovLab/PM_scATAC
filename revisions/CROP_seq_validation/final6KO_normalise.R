# final6KO_normalise.R: the expression processing of spec195 (sections 1-3, default TMM log-CPM), from the counts to every normalised matrix, written as tables. Same cohort
# (TEAD4, TWIST1, PITX1, SOX9, TEAD2, BPTF + NTC; MSTO-211H, NCI-H2452, NCI-H28; 43 samples) and the same code, so the tables are what the final6KO figures use.
#   bulk    STAR gene counts (load_bulk_S91) -> genes passing filterByExpr in every line -> per line: TMM log-CPM (edgeR)
#           -> cell-cycle score S.Score + G2M.Score (Seurat CellCycleScoring on the line's samples) -> regressed out within the line -> module scores
#           (AddModuleScore, top-50 genes of the 20 Cm programs) on the regressed values -> KO profiles = mean of the 6 KO means - KO mean (positive = lower in the KO)
#           -> proliferation z per line = (cc - mean NTC cc) / SD
#   screen  NCI-H2052 CROP-seq cells of NTC + the 6 KOs (unique guide) -> CellCycleScoring -> genes detected in >= 1% of the cells and kept in bulk
#           -> KO means with S/G2M regressed out -> KO profiles as in bulk; vst variable genes; per-cell module scores, cell cycle regressed
# needs ($CROP_OUT): spec128_inputs.rds (SA), spec72_profiles_S70.rds (sample list); STAR counts, the screen object, scRNA_PM_modules.rds
# outputs ($CROP_OUT/final6KO_normalised/): samples.csv (annotation + cc and proliferation z), counts.csv.gz (cohort counts, all genes),
#   genes_kept.txt; tmm_logcpm.csv.gz and tmm_logcpm_ccregressed.csv.gz (genes x samples, each line normalised on its own samples);
#   tmm_logcpm_module_scores.csv (Cm x samples); tmm_logcpm_KO_profiles_<line>.csv and tmm_logcpm_KO_profiles_pooled.csv (genes x KO, lines averaged);
#   screen_KO_means_ccregressed.csv.gz, screen_KO_profiles.csv.gz, screen_variable_genes.txt (vst order), screen_module_scores_by_cell.csv.gz, screen_cells.csv
# check: recomputes page 1 C (proliferation p) and page 1 D (match to the screen) of spec195 from these tables, to compare with final6KO_stats.csv
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages({ library(Seurat); library(limma); library(edgeR); library(data.table) }); options(Seurat.object.assay.version = "v3", width = 200); set.seed(5)
KOS <- c("TEAD4", "TWIST1", "PITX1", "SOX9", "TEAD2", "BPTF"); GR <- c("NTC", KOS); LN <- c("MSTO-211H", "NCI-H2452", "NCI-H28"); SH <- function(l) sub("NCI-|-211H", "", l); TOPN <- 500
ND <- file.path(OUT, "final6KO_normalised"); dir.create(ND, showWarnings = FALSE)
wt <- function(M, f) fwrite(data.table(gene = rownames(M), as.data.frame(M, check.names = FALSE)), file.path(ND, f))

## =============================================================== 1. samples and counts (cohort only) -- as spec195
SA <- readRDS(file.path(OUT, "spec128_inputs.rds"))$SA; SA <- SA[SA$line %in% LN & SA$KO %in% GR, c("sample", "line", "KO", "indel_best")]; SA <- SA[order(match(SA$line, LN), match(SA$KO, GR), SA$sample), ]; rownames(SA) <- SA$sample
stopifnot(nrow(SA) == 43, all(table(SA$KO)[KOS] > 0)); cat("cohort:", nrow(SA), "samples\n"); print(table(factor(SA$KO, GR), factor(SA$line, LN)))
pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); CS <- load_bulk_S91(pr$samples, drop_bad_ntc = TRUE, drop_kos = "RUNX1")$cs[, SA$sample]; CS <- CS[rowSums(CS) > 0, ]
lines_of <- function(cl, set = GR) SA$sample[SA$line == cl & SA$KO %in% set]
wt(CS, "counts.csv.gz")

## =============================================================== 2. bulk processing -- as spec195 (TMM log-CPM)
reg <- function(y, X) { X <- scale(as.matrix(X), scale = FALSE); X <- X[, apply(X, 2, sd) > 0, drop = FALSE]; y - y %*% X %*% solve(crossprod(X)) %*% t(X) }
ccs <- function(y) { y <- y - min(y); rownames(y) <- gsub("_", "-", rownames(y)); o <- CreateSeuratObject(counts = round(y)); o <- suppressWarnings(SetAssayData(o, layer = "data", new.data = y))
  o <- suppressWarnings(CellCycleScoring(o, s.features = intersect(Seurat::cc.genes$s.genes, rownames(y)), g2m.features = intersect(Seurat::cc.genes$g2m.genes, rownames(y)), set.ident = FALSE)); setNames(as.numeric(o$S.Score + o$G2M.Score), colnames(y)) }
ms <- function(y, sets) { y <- y - min(y); sb <- suppressWarnings(SetAssayData(CreateSeuratObject(counts = round(pmax(y, 0))), layer = "data", new.data = y)); sb <- suppressWarnings(AddModuleScore(sb, features = sets, name = "Cm_", seed = 1))
  x <- t(as.matrix(sb@meta.data[colnames(y), paste0("Cm_", seq_along(sets))])); rownames(x) <- names(sets); x }
kprof <- function(y, s) { gm <- sapply(KOS, function(k) rowMeans(y[, s[SA[s, "KO"] == k], drop = FALSE])); rowMeans(gm) - gm }
keep <- Reduce(intersect, lapply(LN, function(cl) { s <- lines_of(cl); rownames(CS)[filterByExpr(DGEList(CS[, s]), group = SA[s, "KO"])] })); cat("bulk genes (filterByExpr in every line):", length(keep), "\n")
writeLines(keep, file.path(ND, "genes_kept.txt")); top50 <- cm_top50(); CM <- names(top50)
B <- list()
for (norm in "tmm_logcpm") { Y <- Yc <- Mo <- list(); cc <- setNames(numeric(nrow(SA)), SA$sample)
  for (cl in LN) { s <- lines_of(cl); y <- normalise(CS[keep, s], "tmm")$y; cc[s] <- ccs(y)[s]
    Y[[cl]] <- y; Yc[[cl]] <- reg(y, cc[s]); Mo[[cl]] <- ms(Yc[[cl]], top50) }
  z <- cc; for (cl in LN) { s <- lines_of(cl); z[s] <- (cc[s] - mean(cc[s[SA[s, "KO"] == "NTC"]])) / sd(cc[s]) }
  SA[[paste0("cc_", norm)]] <- cc[SA$sample]; SA[[paste0("prolif_z_", norm)]] <- z[SA$sample]
  wt(do.call(cbind, Y)[, SA$sample], paste0(norm, ".csv.gz")); wt(do.call(cbind, Yc)[, SA$sample], paste0(norm, "_ccregressed.csv.gz")); wt(do.call(cbind, Mo)[, SA$sample], paste0(norm, "_module_scores.csv"))
  KP <- lapply(setNames(LN, LN), function(cl) kprof(Yc[[cl]], lines_of(cl, KOS))); for (cl in LN) wt(KP[[cl]], sprintf("%s_KO_profiles_%s.csv", norm, SH(cl)))
  wt(Reduce(`+`, KP) / length(KP), paste0(norm, "_KO_profiles_pooled.csv")); B[[norm]] <- list(Yc = Yc, KP = KP) }
fwrite(SA, file.path(ND, "samples.csv"))

## =============================================================== 3. screen processing (cohort cells only) -- as spec195
srt <- load_screen(); grp <- as.character(srt$merged_call); srt <- srt[, grp %in% GR]; srt <- srt[Matrix::rowSums(GetAssayData(srt, layer = "counts") > 0) > 0, ]; g <- as.character(srt$merged_call)
srt <- suppressWarnings(CellCycleScoring(srt, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes, set.ident = FALSE))
Y <- GetAssayData(srt, layer = "data"); Y <- Y[Matrix::rowSums(Y > 0) >= 0.01 * ncol(Y), ]; Gs <- intersect(keep, rownames(Y)); Y <- Y[Gs, ]
Xc <- scale(cbind(srt$S.Score, srt$G2M.Score), scale = FALSE); Bc <- as.matrix(Y %*% Xc) %*% solve(crossprod(Xc)); gm <- sapply(GR, function(k) Matrix::rowMeans(Y[, g == k, drop = FALSE])); gx <- sapply(GR, function(k) colMeans(Xc[g == k, , drop = FALSE]))
SCR <- gm - Bc %*% gx; Sp <- rowMeans(SCR[, KOS]) - SCR[, KOS]; hv <- intersect(VariableFeatures(FindVariableFeatures(srt, selection.method = "vst", nfeatures = nrow(srt), verbose = FALSE)), Gs)
o <- suppressWarnings(AddModuleScore(srt, features = top50, name = "Cm_", seed = 1)); s50 <- as.matrix(o@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(s50) <- CM; s50 <- apply(s50, 2, function(v) resid(lm(v ~ Xc)) + mean(v))
wt(SCR, "screen_KO_means_ccregressed.csv.gz"); wt(Sp, "screen_KO_profiles.csv.gz"); writeLines(hv, file.path(ND, "screen_variable_genes.txt")); wt(s50, "screen_module_scores_by_cell.csv.gz")
fwrite(data.table(cell = colnames(srt), group = g, S.Score = srt$S.Score, G2M.Score = srt$G2M.Score), file.path(ND, "screen_cells.csv"))
cat("screen cells:", ncol(srt), "| genes:", length(Gs), "\n")

## =============================================================== check: spec195 page 1 C and D from these tables
for (norm in names(B)) { z <- SA[[paste0("prolif_z_", norm)]]; names(z) <- SA$sample; ntc <- z[SA$KO == "NTC"]
  klp <- aggregate(z ~ KO + line, data.frame(z = z, SA)[SA$KO != "NTC", ], mean); pmC <- sapply(KOS, function(k) median(klp$z[klp$KO == k])); pC <- suppressWarnings(wilcox.test(pmC, ntc, alternative = "less")$p.value)
  vg <- Gs[order(-rowMeans(sapply(LN, function(cl) matrixStats::rowVars(B[[norm]]$Yc[[cl]][Gs, ]))))]; f <- intersect(hv[1:TOPN], vg[1:TOPN]); PB <- Reduce(`+`, lapply(B[[norm]]$KP, function(x) x[Gs, ])) / length(LN)
  CD <- cor(Sp[f, ], PB[f, ], method = "spearman"); same <- diag(CD); other <- CD[row(CD) != col(CD)]
  cat(sprintf("check %s: proliferation p %.4g | match to screen: same %.3f vs other %.3f, %d genes, p %.4g\n", norm, pC, mean(same), mean(other), length(f), wilcox.test(same, other, alternative = "greater")$p.value)) }
cat("written:", ND, "\n")
