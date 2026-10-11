# Spec 2: diagnostics for the low KO-KO specificity
#   D1 which bulk line resembles the screen (baseline expression)     -> spec_line_resemblance.csv
#   D2 on-target knockdown in the screen and in bulk                    -> spec_ontarget_screen.csv (bulk: spec_ontarget_knockdown.csv)
#   D3 reliability ceilings: screen split-half and bulk replicate-vs-replicate KO specificity -> spec_reliability.csv
#   D4 technical factors in bulk: variance of the top PCs explained by KO, replicate, index, depth, cell-cycle -> spec_technical_factors.csv
# needs: spec1_prepare_stats.R
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
set.seed(42)
b <- load_bulk(); cp10k <- b$cp10k; meta <- b$meta
scr_g <- as.matrix(read.csv(file.path(OUT, "screen_var5000_log2avg_by_KO.csv"), row.names = 1, check.names = FALSE))
genes <- rownames(scr_g)[rownames(scr_g) %in% rownames(cp10k)]
lg <- log2(cp10k[genes, , drop = FALSE] + 1)

## ---------------------------------------------------------------- D1 line resemblance
d1 <- do.call(rbind, lapply(LINES, function(cl) {
  g <- group_log2avg(cp10k, meta, cl)[genes, , drop = FALSE]
  data.frame(cell_line = cl, r_NTC = cor(scr_g[genes, "NTC"], g[, "NTC"]), r_all_KO_mean = cor(rowMeans(scr_g[genes, colnames(scr_g) != "NTC"]), rowMeans(g[, colnames(g) != "NTC"])),
             spearman_NTC = cor(scr_g[genes, "NTC"], g[, "NTC"], method = "spearman"))
}))
write.csv(d1, file.path(OUT, "spec_line_resemblance.csv"), row.names = FALSE); cat("\nD1 screen vs bulk baseline expression (4127 variable genes):\n"); print(format(d1, digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- D2 on-target knockdown in the screen
ss <- readRDS(file.path(OUT, "spec_screen_stats.rds"))
d2 <- do.call(rbind, lapply(colnames(ss$t0), function(k) if (k %in% ss$genes) data.frame(KO = k, n_cells = ss$n[[k]], lfc_lognorm = ss$diff0[k, k], t = ss$t0[k, k], t_cc_regressed = ss$tcc[k, k]) else data.frame(KO = k, n_cells = ss$n[[k]], lfc_lognorm = NA, t = NA, t_cc_regressed = NA)))
write.csv(d2, file.path(OUT, "spec_ontarget_screen.csv"), row.names = FALSE); cat("\nD2 on-target effect of each KO on its own gene in the screen:\n"); print(format(d2, digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- D3 reliability ceilings
metrics <- function(A, B) {                      # A, B: features x KO, same KO order; row = A profile, column = B profile
  C <- cor(A, B); k <- ncol(C); ranks <- sapply(seq_len(k), function(i) sum(C[i, ] > C[i, i]) + 1)
  Z <- t(scale(t(C)))
  c(n_KO = k, matched_r = mean(diag(C)), other_r = mean(C[row(C) != col(C)]), z_matched = mean(diag(Z)), top1_frac = mean(ranks == 1), median_rank = median(ranks), chance_top1 = 1 / k)
}
centre <- function(M) M - rowMeans(M)
res <- list()
# (a) screen: random halves of the cells of every KO and of NTC (NTC halves differ between A and B)
srt <- load_screen(); Y <- GetAssayData(srt, layer = "data")[genes, ]; grp <- srt$merged_call; kos <- setdiff(sort(unique(grp), method = "radix"), "NTC")
acc <- list()
for (s in 1:100) {
  h <- ave(seq_along(grp), grp, FUN = function(i) sample(rep(1:2, length.out = length(i))))
  key <- paste(grp, h); lev <- unique(key); G <- sparseMatrix(i = seq_along(key), j = match(key, lev), x = 1, dims = c(length(key), length(lev)), dimnames = list(NULL, lev))
  M <- as.matrix(Y %*% G); M <- sweep(M, 2, colSums(G), "/")
  A <- M[, paste(kos, 1)] - M[, "NTC 1"]; B <- M[, paste(kos, 2)] - M[, "NTC 2"]
  acc[[s]] <- rbind(raw = metrics(A, B), centred = metrics(centre(A), centre(B)))
}
avg <- Reduce(`+`, acc) / length(acc); res[["screen split-half (cells)"]] <- data.frame(dataset = "CROP-seq screen split-half (100 splits)", variant = rownames(avg), avg, row.names = NULL)
# (b) bulk: replicate 1 vs replicate 2 of every KO that has two replicates
for (cl in LINES) {
  m <- meta[meta$Cell_Line == cl, ]; ntc <- rowMeans(lg[, rownames(m)[m$Gene_Control == "NTC"], drop = FALSE])
  two <- names(which(table(m$Gene_Control[m$Gene_Control != "NTC"]) == 2)); two <- sort(two, method = "radix")
  A <- sapply(two, function(k) lg[, sort(rownames(m)[m$Gene_Control == k])[1]] - ntc); B <- sapply(two, function(k) lg[, sort(rownames(m)[m$Gene_Control == k])[2]] - ntc)
  r <- rbind(raw = metrics(A, B), centred = metrics(centre(A), centre(B)))
  res[[cl]] <- data.frame(dataset = paste(cl, "bulk replicate 1 vs 2"), variant = rownames(r), r, row.names = NULL)
}
d3 <- do.call(rbind, res); rownames(d3) <- NULL; write.csv(d3, file.path(OUT, "spec_reliability.csv"), row.names = FALSE)
cat("\nD3 reliability: does a KO profile match ITS OWN second half/replicate better than other KOs? (variable genes, KO - NTC)\n"); print(format(d3, digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- D4 technical factors in bulk
ccb <- read.csv(file.path(OUT, "bulk_cellcycle_per_sample_RUN1.csv"))
libsize <- colSums(b$counts)
d4 <- do.call(rbind, lapply(LINES, function(cl) {
  m <- meta[meta$Cell_Line == cl, ]; s <- rownames(m); x <- log2(cp10k[, s] + 1)
  v <- apply(x, 1, var); top <- names(sort(v, decreasing = TRUE))[1:3000]; pc <- prcomp(t(x[top, ]), scale. = FALSE); k <- 10
  f <- data.frame(KO = factor(m$Gene_Control), replicate = factor(sub(".*_R", "R", m$Sample_ID)), i5 = factor(m$i5_index_num), i7 = factor(m$i7_index_num),
                  log_depth = log10(libsize[s]), cell_cycle = ccb$cc[match(s, ccb$sample)])
  # variance-weighted adjusted R^2 of each factor over the top-10 PCs (each factor alone)
  w <- pc$sdev[1:k]^2 / sum(pc$sdev[1:k]^2)
  r2 <- sapply(names(f), function(nm) sum(w * sapply(1:k, function(j) summary(lm(pc$x[, j] ~ f[[nm]]))$adj.r.squared)))
  data.frame(cell_line = cl, n = length(s), t(round(pmax(r2, 0), 3)))
}))
write.csv(d4, file.path(OUT, "spec_technical_factors.csv"), row.names = FALSE)
cat("\nD4 variance-weighted adjusted R2 of the top-10 PCs (top-3000 variable genes) explained by each factor, per cell line:\n"); print(d4, row.names = FALSE)
