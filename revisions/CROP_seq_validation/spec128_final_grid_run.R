# Spec 128 (step 2 of 3): the final combination grid. Inputs from spec128_final_grid_build.R (corrected sample set: clean H2052 controls, no H2052_NTC_R3, RUNX2 / SNAI2 added, RUNX1 excluded).
# factors
#   norm      tmm | cp10k | tmmq (quantile) | vst
#   corr      none | cc | tech | decont | ruv2 | cc+tech | decont+cc+ruv2 | i5block (KO profiles centred within the i5 group of the line)
#   filter    all | cut (samples without evidence of cutting dropped) | indel5 / indel10 (KO x line with best-guide indel >= 5 / 10 %) | bestrep (best-edited replicate of each KO) | guide200 (own-guide reads >= 200 per million)
#   stat      mean (difference of group means) | t (limma moderated t; not with bestrep)
#   contrast  NTC (NTC - KO) | KO (mean of the KOs - KO, no control)
#   scale     none | sd (each gene divided by its SD across the KOs of the line / of the screen)
#   features  top 500 / 1000 / 2000 / 3000 variable genes (E1: in both data sets; E2: of the lines compared, plus the top-500 KO-signal genes of the lines not compared) | Cm modules, top-50 or top-25 genes (stat mean)
#   screen    raw | cc (S.Score + G2M.Score regressed out over the cells)                                        method   spearman | pearson
#   lines     E1: each line, and the bulk pooled over all 4 lines or over 3 (each line left out); E2 / E3: all 4 lines or 3
# endpoints
#   E1  CROP-seq vs bulk (8 shared KOs): KO x KO correlation, own - other, exact permutation p of the bulk labels, own correlation and rank per TF; every KO x KO matrix is kept, and a grid-wide
#       p for the best setting of each target is obtained by applying the same label permutations to all settings and taking the maximum standardised statistic
#   E2  bulk across lines (10 KOs): agreement() -- same-KO vs other-KO correlation between lines, KO-label permutation within lines (N_PERM, default 2000), global and per KO
#   E3  Cm lower in KO than NTC (20 Cm x 8 KOs): NTC - KO of the module score per line (z over the samples of the line), combined over the lines (stratified, pooled within-group variance), with the screen effect
# env: NCORE (default 8), N_PERM (2000), TEST=1 (2 tasks only)
# outputs: spec128_grid_crop.csv, spec128_grid_crop_matrices.rds, spec128_grid_crop_fwer.csv, spec128_grid_lines.csv, spec128_grid_Cm.csv.gz, spec128_screen_Cm.csv
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R")); suppressPackageStartupMessages(library(limma)); options(width = 220)
NCORE <- as.integer(Sys.getenv("NCORE", unset = "8")); N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "2000")); TEST <- Sys.getenv("TEST", unset = "0") == "1"
I <- readRDS(file.path(OUT, "spec128_inputs.rds")); SA <- I$SA; IND <- I$IND; KOS10 <- I$KOS10; Gs <- I$Gs; hv <- I$hv; CM <- I$CM; genes <- I$genes; GR <- c("NTC", KOS)
FILT <- c("all", "cut", "indel5", "indel10", "bestrep", "guide200"); NS <- c(500, 1000, 2000, 3000); POOLS <- c(list(pool_all = LINES), setNames(lapply(LINES, function(l) setdiff(LINES, l)), paste0("pool_no_", sub("NCI-|-211H", "", LINES)))); LSETS <- setNames(POOLS, sub("pool_", "", names(POOLS)))
samples_of <- function(cl, flt) { s <- SA$sample[SA$line == cl]; k <- SA[s, "KO"] != "NTC"; keep <- switch(flt, all = rep(TRUE, length(s)), cut = SA[s, "cut_call"] != "no evidence", indel5 = , indel10 = { th <- if (flt == "indel5") 5 else 10; v <- IND[cbind(match(SA[s, "KO"], KOS10), match(cl, LINES))]; !k | (!is.na(v) & v >= th) },
    bestrep = { e <- SA[s, "indel_best"]; e[is.na(e)] <- -1; b <- ave(e, SA[s, "KO"], FUN = function(v) seq_along(v) == which.max(v)) == 1; !k | b }, guide200 = !k | is.na(SA[s, "guide_rpm"]) | SA[s, "guide_rpm"] >= 200); s[keep | !k] }
# KO profiles of one line: features x kos (NA column = KO absent). y = features x samples
prof <- function(y, gl, kos, stat, contrast) { kp <- kos[kos %in% gl]; P <- matrix(NA_real_, nrow(y), length(kos), dimnames = list(rownames(y), kos)); if (length(kp) < 2) return(P)
  if (stat == "t") { lv <- c("NTC", unique(gl[gl != "NTC"])); g <- factor(gl, levels = lv); des <- model.matrix(~ 0 + g); colnames(des) <- lv; fit <- lmFit(y, des)
    L <- sapply(kp, function(k) { v <- setNames(numeric(length(lv)), lv); if (contrast == "NTC") v["NTC"] <- 1 else v[kp] <- 1 / length(kp); v[k] <- v[k] - 1; v }); P[, kp] <- eBayes(contrasts.fit(fit, L))$t
  } else { gm <- sapply(c("NTC", kp), function(k) rowMeans(y[, gl == k, drop = FALSE])); P[, kp] <- (if (contrast == "NTC") gm[, "NTC"] else rowMeans(gm[, kp, drop = FALSE])) - gm[, kp, drop = FALSE] }
  P }
blockc <- function(P, s, gl) { kp <- colnames(P)[colSums(!is.na(P)) > 0]; blk <- sapply(kp, function(k) names(which.max(table(SA[s[gl == k], "i5"])))); for (b in unique(blk)) { w <- kp[blk == b]; P[, w] <- if (length(w) > 1) P[, w, drop = FALSE] - rowMeans(P[, w, drop = FALSE]) else NA }; P }
sdsc <- function(P) { ok <- colSums(!is.na(P)) > 0; s <- matrixStats::rowSds(P[, ok, drop = FALSE]); P / ifelse(s > 0, s, Inf) }
poolm <- function(PL, lines) { A <- simplify2array(PL[lines]); r <- rowMeans(A, dims = 2, na.rm = TRUE); r[is.nan(r)] <- NA; r }
PMl <- lapply(1:8, function(n) { pf <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(pf(v[-i]), function(r) c(v[i], r)))); do.call(rbind, pf(seq_len(n))) })
e1 <- function(Sp, Bp, f, method) { ok <- which(colSums(!is.na(Bp)) > 0); n <- length(ok); if (n < 4 || length(f) < 15) return(NULL); C <- suppressWarnings(cor(Sp[f, ok], Bp[f, ok], method = method)); if (anyNA(C)) return(NULL)
  dm <- mean(diag(C)); P <- PMl[[n]]; nul <- rowMeans(matrix(C[cbind(as.vector(col(P)), as.vector(P))], nrow(P))); own <- rk <- rep(NA_real_, 8); own[ok] <- diag(C); rk[ok] <- sapply(seq_len(n), function(i) rank(-C[i, ])[i]); Cf <- matrix(NA_real_, 8, 8); Cf[ok, ok] <- C
  list(n = n, dmo = dm - mean(C[row(C) != col(C)]), p = mean(nul >= dm - 1e-12), own = own, rank = rk, C = Cf) }
scr_prof <- function(M, contrast) { M <- M[, GR]; (if (contrast == "NTC") M[, "NTC"] else rowMeans(M[, KOS])) - M[, KOS] }
SCRG <- list(); for (sm in c("raw", "cc")) for (ct in c("NTC", "KO")) SCRG[[sm]][[ct]] <- scr_prof(I$screen[[sm]], ct)
SCRM <- list(); for (tp in c("top50", "top25")) for (sm in c("raw", "cc")) { gmn <- sapply(GR, function(k) colMeans(I$screen_scores[[tp]][[sm]][I$screen_groups == k, , drop = FALSE])); for (ct in c("NTC", "KO")) SCRM[[tp]][[sm]][[ct]] <- scr_prof(gmn, ct) }
permsL <- make_perms(N_PERM, K = length(KOS10))

## ---------------------------------------------------------------- one task = one normalisation x correction
task <- function(nm, cr) { mcr <- if (cr == "i5block") "none" else cr; R1 <- list(); C1 <- list(); R2 <- list(); R3 <- list(); VAR <- sapply(LINES, function(cl) matrixStats::rowVars(I$MAT[[nm]][[cl]][[mcr]])); rownames(VAR) <- genes
  for (flt in FILT) { SS <- lapply(setNames(LINES, LINES), samples_of, flt = flt); GL <- lapply(SS, function(s) SA[s, "KO"]); Yl <- lapply(setNames(LINES, LINES), function(cl) I$MAT[[nm]][[cl]][[mcr]][, SS[[cl]], drop = FALSE])
    # KO signal per line for the E2 gene choice: max |t| of NTC - KO (bestrep: |mean difference| / SD)
    SIG <- sapply(LINES, function(cl) { P <- if (flt == "bestrep") prof(Yl[[cl]], GL[[cl]], KOS10, "mean", "NTC") / pmax(matrixStats::rowSds(Yl[[cl]]), 1e-6) else prof(Yl[[cl]], GL[[cl]], KOS10, "t", "NTC"); matrixStats::rowMaxs(abs(P), na.rm = TRUE) }); rownames(SIG) <- genes; SIG[!is.finite(SIG)] <- 0
    for (st in (if (flt == "bestrep") "mean" else c("mean", "t"))) for (ct in (if (cr == "i5block") "KO" else c("NTC", "KO"))) {
      fin <- function(P, cl) if (cr == "i5block") blockc(P, SS[[cl]], GL[[cl]]) else P
      P8 <- lapply(setNames(LINES, LINES), function(cl) fin(prof(Yl[[cl]][Gs, ], GL[[cl]], KOS, st, ct), cl)); P10 <- lapply(setNames(LINES, LINES), function(cl) fin(prof(Yl[[cl]], GL[[cl]], KOS10, st, ct), cl))
      M8 <- M10 <- NULL; if (st == "mean") { mm <- function(kos) lapply(c(top50 = "top50", top25 = "top25"), function(tp) lapply(setNames(LINES, LINES), function(cl) fin(prof(I$MOD[[nm]][[cl]][[mcr]][[tp]][, SS[[cl]], drop = FALSE], GL[[cl]], kos, "mean", ct), cl))); M8 <- mm(KOS); M10 <- mm(KOS10) }
      ## E1: targets = lines and pools
      TG <- c(setNames(as.list(LINES), LINES), POOLS)
      for (tg in names(TG)) { ls <- TG[[tg]]; vb <- Gs[order(-rowMeans(VAR[Gs, ls, drop = FALSE]))]
        for (sc in c("none", "sd")) { tr <- function(P) if (sc == "sd") sdsc(P) else P; Bg <- if (length(ls) == 1) tr(P8[[ls]]) else poolm(lapply(P8, tr), ls)
          feats <- c(setNames(lapply(NS, function(n) intersect(hv[seq_len(n)], vb[seq_len(n)])), paste0("top", NS)), if (st == "mean") list(Cm50 = CM, Cm25 = CM))
          for (ft in names(feats)) { Bx <- if (ft == "Cm50" || ft == "Cm25") { tp <- if (ft == "Cm50") "top50" else "top25"; if (length(ls) == 1) tr(M8[[tp]][[ls]]) else poolm(lapply(M8[[tp]], tr), ls) } else Bg
            for (sm in c("raw", "cc")) { Sx <- tr(if (ft == "Cm50") SCRM$top50[[sm]][[ct]] else if (ft == "Cm25") SCRM$top25[[sm]][[ct]] else SCRG[[sm]][[ct]])
              for (me in c("spearman", "pearson")) { r <- e1(Sx, Bx, feats[[ft]], me); if (is.null(r)) next
                R1[[length(R1) + 1]] <- c(list(target = tg, norm = nm, corr = cr, filter = flt, stat = st, contrast = ct, scale = sc, feat = ft, screen = sm, method = me, n_KO = r$n, n_feat = length(feats[[ft]]), dmo = r$dmo, p = r$p), setNames(as.list(r$own), paste0("own_", KOS)), setNames(as.list(r$rank), paste0("rank_", KOS))); C1[[length(C1) + 1]] <- r$C } } } } }
      ## E2: across lines
      for (lsn in names(LSETS)) { ls <- LSETS[[lsn]]; prs <- Filter(function(p) all(p %in% ls), PAIRS); vr <- genes[order(-rowMeans(VAR[, ls, drop = FALSE]))]
        for (sc in c("none", "sd")) { tr <- function(P) if (sc == "sd") sdsc(P) else P; Pg <- lapply(P10, tr); feats <- c(list(signal500 = "signal"), setNames(lapply(c(1000, 2000, 3000), function(n) vr[seq_len(n)]), paste0("var", c(1000, 2000, 3000))), if (st == "mean") list(Cm50 = CM, Cm25 = CM))
          for (ft in names(feats)) for (me in c("spearman", "pearson")) { CLs <- lapply(prs, function(p) { if (ft == "signal500") { oth <- setdiff(ls, p); f <- genes[order(-rowMeans(SIG[, oth, drop = FALSE]))[1:500]]; Pa <- Pg[[p[1]]][f, ]; Pb <- Pg[[p[2]]][f, ] } else if (ft %in% c("Cm50", "Cm25")) { tp <- if (ft == "Cm50") "top50" else "top25"; Pa <- tr(M10[[tp]][[p[1]]]); Pb <- tr(M10[[tp]][[p[2]]]) } else { Pa <- Pg[[p[1]]][feats[[ft]], ]; Pb <- Pg[[p[2]]][feats[[ft]], ] }
              Cc <- pair_cor(Pa, Pb, me); dimnames(Cc) <- list(KOS10, KOS10); Cc }); a <- tryCatch(agreement(CLs, permsL, pairs = prs), error = function(e) NULL); if (is.null(a)) next; k <- a$per_KO[match(KOS10, a$per_KO$KO), ]
            R2[[length(R2) + 1]] <- c(list(lines = lsn, norm = nm, corr = cr, filter = flt, stat = st, contrast = ct, scale = sc, feat = ft, method = me, diag_mean = a$global[["diag_mean"]], off_mean = a$global[["off_mean"]], perm_p = a$global[["perm_p"]], perm_z = a$global[["perm_z"]]), setNames(as.list(k$mean_intra), paste0("r_", KOS10)), setNames(as.list(k$perm_p), paste0("p_", KOS10))) } } } }
    ## E3: Cm lower in KO than NTC (not for i5block)
    if (cr != "i5block") for (tp in c("top50", "top25")) { Dl <- array(NA_real_, c(20, 8, 4), dimnames = list(CM, KOS, LINES)); Wl <- SSl <- Nl <- array(0, c(20, 8, 4));
      for (cl in LINES) { S <- I$MOD[[nm]][[cl]][[mcr]][[tp]][, SS[[cl]], drop = FALSE]; z <- (S - rowMeans(S)) / matrixStats::rowSds(S); gl <- GL[[cl]]; a <- z[, gl == "NTC", drop = FALSE]; li <- match(cl, LINES)
        for (k in KOS[KOS %in% gl]) { b <- z[, gl == k, drop = FALSE]; ki <- match(k, KOS); Dl[, ki, li] <- rowMeans(a) - rowMeans(b); Wl[, ki, li] <- ncol(a) * ncol(b) / (ncol(a) + ncol(b)); SSl[, ki, li] <- rowSums((a - rowMeans(a))^2) + rowSums((b - rowMeans(b))^2); Nl[, ki, li] <- ncol(a) + ncol(b) } }
      for (lsn in names(LSETS)) { li <- match(LSETS[[lsn]], LINES); W <- Wl[, , li, drop = FALSE]; D <- Dl[, , li, drop = FALSE]; D0 <- D; D0[is.na(D0)] <- 0; sw <- apply(W, 1:2, sum); est <- apply(W * D0, 1:2, sum) / sw; df <- apply(Nl[, , li, drop = FALSE], 1:2, sum) - 2 * apply(W > 0, 1:2, sum); s2 <- apply(SSl[, , li, drop = FALSE], 1:2, sum) / pmax(df, 1)
        tt <- est / sqrt(s2 / sw); pp <- ifelse(df >= 1, 2 * pt(-abs(tt), pmax(df, 1)), NA); R3[[length(R3) + 1]] <- data.frame(lines = lsn, norm = nm, corr = cr, filter = flt, modules = tp, Cm = rep(CM, 8), KO = rep(KOS, each = 20), est = as.vector(est), p = as.vector(pp), n_lines = as.vector(apply(W > 0, 1:2, sum)), n_lower = as.vector(apply(D > 0, 1:2, sum, na.rm = TRUE)),
                                                                                                           d_H2052 = as.vector(Dl[, , 1]), d_MSTO = as.vector(Dl[, , 2]), d_H2452 = as.vector(Dl[, , 3]), d_H28 = as.vector(Dl[, , 4]), stringsAsFactors = FALSE) } } }
  list(R1 = data.table::rbindlist(R1), C1 = simplify2array(C1), R2 = data.table::rbindlist(R2), R3 = data.table::rbindlist(R3)) }

TASKS <- expand.grid(norm = I$NORMS, corr = c(I$CORRS, "i5block"), stringsAsFactors = FALSE); if (TEST) TASKS <- TASKS[TASKS$norm == "tmm" & TASKS$corr %in% c("none", "cc"), ]; cat("tasks:", nrow(TASKS), "| cores:", NCORE, "| permutations (E2):", N_PERM, "\n"); t0 <- Sys.time()
RES <- parallel::mclapply(seq_len(nrow(TASKS)), function(i) { r <- tryCatch(task(TASKS$norm[i], TASKS$corr[i]), error = function(e) { message("task ", TASKS$norm[i], " ", TASKS$corr[i], ": ", conditionMessage(e)); NULL }); r }, mc.cores = NCORE, mc.preschedule = FALSE)
bad <- sapply(RES, function(r) is.null(r) || inherits(r, "try-error")); cat("tasks failed:", sum(bad), "| elapsed:", format(Sys.time() - t0), "\n"); RES <- RES[!bad]; TAG <- if (TEST) "_test" else ""
G1 <- data.table::rbindlist(lapply(RES, `[[`, "R1")); G1$id <- seq_len(nrow(G1)); CC <- do.call(abind::abind, c(lapply(RES, `[[`, "C1"), list(along = 3))); G2 <- data.table::rbindlist(lapply(RES, `[[`, "R2")); G3 <- data.table::rbindlist(lapply(RES, `[[`, "R3"))
cat("settings: E1", nrow(G1), "| E2", nrow(G2), "| E3", nrow(G3) / 160, "x 160 pairs\n")

## ---------------------------------------------------------------- E1: grid-wide (max-statistic) p per target
set.seed(1); NP <- 2000; PERM <- t(replicate(NP, sample.int(8))); fw <- function(idx) { z <- matrix(NA_real_, length(idx), NP); zo <- numeric(length(idx))
  for (j in seq_along(idx)) { C <- CC[, , idx[j]]; ok <- which(!is.na(diag(C))); n <- length(ok); Ps <- t(apply(PERM[, ok, drop = FALSE], 1, rank)); D <- rowMeans(matrix(C[ok, ok][cbind(rep(seq_len(n), each = NP), as.vector(Ps))], NP)); mu <- mean(C[ok, ok]); sdv <- sd(D); z[j, ] <- (D - mu) / sdv; zo[j] <- (mean(diag(C[ok, ok])) - mu) / sdv }
  list(zo = zo, zmax = matrixStats::colMaxs(z)) }
tgs <- unique(G1$target); FW <- parallel::mclapply(tgs, function(tg) { idx <- G1$id[G1$target == tg]; ch <- split(idx, ceiling(seq_along(idx) / 4000)); r <- lapply(ch, fw); zo <- unlist(lapply(r, `[[`, "zo")); zm <- do.call(pmax, lapply(r, `[[`, "zmax")); list(target = tg, id = idx, zo = zo, zmax = zm) }, mc.cores = min(NCORE, length(tgs)))
G1$z <- NA_real_; for (r in FW) G1$z[match(r$id, G1$id)] <- r$zo; allmax <- do.call(pmax, lapply(FW, `[[`, "zmax"))
FWT <- do.call(rbind, lapply(FW, function(r) { b <- r$id[which.max(r$zo)]; data.frame(target = r$target, n_settings = length(r$id), frac_p05 = mean(G1$p[match(r$id, G1$id)] < 0.05), best_id = b, best_z = max(r$zo), fwer_p_target = mean(r$zmax >= max(r$zo)), fwer_p_all_targets = mean(allmax >= max(r$zo))) }))
data.table::fwrite(G1, file.path(OUT, sprintf("spec128_grid_crop%s.csv", TAG))); saveRDS(CC, file.path(OUT, sprintf("spec128_grid_crop_matrices%s.rds", TAG))); write.csv(FWT, file.path(OUT, sprintf("spec128_grid_crop_fwer%s.csv", TAG)), row.names = FALSE)
data.table::fwrite(G2, file.path(OUT, sprintf("spec128_grid_lines%s.csv", TAG))); data.table::fwrite(G3, file.path(OUT, sprintf("spec128_grid_Cm%s.csv.gz", TAG)))
## screen Cm effects (z over the NTC + KO cells): NTC - KO, Wilcoxon p
SCM <- do.call(rbind, lapply(c("top50", "top25"), function(tp) do.call(rbind, lapply(c("raw", "cc"), function(sm) { S <- scale(I$screen_scores[[tp]][[sm]]); g <- I$screen_groups; do.call(rbind, lapply(KOS, function(k) data.frame(modules = tp, screen = sm, Cm = CM, KO = k, screen_est = colMeans(S[g == "NTC", ]) - colMeans(S[g == k, , drop = FALSE]), screen_p = sapply(CM, function(m) wilcox.test(S[g == "NTC", m], S[g == k, m])$p.value)))) }))))
write.csv(SCM, file.path(OUT, sprintf("spec128_screen_Cm%s.csv", TAG)), row.names = FALSE)
cat("\nE1 calibration and grid-wide p per target:\n"); print(FWT, digits = 3, row.names = FALSE); cat("\nE2: settings with global p < 0.05:", sprintf("%.1f%%", 100 * mean(G2$perm_p < 0.05)), "| median z", round(median(G2$perm_z), 2), "\n")
b <- G1[G1$target == "NCI-H2052" & G1$norm == "tmm" & G1$corr == "none" & G1$filter == "all" & G1$stat == "mean" & G1$contrast == "KO" & G1$scale == "none" & G1$feat == "top2000" & G1$screen == "raw" & G1$method == "spearman", ]; cat("baseline check E1 (H2052, tmm, none, all, mean, KO, top2000, raw, spearman): dmo", round(b$dmo, 3), "p", round(b$p, 2), "| TEAD4 own", round(b$own_TEAD4, 2), "\n")
b <- G2[G2$lines == "all" & G2$norm == "tmm" & G2$corr == "none" & G2$filter == "all" & G2$stat == "t" & G2$contrast == "KO" & G2$scale == "none" & G2$feat == "signal500" & G2$method == "spearman", ]; cat("baseline check E2 (all lines, tmm, none, all, t, KO-centred, signal500, spearman): z", round(b$perm_z, 2), "p", signif(b$perm_p, 2), "\n"); cat("done in", format(Sys.time() - t0), "\n")
