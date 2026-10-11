# Shared helpers for the cross-line agreement analysis (spec72-76). Sourced after repro_common.R and profiles.R.
# Direction: every profile is NTC - KO (positive = lower after KO), DIFF_DIR=KO_NTC flips the sign (SGN).
suppressPackageStartupMessages({ library(edgeR); library(limma) })
DIFF_DIR <- Sys.getenv("DIFF_DIR", unset = "NTC_KO"); SGN <- if (DIFF_DIR == "NTC_KO") 1 else -1
DLAB <- if (SGN == 1) "NTC - KO" else "KO - NTC"; DSUF <- if (SGN == 1) "" else "_KOminusNTC"
KOS <- ORDER; CM_RDS <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds"
PAIRS <- combn(LINES, 2, simplify = FALSE); PLAB <- sapply(PAIRS, function(p) paste(sub("NCI-|-211H", "", p), collapse = "-"))
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

# sample audit table of spec71 (cc, % assigned, depth per sample) and the sample sets
load_audit <- function() { q <- read.csv(file.path(OUT, "spec71_sample_audit.csv"), stringsAsFactors = FALSE); rownames(q) <- q$sample; q }
# RUN3 (2026-09-29): S70 minus samples removed after the batch / QC review (spec77 PCA + audit): whole KO groups where both replicates look off (MSTO BPTF, MSTO HMGA1,
# H2452 TEAD4) and single replicates. With these removed, BPTF and HMGA1 are absent from MSTO-211H and TEAD4 from NCI-H2452 (NA in the profiles, see ko_effects / agreement).
RUN3_EXCLUDE <- c("MSTO_BPTF_R1", "MSTO_BPTF_R3", "MSTO_HMGA1_R1", "MSTO_HMGA1_R2", "MSTO_SOX9_R3", "MSTO_TWIST1_R2", "H2452_BPTF_R1", "H2452_TWIST1_R2", "H2452_TEAD4_R1", "H2452_TEAD4_R3",
                  "MIXED_H28_TCF3_R1_Huh7", "H28_TEAD4_R3", "H28_HMGA1_R2", "H28_NTC_R1", "H28_NTC_R2", "H28_NTC_R3")
cm_top50 <- function() { m <- readRDS(CM_RDS); m <- lapply(m, head, 50); m[order(as.integer(sub("Cm", "", names(m))))] }

# covariate matrix for one line (columns centred); spec = "none" | "cc" | "cc_pa" | "cc_pa_depth"
covar_matrix <- function(q, spec) { if (spec == "none") return(NULL)
  X <- switch(spec, cc = cbind(cc = q$cc), cc_pa = cbind(cc = q$cc, pa = q$pct_assigned), cc_pa_depth = cbind(cc = q$cc, pa = q$pct_assigned, depth = log10(q$assigned_M)))
  scale(X, scale = FALSE) }

# normalised expression of one line: list(y = genes x samples, d = DGEList (TMM) for voom). counts = symbol-level counts of the line's samples, already gene-filtered
normalise <- function(counts, norm) { d <- calcNormFactors(DGEList(counts))
  y <- switch(norm, cp10k = log2(sweep(counts, 2, colSums(counts), "/") * 1e4 + 1), tmm = cpm(d, log = TRUE, prior.count = 2), voomw = cpm(d, log = TRUE, prior.count = 2),
              vst = { suppressPackageStartupMessages(library(DESeq2)); x <- round(counts); storage.mode(x) <- "integer"; SummarizedExperiment::assay(vst(DESeqDataSetFromMatrix(x, data.frame(a = rep(1, ncol(x))), ~ 1), blind = TRUE)) })
  list(y = y, d = d) }

# NTC - KO effect of every KO in one line. y = features x samples (genes, or module scores); returns list of features x KO matrices:
#   med  = median(NTC) - median(KO) after removing the covariates (group effects protected)      [legacy statistic]
#   coef = limma coefficient NTC - KO;  t = moderated t of the same contrast
#   t_rest = moderated t of (mean of the OTHER KOs) - KO: NTC is not used, so the noise of the shared NTC does not enter
# weights: NULL | "voomw" (voomWithQualityWeights on d; genes only) | numeric sample weights (module scores)
ko_effects <- function(y, group, X = NULL, d = NULL, weights = NULL, kos = KOS) {
  all_kos <- kos; kos <- kos[kos %in% group]                                        # KOs with no sample in this line (RUN3) are returned as NA columns
  if (length(kos) < length(all_kos)) { e <- ko_effects(y, group, X, d, weights, kos); pad <- function(m) { o <- matrix(NA_real_, nrow(m), length(all_kos), dimnames = list(rownames(m), all_kos)); o[, kos] <- m; o }
    out <- lapply(e[c("med", "coef", "t", "t_rest", "coef_rest")], pad); attr(out, "sample_weights") <- attr(e, "sample_weights"); return(out) }
  g <- factor(group, levels = c("NTC", kos)); G <- model.matrix(~ 0 + g); colnames(G) <- levels(g); des <- cbind(G, X)
  yr <- if (is.null(X)) y else removeBatchEffect(y, covariates = X, design = G)
  md <- function(k) matrixStats::rowMedians(yr[, g == k, drop = FALSE])
  med <- sapply(kos, function(k) md("NTC") - md(k)); rownames(med) <- rownames(y)
  v <- if (identical(weights, "voomw")) voomWithQualityWeights(d, des, plot = FALSE) else NULL
  fit <- if (!is.null(v)) lmFit(v, des) else if (is.numeric(weights)) lmFit(y, des, weights = matrix(weights, nrow(y), ncol(y), byrow = TRUE)) else lmFit(y, des)
  cn <- function(L) { f <- eBayes(contrasts.fit(fit, L)); list(coef = f$coefficients, t = f$t) }
  L1 <- sapply(kos, function(k) { v <- setNames(numeric(ncol(des)), colnames(des)); v["NTC"] <- 1; v[k] <- -1; v })
  L2 <- sapply(kos, function(k) { v <- setNames(numeric(ncol(des)), colnames(des)); v[setdiff(kos, k)] <- 1 / (length(kos) - 1); v[k] <- -1; v })
  a <- cn(L1); b <- cn(L2); sw <- if (!is.null(v)) setNames(v$targets$sample.weights, colnames(y)) else NULL
  out <- list(med = med, coef = a$coef, t = a$t, t_rest = b$t, coef_rest = b$coef); out <- lapply(out, function(m) { dimnames(m) <- list(rownames(y), kos); SGN * m }); attr(out, "sample_weights") <- sw; out }

# KO x KO correlation between two lines; Pa, Pb = features x KO. centre: subtract each feature's mean over the KOs (removes what all KOs share, incl. the NTC term);
# zscore: also scale each feature across the KOs. Done per line.
prep <- function(P, centre, zscore) { if (centre || zscore) P <- P - rowMeans(P, na.rm = TRUE); if (zscore) { s <- matrixStats::rowSds(P, na.rm = TRUE); P <- P / ifelse(s > 0, s, Inf) }; P }
pair_cor <- function(Pa, Pb, method) { C <- suppressWarnings(cor(Pa, Pb, method = method)); C[is.na(C)] <- 0
  C[colSums(!is.na(Pa)) == 0, ] <- NA; C[, colSums(!is.na(Pb)) == 0] <- NA; C }                                 # KO absent from a line (RUN3): NA, not 0

# agreement from the list of 6 KO x KO matrices (rows = KOs of line a, columns = KOs of line b). Per KO: intra = the 6 same-KO values, cross = that KO vs every other KO
# (both directions). p-values: one-sided Wilcoxon intra > cross (pairs share samples -> descriptive), and a label permutation: KO labels permuted independently in each line,
# statistic = mean same-label correlation over the pairs (per KO), and mean over KOs (global). perms = list (one per line of LINES) of n_perm x 8 permuted indices (make_perms());
# pairs = the line pairs CL refers to (default all 6).
make_perms <- function(n, K = length(KOS), L = length(LINES), seed = 1) { set.seed(seed); lapply(seq_len(L), function(l) t(replicate(n, sample.int(K)))) }
# NA-aware: a KO absent from a line (RUN3) has NA in the pairs involving that line and is averaged over the pairs where it exists (identical results when there is no NA).
agreement <- function(CL, perms = NULL, pairs = PAIRS) { K <- nrow(CL[[1]]); kos <- rownames(CL[[1]])
  intra <- sapply(CL, diag); rownames(intra) <- kos                                   # KO x pair
  cross <- lapply(seq_len(K), function(k) { v <- unlist(lapply(CL, function(C) c(C[k, -k], C[-k, k]))); v[!is.na(v)] })
  res <- data.frame(KO = kos, n_pairs = rowSums(!is.na(intra)), mean_intra = rowMeans(intra, na.rm = TRUE), median_intra = apply(intra, 1, median, na.rm = TRUE), median_cross = sapply(cross, median), n_pos_pairs = rowSums(intra > 0, na.rm = TRUE),
                    wilcox_p = sapply(seq_len(K), function(k) { v <- intra[k, !is.na(intra[k, ])]; if (length(v) < 2) NA else suppressWarnings(wilcox.test(v, cross[[k]], alternative = "greater")$p.value) }))
  res$sep <- res$median_intra - res$median_cross
  t1 <- function(C) { r <- sapply(seq_len(K), function(k) if (all(is.na(C[k, ])) || is.na(C[k, k])) NA else which.max(C[k, ]) == k); c2 <- sapply(seq_len(K), function(k) if (all(is.na(C[, k])) || is.na(C[k, k])) NA else which.max(C[, k]) == k); mean(c(r, c2), na.rm = TRUE) }
  off <- mean(sapply(CL, function(C) mean(C[row(C) != col(C)], na.rm = TRUE))); dm <- mean(intra, na.rm = TRUE)
  glob <- c(diag_mean = dm, off_mean = off, dmo = dm - off, top1 = mean(sapply(CL, t1)))
  if (!is.null(perms)) { li <- lapply(pairs, function(p) match(p, LINES)); n <- nrow(perms[[1]])
    V <- lapply(seq_along(CL), function(i) { C <- CL[[i]]; a <- perms[[li[[i]][1]]]; b <- perms[[li[[i]][2]]]; matrix(C[cbind(as.vector(a), as.vector(b))], n, K) })   # n_perm x K (label slot) per pair
    S <- Reduce(`+`, lapply(V, function(v) { v[is.na(v)] <- 0; v })); N <- Reduce(`+`, lapply(V, function(v) !is.na(v))); nul <- S / N           # mean over the pairs where the permuted slot exists
    res$perm_p <- sapply(seq_len(K), function(k) { x <- nul[, k]; x <- x[is.finite(x)]; (1 + sum(x >= res$mean_intra[k] - 1e-12)) / (length(x) + 1) }); res$perm_fdr <- p.adjust(res$perm_p, "BH")
    res$perm_z <- (res$mean_intra - colMeans(nul, na.rm = TRUE)) / apply(nul, 2, sd, na.rm = TRUE)
    gn <- rowSums(S) / rowSums(N)
    glob <- c(glob, perm_p = (1 + sum(gn >= dm - 1e-12)) / (n + 1), perm_z = (dm - mean(gn)) / sd(gn)) }
  list(per_KO = res, global = glob, intra = intra, cross = cross) }

# list of KO x KO matrices (one per line pair) from an array A (features x KO x line). feats = NULL (all features) | character vector | "signal" with sig = genes x line matrix of
# KO signal (max |t| over KOs): for every pair the n_sig genes with the largest mean signal in the lines NOT in the pair (selection never sees the two lines being compared).
cl_from_array <- function(A, feats = NULL, sig = NULL, n_sig = 500, centre = TRUE, zscore = FALSE, method = "spearman", pairs = PAIRS, lines_avail = LINES) {
  lapply(pairs, function(p) { f <- if (identical(feats, "signal")) rownames(sig)[order(rowMeans(sig[, setdiff(lines_avail, p), drop = FALSE]), decreasing = TRUE)[seq_len(n_sig)]] else if (is.null(feats)) dimnames(A)[[1]] else feats
    C <- pair_cor(prep(A[f, , p[1]], centre, zscore), prep(A[f, , p[2]], centre, zscore), method); dimnames(C) <- list(dimnames(A)[[2]], dimnames(A)[[2]]); C }) }

# gene-level NTC - KO moderated t for every line from symbol-level counts (TMM logCPM, limma), optional covariates; returns list(t = genes x KO x line, sig = genes x line)
t_profiles <- function(counts, meta, q, ids, genes, cov = "none", lines = LINES) { A <- array(NA_real_, c(length(genes), length(KOS), length(lines)), dimnames = list(genes, KOS, lines)); sig <- matrix(NA_real_, length(genes), length(lines), dimnames = list(genes, lines))
  for (cl in lines) { s <- ids[meta[ids, "Cell_Line"] == cl]; nz <- normalise(counts[genes, s], "tmm"); e <- ko_effects(nz$y, meta[s, "Gene_Control"], covar_matrix(q[s, ], cov)); A[, , cl] <- e$t; sig[, cl] <- matrixStats::rowMaxs(abs(e$t), na.rm = TRUE) }
  list(t = A, sig = sig) }

## ---------------------------------------------------------------- sample set S91 (2026-10-06, after the guide / editing audit of spec118-119)
# H2052_NTC_R3 (the only H2052-genotype NTC of S70) carries BPTF guides and BPTF indels -> relabelled BPTF; the H2052 controls are the 3 H2052-genotype NTCs labelled H28 (no guide reads).
# The RUNX1, RUNX2 and SNAI2 KO samples are added (not in the screen; not in the spec71 audit: cell line taken from the sample name, genotype not verified; MIXED_* samples excluded).
# Caveat (spec119): H2052_RUNX1_R1/R3 carry PITX1 guides and PITX1 indels, MSTO_RUNX1 carry no detectable guide; runx1_h2052_as_pitx1 = TRUE relabels the two H2052 samples PITX1.
S91_NTC_AS_BPTF <- "H2052_NTC_R3"; S91_H2052_NTC <- c("H28_NTC_R1", "H28_NTC_R2", "H28_NTC_R3"); EXTRA_KOS <- c("RUNX1", "RUNX2", "SNAI2"); KOS11 <- c(KOS, EXTRA_KOS)
# returns cs (symbol-level counts), meta (Cell_Line, Gene_Control, label_original), ids, cc (S.Score + G2M.Score scored over the samples of the set). s70 = the S70 sample ids; extra_kos = FALSE keeps the 8 KOs only
# drop_bad_ntc = TRUE removes H2052_NTC_R3 instead of relabelling it; drop_kos removes whole KOs (e.g. "RUNX1": its samples do not carry RUNX1 guides)
load_bulk_S91 <- function(s70, extra_kos = TRUE, runx1_h2052_as_pitx1 = FALSE, drop_bad_ntc = FALSE, drop_kos = character(0)) {
  counts <- read.delim(file.path(PROC, "counts/gene_counts_matrix_stranded_reverse.tsv"), row.names = 1, check.names = FALSE); sym <- read.delim(file.path(PROC, "counts/gene_id_to_symbol.tsv"), header = FALSE, col.names = c("gene_id", "symbol"))
  sym <- sym[sym$gene_id %in% rownames(counts) & !is.na(sym$symbol) & sym$symbol != "", ]; cs <- rowsum(as.matrix(counts[sym$gene_id, ]), group = sym$symbol)
  m0 <- load_full_audited_bulk()$meta; PFX <- c(H2052 = "NCI-H2052", MSTO = "MSTO-211H", H2452 = "NCI-H2452", H28 = "NCI-H28")
  ex <- grep(sprintf("^(H2052|MSTO|H2452|H28)_(%s)_R[0-9]+$", paste(EXTRA_KOS, collapse = "|")), colnames(cs), value = TRUE); ids <- c(s70, S91_H2052_NTC, if (extra_kos) ex)
  meta <- data.frame(row.names = ids, Cell_Line = ifelse(ids %in% rownames(m0), m0[ids, "Cell_Line"], PFX[sub("_.*", "", ids)]), Gene_Control = ifelse(ids %in% rownames(m0), m0[ids, "Gene_Control"], sub("^[^_]+_([^_]+)_.*", "\\1", ids)), stringsAsFactors = FALSE)
  meta$label_original <- meta$Gene_Control; meta[S91_NTC_AS_BPTF, "Gene_Control"] <- "BPTF"; if (runx1_h2052_as_pitx1) meta[rownames(meta) %in% c("H2052_RUNX1_R1", "H2052_RUNX1_R3"), "Gene_Control"] <- "PITX1"
  if (drop_bad_ntc) ids <- setdiff(ids, S91_NTC_AS_BPTF); ids <- ids[!meta[ids, "label_original"] %in% drop_kos]; meta <- meta[ids, ]
  obj <- Seurat::NormalizeData(Seurat::CreateSeuratObject(counts = cs[, ids], min.cells = 0, min.features = 0), verbose = FALSE); obj <- Seurat::CellCycleScoring(obj, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes, set.ident = FALSE)
  list(cs = cs[, ids], meta = meta, ids = ids, cc = setNames(as.numeric(obj$S.Score + obj$G2M.Score), ids)) }
