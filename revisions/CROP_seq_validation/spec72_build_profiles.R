# Spec 72: NTC - KO profiles of every KO in every line, for every combination of normalisation x covariates x statistic, on gene level and on the 20 Cm modules.
#   sample set   : S70 (baseline, load_bulk() screen KOs) | S70N4 (S70 + the 3 NTC reassigned to NCI-H2052: H28_NTC_R1/2/3) | S75 (all 76 minus the 0.36 M-read H2452_TWIST1_R2)
#                  | RUN3 (S70 minus RUN3_EXCLUDE, 60 samples; build with SETS=RUN3)
#   normalisation: cp10k (log2 CP10K+1, legacy) | tmm (edgeR TMM logCPM) | vst (DESeq2) | voomw (limma voom with sample quality weights: noisy samples are down-weighted, not dropped)
#   covariates   : none | cc | cc_pa (cc + % reads assigned to genes) | cc_pa_depth (+ log10 assigned reads)      [% assigned is the main technical factor in the PCA, spec71]
#   statistic    : med (median NTC - median KO, legacy) | coef (limma) | t (moderated t) | t_rest (moderated t of mean(other KOs) - KO; NTC not used)
# Genes: filterByExpr in every line (intersection). Cm scores: Seurat::AddModuleScore, top-50 genes per module (as spec45), over all samples of the set.
# needs: spec71_sample_audit_76.R (spec71_sample_audit.csv), git_repo/files/scRNA_PM_modules.rds
# outputs: spec72_profiles_<set>.rds (genes, P = gene profiles, M = module profiles, gene_var, gene_signal, cm_scores, sample weights, n per group)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3")
NORMS <- c("cp10k", "tmm", "vst", "voomw"); COVS <- c("none", "cc", "cc_pa", "cc_pa_depth"); STATS <- c("med", "coef", "t", "t_rest")
fab <- load_full_audited_bulk(); meta <- fab$meta; cs <- fab$cs; q <- load_audit(); top50 <- cm_top50(); CM <- names(top50)
SETS <- list(S70 = q$sample[q$in_S70], S70N4 = q$sample[q$in_S70 | q$sample %in% c("H28_NTC_R1", "H28_NTC_R2", "H28_NTC_R3")], S75 = q$sample[!q$low_depth],
             RUN3 = q$sample[q$in_S70 & !q$sample %in% RUN3_EXCLUDE])   # RUN3: S70 minus RUN3_EXCLUDE (agree_common.R); not built by default
SETS <- SETS[strsplit(Sys.getenv("SETS", unset = "S70,S70N4,S75"), ",")[[1]]]

for (sn in names(SETS)) { ids <- SETS[[sn]]; m <- meta[ids, ]; message("== ", sn, ": ", length(ids), " samples")
  keep <- Reduce(intersect, lapply(LINES, function(cl) { s <- ids[m$Cell_Line == cl]; rownames(cs)[filterByExpr(DGEList(cs[, s]), group = m[s, "Gene_Control"])] }))
  message("genes expressed in all 4 lines: ", length(keep))
  # module scores over all samples of the set (cp10k and TMM logCPM input)
  score <- function(y) { y <- y[rowSums(cs[, ids]) > 0, ]; y <- y - min(y)     # undetected genes dropped, scale shifted to >= 0 (AddModuleScore bins genes by mean; scores are differences, so unchanged)
    sb <- CreateSeuratObject(counts = cs[rownames(y), ids]); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(y)); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1))
    x <- t(as.matrix(sb@meta.data[ids, paste0("Cm_", seq_along(CM))])); rownames(x) <- CM; x }
  cm_scores <- list(cp10k = score(log2(sweep(cs[, ids], 2, colSums(cs[, ids]), "/") * 1e4 + 1)), tmm = score(cpm(calcNormFactors(DGEList(cs[, ids])), log = TRUE, prior.count = 2)))
  arr <- function(f) array(NA_real_, c(length(f), length(KOS), length(LINES)), dimnames = list(f, KOS, LINES))
  P <- M <- list(); SW <- list(); gene_var <- gene_signal <- matrix(NA, length(keep), length(LINES), dimnames = list(keep, LINES))
  for (nm in NORMS) for (cv in COVS) { P[[nm]][[cv]] <- setNames(lapply(STATS, function(s) arr(keep)), STATS); M[[nm]][[cv]] <- setNames(lapply(STATS, function(s) arr(CM)), STATS) }
  for (cl in LINES) { s <- ids[m$Cell_Line == cl]; grp <- m[s, "Gene_Control"]; qq <- q[s, ]
    for (nm in NORMS) { nz <- normalise(cs[keep, s], nm)
      if (nm == "tmm") gene_var[, cl] <- matrixStats::rowVars(nz$y)
      for (cv in COVS) { X <- covar_matrix(qq, cv); if (!is.null(X) && ncol(nz$y) - length(unique(grp)) - ncol(X) < 2) X <- X[, seq_len(max(1, ncol(nz$y) - length(unique(grp)) - 2)), drop = FALSE]
        e <- ko_effects(nz$y, grp, X, d = nz$d, weights = if (nm == "voomw") "voomw" else NULL); sw <- attr(e, "sample_weights")
        for (st in STATS) P[[nm]][[cv]][[st]][, , cl] <- e[[st]]
        if (nm == "tmm" && cv == "none") gene_signal[, cl] <- matrixStats::rowMaxs(abs(e$t), na.rm = TRUE)
        if (nm == "voomw" && cv == "cc_pa") SW[[cl]] <- sw
        ms <- if (nm == "cp10k") cm_scores$cp10k[, s] else cm_scores$tmm[, s]       # vst / voomw: TMM-based scores (voomw adds the sample quality weights)
        em <- ko_effects(ms, grp, X, weights = if (nm == "voomw") as.numeric(sw) else NULL)
        for (st in STATS) M[[nm]][[cv]][[st]][, , cl] <- em[[st]] } }
    message("  ", cl, ": ", length(s), " samples, NTC n=", sum(grp == "NTC")) }
  saveRDS(list(set = sn, samples = ids, genes = keep, P = P, M = M, gene_var = gene_var, gene_signal = gene_signal, cm_scores = cm_scores, sample_weights = SW,
               n = table(m$Cell_Line, m$Gene_Control)), file.path(OUT, sprintf("spec72_profiles_%s%s.rds", sn, DSUF)))
  cat(sprintf("\n%s voom sample quality weights (cc_pa design), lowest 8:\n", sn)); print(round(head(sort(unlist(unname(SW))), 8), 2)) }
