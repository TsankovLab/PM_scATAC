#!/usr/bin/env Rscript
# R1_Q3_SOX9_chromVAR_motif_cooccurrence.R
#
# Which ArchR/chromVAR motifs co-occur in the same peaks as SOX9?
# Uses the ArchR motif-match matrix (Annotations/Motif-Matches-In-Peaks.rds,
# peaks x motifs binary) from the tumor_compartment SOX9 P23 project.
#
# For every motif M: peak-level co-occurrence with SOX9 across the full peakset
#   contingency (SOX9+ x M+), hypergeometric enrichment p + odds ratio.

suppressMessages({ library(SummarizedExperiment); library(Matrix); library(data.table) })

P   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_P23"
OUT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
QUERY <- "SOX9"

message("Loading Motif-Matches-In-Peaks.rds ...")
obj <- readRDS(file.path(P, "Annotations", "Motif-Matches-In-Peaks.rds"))
M <- tryCatch(assays(obj)[["matches"]], error = function(e) assay(obj, 1))
M <- as(M, "CsparseMatrix")                 # peaks x motifs, logical/numeric
storage.mode_ok <- TRUE
mot <- colnames(M)
n   <- nrow(M)
message("peaks=", n, "  motifs=", ncol(M))

qcol <- grep(QUERY, mot, value = TRUE, ignore.case = TRUE)
message("Query motif column(s): ", paste(qcol, collapse = ", "))
stopifnot(length(qcol) >= 1)
sox9 <- as.numeric(rowSums(M[, qcol, drop = FALSE]) > 0)   # SOX9+ peak indicator
s <- sum(sox9)
message(QUERY, "+ peaks: ", s, " (", round(100 * s / n, 2), "%)")

# vectorized co-occurrence: a = #peaks with SOX9 & motif j
a <- as.numeric(crossprod(M, sox9))         # length = n motifs
m <- as.numeric(colSums(M))                 # motif+ peaks
# 2x2: a (both), m-a (M only), s-a (SOX9 only), n-m-s+a (neither)
OR <- (a * (n - m - s + a)) / ((m - a) * (s - a))
# hypergeometric enrichment P(X >= a)
pval <- phyper(a - 1, m, n - m, s, lower.tail = FALSE)

res <- data.table(
  motif            = mot,
  n_motif_peaks    = as.integer(m),
  cooccur_SOX9     = as.integer(a),
  pct_of_SOX9peaks = round(100 * a / s, 2),      # % of SOX9+ peaks that also have M
  pct_of_Mpeaks    = round(100 * a / m, 2),      # % of M+ peaks that also have SOX9
  baseline_pct     = round(100 * m / n, 2),      # genome-wide M+ peak rate
  odds_ratio       = round(OR, 3),
  hyper_p          = pval
)
res <- res[!motif %in% qcol]                      # drop SOX9 self
res[, padj := p.adjust(hyper_p, "BH")]
setorder(res, -odds_ratio)
fwrite(res, file.path(OUT, "SOX9_chromVAR_motif_cooccurrence.csv"))
message("Wrote ", file.path(OUT, "SOX9_chromVAR_motif_cooccurrence.csv"))

cat("\n#### TOP 30 motifs co-occurring with SOX9 (by odds ratio; padj<0.05) ####\n")
print(head(res[padj < 0.05 & odds_ratio > 1][order(-odds_ratio)], 30))
cat("\n#### TOP 20 by % of SOX9+ peaks also containing the motif ####\n")
print(head(res[order(-pct_of_SOX9peaks)], 20)[,
       .(motif, n_motif_peaks, cooccur_SOX9, pct_of_SOX9peaks, baseline_pct, odds_ratio, padj = signif(padj,2))])
