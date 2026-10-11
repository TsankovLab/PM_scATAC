# Spec 71b (step 3 of 3): cell-line identity of every bulk RNA-seq sample from its SNP genotypes (independent of expression, so not affected by the KO or by batch).
# Genotype of a sample at a site = 0/1/2 (ref/het/alt) from the allele depths, only where the sample has >= MIN_DP reads (default 10): alt fraction < 0.1 -> 0, > 0.9 -> 2, 0.25-0.75 -> 1,
# otherwise no call. Pairwise concordance = fraction of sites with the same genotype, over sites called in both samples and variable in the cohort. Reference genotype of a line = majority
# genotype of the samples NAMED after that line, leaving the tested sample out. Assigned line = the reference the sample agrees with most.
# A mixture of two lines shows as a low best concordance and a shifted alt-allele fraction at sites where its own line is homozygous (contamination column).
# needs: spec71b_genotype_identity.lsf (calls_GT_AD_DP.tsv.gz, calls_samples.txt)
# outputs: spec71_genotype_concordance.csv (per sample), spec71_genotype_pairwise.csv (matrix), fig7E_genotype_identity.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R"))
suppressPackageStartupMessages(library(data.table))
GENO_DIR <- Sys.getenv("GENO_DIR", unset = "/sc/arion/scratch/giottb01/CROPseq_validation_bulkRNA/work/genotype"); MIN_DP <- as.integer(Sys.getenv("MIN_DP", unset = "10"))
smp <- readLines(file.path(GENO_DIR, "calls_samples.txt")); x <- fread(file.path(GENO_DIR, "calls_GT_AD_DP.tsv.gz"), header = FALSE, sep = "\t")
REF <- ALT <- matrix(NA_real_, nrow(x), length(smp), dimnames = list(NULL, smp))
for (j in seq_along(smp)) { p <- tstrsplit(x[[4 + j]], ":", fixed = TRUE); a <- tstrsplit(p[[2]], ",", fixed = TRUE); REF[, j] <- suppressWarnings(as.numeric(a[[1]])); ALT[, j] <- suppressWarnings(as.numeric(a[[2]])) }
DP <- REF + ALT; AF <- ALT / DP; G <- matrix(NA_integer_, nrow(AF), ncol(AF), dimnames = dimnames(AF)); ok <- !is.na(DP) & DP >= MIN_DP
G[ok & AF < 0.1] <- 0L; G[ok & AF > 0.9] <- 2L; G[ok & AF >= 0.25 & AF <= 0.75] <- 1L
keep <- rowSums(!is.na(G)) >= 0.5 * ncol(G) & apply(G, 1, function(g) length(unique(g[!is.na(g)])) > 1); G <- G[keep, ]; AF <- AF[keep, ]; DP <- DP[keep, ]
message(nrow(G), " informative sites; median sites called per sample: ", median(colSums(!is.na(G))))

PFX <- c(H2052 = "NCI-H2052", MSTO = "MSTO-211H", H2452 = "NCI-H2452", H28 = "NCI-H28"); named <- PFX[sub("_.*", "", sub("^MIXED_", "", smp))]; names(named) <- smp
conc <- function(a, b) { i <- !is.na(a) & !is.na(b); c(mean(a[i] == b[i]), sum(i)) }
PW <- sapply(smp, function(a) sapply(smp, function(b) conc(G[, a], G[, b])[1])); write.csv(round(PW, 4), file.path(OUT, "spec71_genotype_pairwise.csv"))
majority <- function(M) apply(M, 1, function(g) { g <- g[!is.na(g)]; if (length(g) < 3) return(NA_integer_); t <- table(g); if (max(t) / length(g) < 0.7) return(NA_integer_); as.integer(names(t)[which.max(t)]) })
res <- do.call(rbind, lapply(smp, function(s) { cc <- sapply(LINES, function(l) { ref <- majority(G[, setdiff(smp[named == l], s), drop = FALSE]); conc(G[, s], ref)[1] }); best <- names(which.max(cc))
  ref <- majority(G[, setdiff(smp[named == best], s), drop = FALSE]); hom <- !is.na(ref) & ref != 1 & DP[, s] >= 20 & !is.na(DP[, s])
  dev <- ifelse(ref[hom] == 0, AF[hom, s], 1 - AF[hom, s])        # fraction of reads carrying the allele the assigned line does not have
  data.frame(sample = s, named_line = named[[s]], genotype_line = best, conc_best = max(cc), conc_second = sort(cc, decreasing = TRUE)[2], t(setNames(cc, paste0("conc_", sub("NCI-|-211H", "", LINES)))),
             n_sites = sum(!is.na(G[, s])), foreign_allele_fraction = mean(dev), frac_sites_foreign_gt5pct = mean(dev > 0.05)) }))
res$agrees_with_name <- res$named_line == res$genotype_line; aud <- tryCatch(read.csv(file.path(OUT, "spec71_sample_audit.csv")), error = function(e) NULL)
if (!is.null(aud)) { res$line_in_analysis <- aud$line[match(res$sample, aud$sample)]; res$in_S70 <- aud$in_S70[match(res$sample, aud$sample)]; res$agrees_with_analysis_line <- res$line_in_analysis == res$genotype_line }
res$possible_mixture <- res$frac_sites_foreign_gt5pct > 3 * median(res$frac_sites_foreign_gt5pct) | res$conc_best < 0.95
write.csv(res, file.path(OUT, "spec71_genotype_concordance.csv"), row.names = FALSE)
options(width = 250); cat("\n=== named line vs genotype line ===\n"); print(table(named = res$named_line, genotype = res$genotype_line))
cat("\n=== samples whose genotype disagrees with their NAME, or flagged as possible mixture ===\n"); print(res[!res$agrees_with_name | res$possible_mixture, ], digits = 3, row.names = FALSE)
if (!is.null(aud)) { cat("\n=== samples of the 76 whose genotype disagrees with the line used in the ANALYSIS ===\n"); print(res[which(!res$agrees_with_analysis_line), ], digits = 3, row.names = FALSE) }

save_fig("fig7E_genotype_identity", 11, 10, function() {
  o <- order(match(res$genotype_line, LINES), res$named_line, res$sample); M <- PW[o, o]; n <- nrow(M); par(mar = c(1, 9, 4, 1), family = "sans")
  image(seq_len(n), seq_len(n), t(M)[, n:1], col = colorRampPalette(c("#f0efec", "#2a78d6", "#0b2d5c"))(100), zlim = c(min(M, na.rm = TRUE), 1), axes = FALSE, xlab = "", ylab = "")
  lab <- rownames(M); bad <- !res$agrees_with_name[o]; axis(2, at = n:1, labels = FALSE, tick = FALSE)
  text(par("usr")[1] - 0.5, n:1, lab, adj = 1, xpd = NA, cex = 0.42, col = ifelse(bad, "#c0392b", INK2), font = ifelse(bad, 2, 1))
  rect(par("usr")[1] - 0.2, (n:1) - 0.5, par("usr")[1], (n:1) + 0.5, col = LINE_COL[res$named_line[o]], border = NA, xpd = NA)
  b <- cumsum(table(factor(res$genotype_line[o], LINES))); abline(v = b + 0.5, h = n - b + 0.5, col = "white", lwd = 2)
  text((c(0, b[-4]) + b) / 2, n + 1.5, LINES, xpd = NA, cex = 0.8, font = 2)
  mtext("SNP genotype concordance between bulk RNA-seq samples (blocks = line assigned by genotype)", side = 3, adj = 0, line = 2.6, font = 2, cex = 0.95)
  mtext(sprintf("%d expressed exonic SNPs, >= %d reads; side strip = line the sample is NAMED after; bold red = name disagrees with genotype", nrow(G), MIN_DP), side = 3, adj = 0, line = 1.6, cex = 0.6, col = INK2)
})
