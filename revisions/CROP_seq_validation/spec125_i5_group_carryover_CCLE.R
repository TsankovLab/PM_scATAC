# Spec 125: where do the trace liver transcripts in some i5 index groups come from? Transcriptome and SNV comparison with CCLE, and a test of read carry-over between samples that share an i5 index.
# The run has 106 samples: 4 mesothelioma lines, 7 liver lines (4 pure samples, 3 samples mixed with NCI-H28) and the NR4A2 / SCRAMBLE samples. i5 (index2) defines 5 groups; i7 is unique within a group.
#   1  identity: Spearman of each sample (logCPM) with the 1019 CCLE lines (log2 TPM) over the 3000 most variable CCLE genes -> best match
#   2  excess profile of the new H2052 controls (H28_NTC_R1-R3, liver-containing i5 groups) over the H2052 KO samples (other i5 groups), correlated with the lineage-centred profile of every CCLE line
#   3  marker carry-over: for a recipient sample of mesothelioma line L and a donor line M, markers = genes off in L (CCLE TPM < 0.2) and on in M (CCLE TPM > 20);
#      marker signal (sum of CPM) in recipients that share an i5 group with a sample of M vs recipients that do not; carry-over rate from the 10 classical liver genes = recipient CPM above background / pooled CPM of its i5 group
#   4  guide carry-over (spec119): guide reads of a KO gene in samples of another label, when a KO sample of that gene is in the same i5 group vs not
#   5  SNVs: private somatic SNVs of the 11 lines (CCLE; spec125a): fraction of alternative-allele reads per sample and line (identity), and at donor-line sites in recipients sharing / not sharing the i5 group
# needs: spec125a, spec119, spec71 audit, sample sheet, CCLE_RNAseq_rsem_genes_tpm.rds
# outputs: spec125_identity.csv, spec125_excess_vs_CCLE.csv, spec125_marker_carryover.csv, spec125_guide_carryover.csv, spec125_SNV_fraction.csv, fig42_i5_group_carryover.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R")); options(width = 230)
CCLE <- Sys.getenv("CCLE_DIR", unset = "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Public_data/CCLE"); fp <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2)))
CC <- c(`NCI-H2052` = "NCIH2052_PLEURA", `MSTO-211H` = "MSTO211H_PLEURA", `NCI-H2452` = "NCIH2452_PLEURA", `NCI-H28` = "NCIH28_PLEURA", `JHH-7` = "JHH7_LIVER", `Li-7` = "LI7_LIVER", `SK-HEP-1` = "SKHEP1_LIVER", `SNU-398` = "SNU398_LIVER", `Hep G2` = "HEPG2_LIVER", `Hep 3B` = "HEP3B217_LIVER", `HuH-7` = "HUH7_LIVER")
MESO <- names(CC)[1:4]; LIVER <- names(CC)[5:11]

## ---------------------------------------------------------------- samples: line, i5 group
ct <- read.delim(file.path(PROC, "counts/gene_counts_matrix_stranded_reverse.tsv"), row.names = 1, check.names = FALSE); sy <- read.delim(file.path(PROC, "counts/gene_id_to_symbol.tsv"), header = FALSE, col.names = c("id", "symbol")); sy <- sy[sy$id %in% rownames(ct) & !is.na(sy$symbol) & sy$symbol != "", ]
cs <- rowsum(as.matrix(ct[sy$id, ]), group = sy$symbol); depth <- colSums(cs); cpm <- sweep(cs, 2, depth / 1e6, "/"); ss <- readLines(file.path(PROC, "sample_sheet/SampleSheet.csv")); ss <- ss[(grep("^\\[BCLConvert_Data\\]", ss) + 2):length(ss)]; ss <- do.call(rbind, strsplit(ss[nzchar(ss)], ","))
q <- load_audit(); PFX <- c(H2052 = "NCI-H2052", MSTO = "MSTO-211H", H2452 = "NCI-H2452", H28 = "NCI-H28"); LV <- c(Jhh7 = "JHH-7", Li7 = "Li-7", Skhep1 = "SK-HEP-1", Snu398 = "SNU-398"); MIX <- c(MIXED_H28_SNAI2_R1_HepG2 = "Hep G2", MIXED_H28_SNAI2_R2_HepB3 = "Hep 3B", MIXED_H28_TCF3_R1_Huh7 = "HuH-7")
SM <- data.frame(sample = colnames(cs), i5 = setNames(ss[, 4], ss[, 2])[colnames(cs)], stringsAsFactors = FALSE); SM$line <- ifelse(SM$sample %in% names(LV), LV[SM$sample], ifelse(SM$sample %in% names(MIX), "NCI-H28", PFX[sub("_.*", "", SM$sample)]))
i <- match(SM$sample, q$sample); SM$line[!is.na(i) & !SM$sample %in% names(MIX)] <- q$line[i[!is.na(i) & !SM$sample %in% names(MIX)]]; SM$extra <- ifelse(SM$sample %in% names(MIX), MIX[SM$sample], NA); SM$kind <- ifelse(SM$sample %in% names(LV), "liver", ifelse(SM$sample %in% names(MIX), "mixed", ifelse(is.na(SM$line), "other", "meso")))
rownames(SM) <- SM$sample; cat("samples per i5 group and kind:\n"); print(table(SM$i5, SM$kind)); present <- function(g) unique(na.omit(c(SM$line[SM$i5 == g], SM$extra[SM$i5 == g])))
cat("lines present per i5 group:\n"); for (g in sort(unique(SM$i5))) cat(" ", g, ":", paste(present(g), collapse = ", "), "\n")

## ---------------------------------------------------------------- 1  identity vs CCLE
E <- readRDS(file.path(CCLE, "CCLE_RNAseq_rsem_genes_tpm.rds")); E <- as.matrix(E); G <- intersect(rownames(E), rownames(cs)); L <- log2(E[G, ] + 1); vg <- G[order(matrixStats::rowVars(L), decreasing = TRUE)][1:3000]; lc <- log2(cpm[vg, ] + 1)
R1 <- cor(lc, L[vg, ], method = "spearman"); ID <- data.frame(sample = colnames(lc), line = SM[colnames(lc), "line"], kind = SM[colnames(lc), "kind"], best = colnames(L)[apply(R1, 1, which.max)], rho_best = apply(R1, 1, max))
ID$expected <- CC[ID$line]; ID$rho_expected <- sapply(seq_len(nrow(ID)), function(k) if (is.na(ID$expected[k])) NA else R1[k, ID$expected[k]]); ID$rank_expected <- sapply(seq_len(nrow(ID)), function(k) if (is.na(ID$expected[k])) NA else rank(-R1[k, ])[ID$expected[k]]); write.csv(ID, file.path(OUT, "spec125_identity.csv"), row.names = FALSE)
cat("\n=== 1  best CCLE match ===\n"); x <- ID[ID$kind %in% c("meso", "liver"), ]; print(table(expected = x$line, top_match_is_expected = x$rank_expected == 1)); print(ID[ID$kind != "meso" | ID$rank_expected > 1, c("sample", "line", "best", "rho_best", "rho_expected", "rank_expected")], digits = 3, row.names = FALSE)

## ---------------------------------------------------------------- 2  excess profile of the new H2052 controls vs CCLE
new <- S91_H2052_NTC; koH <- SM$sample[SM$line == "NCI-H2052" & SM$kind == "meso" & !SM$sample %in% new & !grepl("NTC", SM$sample)]; lcA <- log2(cpm[G, ] + 1); ex <- rowMeans(lcA[, new]) - rowMeans(lcA[, koH]); Lc <- L - rowMeans(L)
EX <- data.frame(ccle = colnames(L), rho = as.numeric(cor(ex[vg], Lc[vg, ], method = "spearman")), stringsAsFactors = FALSE); EX$lineage <- sub("^[^_]+_", "", EX$ccle); EX <- EX[order(-EX$rho), ]; EX$rank <- seq_len(nrow(EX)); write.csv(EX, file.path(OUT, "spec125_excess_vs_CCLE.csv"), row.names = FALSE)
cat("\n=== 2  new H2052 controls minus H2052 KO samples: most similar CCLE lines (lineage-centred) ===\n"); print(head(EX, 12), digits = 3, row.names = FALSE); cat("lineages of the top 30:", paste(names(sort(table(EX$lineage[1:30]), decreasing = TRUE)), sort(table(EX$lineage[1:30]), decreasing = TRUE), collapse = ", "), "\n")
cat("rank of the liver lines of this run:", paste(sub("_LIVER", "", CC[LIVER]), EX$rank[match(CC[LIVER], EX$ccle)], collapse = ", "), "of", nrow(EX), "| liver lineage median rank", median(EX$rank[EX$lineage == "LIVER"]), "\n")
top <- names(sort(ex[vg], decreasing = TRUE))[1:40]; cat("top excess genes:", paste(top, collapse = " "), "\n")

## ---------------------------------------------------------------- 3  marker carry-over between lines within an i5 group
rec <- SM[SM$kind == "meso", ]; MK <- list()
for (r in rec$sample) { Lr <- rec[r, "line"]; g <- rec[r, "i5"]
  for (M in setdiff(names(CC), Lr)) { mk <- G[E[G, CC[[Lr]]] < 0.2 & E[G, CC[[M]]] > 20]; if (length(mk) < 5) next; don <- SM$sample[which(SM$i5 == g & SM$sample != r & ((!is.na(SM$line) & SM$line == M) | (!is.na(SM$extra) & SM$extra == M)))]; same <- length(don) > 0
    pool <- if (same) 1e6 * sum(cs[mk, don, drop = FALSE]) / sum(depth[SM$sample[SM$i5 == g]]) else NA                              # donor marker reads per million reads of the whole i5 group
    MK[[length(MK) + 1]] <- data.frame(recipient = r, recipient_line = Lr, i5 = g, donor_line = M, donor_kind = ifelse(M %in% LIVER, "liver", "meso"), n_markers = length(mk), same_group = same, n_donor_samples = length(don), marker_cpm = sum(cpm[mk, r]), pool_cpm = pool) } }
MK <- do.call(rbind, MK); MK$rate <- MK$marker_cpm / MK$pool_cpm; write.csv(MK, file.path(OUT, "spec125_marker_carryover.csv"), row.names = FALSE)
cat("\n=== 3  marker carry-over: marker CPM in recipients sharing / not sharing an i5 group with the donor line ===\n")
T3 <- do.call(rbind, lapply(split(MK, list(MK$recipient_line, MK$donor_line), drop = TRUE), function(d) { a <- d$marker_cpm[d$same_group]; b <- d$marker_cpm[!d$same_group]; if (!length(a) || !length(b)) return(NULL)
  data.frame(recipient_line = d$recipient_line[1], donor_line = d$donor_line[1], n_markers = d$n_markers[1], n_same = length(a), n_other = length(b), cpm_same = median(a), cpm_other = median(b), ratio = median(a) / max(median(b), 0.01), p = if (min(length(a), length(b)) >= 3) wilcox.test(a, b)$p.value else NA, rate_pct = 100 * median(d$rate[d$same_group])) }))
print(T3[order(T3$donor_line, T3$recipient_line), ], digits = 3, row.names = FALSE)
# CCLE-defined "off" genes still have a background of 20-50 CPM in these cultures, so the marker rates above are not usable as a rate; the 10 classical liver genes are (background 7 CPM):
LIVG <- intersect(c("ALB", "APOA1", "APOA2", "APOB", "FGB", "FGG", "TF", "TTR", "SERPINA1", "AFP"), rownames(cs)); livc <- colSums(cs[LIVG, ]); lg <- names(which(tapply(SM$kind %in% c("liver", "mixed"), SM$i5, any)))
base <- median((1e6 * livc / depth)[SM$sample[SM$kind == "meso" & !SM$i5 %in% lg]]); LR <- do.call(rbind, lapply(lg, function(g) { ss_ <- SM$sample[SM$i5 == g]; pool <- 1e6 * sum(livc[ss_]) / sum(depth[ss_]); r <- SM$sample[SM$i5 == g & SM$kind == "meso"]
  data.frame(i5 = g, recipient = r, liver_cpm = 1e6 * livc[r] / depth[r], pool_cpm = pool, rate = (1e6 * livc[r] / depth[r] - base) / pool) })); write.csv(LR, file.path(OUT, "spec125_liver_gene_carryover.csv"), row.names = FALSE)
cat(sprintf("liver genes (10): background %.1f CPM in the groups without liver samples; carry-over rate (recipient CPM above background / pooled CPM of its i5 group):\n", base)); print(round(100 * tapply(LR$rate, LR$i5, median), 2)); rt <- LR

## ---------------------------------------------------------------- 4  guide carry-over
gr <- read.csv(file.path(OUT, "spec119_guide_reads.csv"), stringsAsFactors = FALSE); gr$gene <- sub("_g.*", "", gr$guide); gc <- read.csv(file.path(OUT, "spec119_guide_calls.csv"), stringsAsFactors = FALSE); carr <- setNames(ifelse(gc$guide_reads >= 300 & gc$top_share >= 0.8, gc$top_target, NA), gc$sample)
GA <- aggregate(reads ~ sample + gene, gr, sum); um <- tapply(gr$unmapped_reads, gr$sample, function(v) v[1]); GA$rpm <- 1e6 * GA$reads / um[GA$sample]; GA <- GA[GA$sample %in% rec$sample & !GA$gene %in% c("NTC", "NR4A2enh"), ]; GA$own <- carr[GA$sample]; GA <- GA[is.na(GA$own) | GA$gene != GA$own, ]
GA$i5 <- SM[GA$sample, "i5"]; GA$same_group <- mapply(function(s, k) any(SM$i5 == SM[s, "i5"] & SM$sample != s & !is.na(carr[SM$sample]) & carr[SM$sample] == k, na.rm = TRUE), GA$sample, GA$gene)
GA$pool <- mapply(function(s, k) { d <- SM$sample[which(SM$i5 == SM[s, "i5"] & SM$sample != s & !is.na(carr[SM$sample]) & carr[SM$sample] == k)]; if (!length(d)) NA else 1e6 * sum(gr$reads[gr$sample %in% d & gr$gene == k]) / sum(um[intersect(SM$sample[SM$i5 == SM[s, "i5"]], names(um))]) }, GA$sample, GA$gene); GA$rate <- GA$rpm / GA$pool
write.csv(GA, file.path(OUT, "spec125_guide_carryover.csv"), row.names = FALSE); w4 <- wilcox.test(GA$rpm[GA$same_group], GA$rpm[!GA$same_group])
cat(sprintf("\n=== 4  guide reads of another KO gene (per million unmapped): donor KO in the same i5 group: median %.2f (n = %d) | not in the group: median %.2f (n = %d) | Wilcoxon p %s | carry-over rate median %.2f%%\n", median(GA$rpm[GA$same_group]), sum(GA$same_group), median(GA$rpm[!GA$same_group]), sum(!GA$same_group), fp(w4$p.value), 100 * median(GA$rate, na.rm = TRUE)))

## ---------------------------------------------------------------- 5  private SNVs
S <- read.csv(file.path(OUT, "spec125_private_SNVs.csv"), stringsAsFactors = FALSE); P <- read.csv(file.path(OUT, "spec125_private_SNV_counts.csv"), stringsAsFactors = FALSE); P$line <- S$line[P$site]; P <- P[P$depth >= 8 & P$sample %in% SM$sample, ]
VF <- do.call(rbind, lapply(split(P, list(P$sample, P$line), drop = TRUE), function(d) data.frame(sample = d$sample[1], snv_line = d$line[1], n_sites = nrow(d), alt = sum(d$alt), depth = sum(d$ref + d$alt), frac = sum(d$alt) / max(1, sum(d$ref + d$alt)), n_sites_alt = sum(d$alt >= 2 & d$alt / pmax(1, d$ref + d$alt) >= 0.1))))
VF$line <- SM[VF$sample, "line"]; VF$extra <- SM[VF$sample, "extra"]; VF$kind <- SM[VF$sample, "kind"]; VF$i5 <- SM[VF$sample, "i5"]; write.csv(VF, file.path(OUT, "spec125_SNV_fraction.csv"), row.names = FALSE)
VM <- tapply(VF$frac, list(VF$sample, VF$snv_line), mean)[, names(CC)]; best <- colnames(VM)[apply(VM, 1, which.max)]; names(best) <- rownames(VM); x <- SM[rownames(VM), ]
cat("\n=== 5  private SNVs: line with the highest alternative-allele fraction vs the sample's line ===\n"); print(table(line = x$line[x$kind %in% c("meso", "liver")], snv_best = best[x$kind %in% c("meso", "liver")]))
cat("mixed samples (alt fraction at the sites of H28 and of the second line):\n"); for (s in names(MIX)) if (s %in% rownames(VM)) cat(sprintf("  %s: NCI-H28 %.2f, %s %.2f\n", s, VM[s, "NCI-H28"], MIX[[s]], VM[s, MIX[[s]]]))
V5 <- VF[VF$kind == "meso" & VF$snv_line != VF$line, ]; V5$same_group <- mapply(function(s, M) M %in% present(SM[s, "i5"]), V5$sample, V5$snv_line)
t5 <- aggregate(cbind(alt, depth) ~ same_group + I(snv_line %in% LIVER), V5, sum); colnames(t5)[2] <- "donor_is_liver"; t5$alt_frac_pct <- 100 * t5$alt / t5$depth; cat("alternative-allele reads at another line's private sites, pooled over the mesothelioma samples:\n"); print(t5, digits = 3, row.names = FALSE)
a <- V5[V5$same_group, ]; b <- V5[!V5$same_group, ]; ft <- fisher.test(matrix(c(sum(a$alt), sum(a$depth) - sum(a$alt), sum(b$alt), sum(b$depth) - sum(b$alt)), 2)); cat(sprintf("same group %.3f%% vs other %.3f%% of reads carry the other line's allele: odds ratio %.1f, Fisher p %s\n", 100 * sum(a$alt) / sum(a$depth), 100 * sum(b$alt) / sum(b$depth), ft$estimate, fp(ft$p.value)))

## ---------------------------------------------------------------- figure
I5 <- sort(unique(SM$i5[SM$kind == "meso"])); I5C <- setNames(c("#2a78d6", "#1baf7a", "#eb6834", "#eda100", "#4a3aa7"), I5); SHL <- function(l) sub("NCI-|-211H", "", l)
save_fig("fig42_i5_group_carryover", 17, 10.4, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 4.2, 0.4), family = "sans", mar = c(5.4, 4.8, 3.6, 1)); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.6, font = 2, cex = 0.8); mtext(b, side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2) }
  ## A  excess profile vs CCLE
  col <- ifelse(EX$lineage == "LIVER", "#c8322b", "grey75"); plot(EX$rank, EX$rho, pch = 16, cex = ifelse(EX$lineage == "LIVER", 0.9, 0.4), col = col, bty = "n", las = 1, xlab = "CCLE cell lines, ranked", ylab = "Spearman with the excess profile of the new H2052 controls", cex.axis = 0.8); abline(h = 0, lty = 2, col = GRID)
  w <- match(CC[LIVER], EX$ccle); points(EX$rank[w], EX$rho[w], pch = 21, bg = "#c8322b", col = "black", cex = 1.3); text(EX$rank[w], EX$rho[w], LIVER, pos = 4, cex = 0.55, xpd = NA); legend("topright", c("liver lines in CCLE", "liver lines of this run", "other lineages"), pch = c(16, 21, 16), col = c("#c8322b", "black", "grey75"), pt.bg = "#c8322b", bty = "n", cex = 0.65)
  ttl("A  is the excess of the new H2052 controls liver-like overall?", sprintf("H28_NTC_R1-R3 minus H2052 KO samples vs %d CCLE lines; liver lineage median rank %d; top lineage of the first 30: %s", nrow(EX), as.integer(median(EX$rank[EX$lineage == "LIVER"])), tolower(names(sort(table(EX$lineage[1:30]), decreasing = TRUE))[1])))
  ## B  liver marker CPM by i5 group
  ms <- SM$sample[SM$kind == "meso"]; a <- data.frame(recipient = ms, i5 = factor(SM[ms, "i5"], levels = I5), marker_cpm = 1e6 * livc[ms] / depth[ms]); boxplot(log10(marker_cpm + 0.1) ~ i5, a, col = adjustcolor(I5C, 0.5), border = INK2, outline = FALSE, las = 2, xlab = "", ylab = "10 liver genes (ALB, APOA1, ...), log10 CPM", cex.axis = 0.62, ylim = range(log10(a$marker_cpm + 0.1)))
  points(jitter(as.integer(a$i5), 0.7), log10(a$marker_cpm + 0.1), pch = 16, cex = 0.5, col = "grey35"); nl <- sapply(I5, function(g) sum(SM$i5 == g & SM$kind %in% c("liver", "mixed"))); mtext(sprintf("%d liver", nl), side = 3, at = seq_along(I5), line = -1, cex = 0.58, col = INK2)
  ttl("B  liver transcripts in the mesothelioma samples, by i5 group", "'n liver' = liver-line or mixed samples that share the i5 index")
  ## C  marker carry-over between mesothelioma lines
  d <- MK[MK$donor_kind == "meso", ]; d$grp <- ifelse(d$same_group, "donor line in the same i5 group", "donor line not in the group"); boxplot(log10(marker_cpm + 0.1) ~ grp, d, col = c("#f0a58f", "#d5d4cf"), border = INK2, outline = FALSE, las = 1, xlab = "", ylab = "markers of another mesothelioma line, log10 summed CPM", cex.axis = 0.7)
  w3 <- wilcox.test(d$marker_cpm[d$same_group], d$marker_cpm[!d$same_group]); ttl("C  markers of another mesothelioma line (CCLE-defined)", sprintf("recipient x donor-line pairs: same group n = %d (median %.1f CPM), other n = %d (median %.1f); Wilcoxon p %s", sum(d$same_group), median(d$marker_cpm[d$same_group]), sum(!d$same_group), median(d$marker_cpm[!d$same_group]), fp(w3$p.value)))
  ## D  guides
  GA$grp <- ifelse(GA$same_group, "KO of that gene in the same i5 group", "not in the group"); boxplot(log10(rpm + 0.01) ~ grp, GA, col = c("#f0a58f", "#d5d4cf"), border = INK2, outline = FALSE, las = 1, xlab = "", ylab = "guide reads of another KO gene, log10 per million unmapped", cex.axis = 0.7)
  ttl("D  guide reads follow the same rule", sprintf("sample x foreign KO gene: same group median %.2f, other %.2f per million; Wilcoxon p %s", median(GA$rpm[GA$same_group]), median(GA$rpm[!GA$same_group]), fp(w4$p.value)))
  ## E  carry-over rate
  r1 <- 100 * pmax(rt$rate, 0); r2 <- 100 * GA$rate[!is.na(GA$rate) & is.finite(GA$rate)]; boxplot(list(`liver genes` = log10(r1 + 1e-3), `guide reads` = log10(r2 + 1e-3)), col = c("#8fb8ea", "#f0a58f"), border = INK2, outline = FALSE, las = 1, ylab = "carry-over rate, log10 %", cex.axis = 0.8)
  abline(h = log10(c(0.1, 1)), lty = 3, col = MUTED); text(2.45, log10(c(0.1, 1)), c("0.1%", "1%"), cex = 0.6, pos = 3, xpd = NA); ttl("E  size of the carry-over", sprintf("recipient signal / pooled signal of its i5 group: liver genes median %.2f%%, guides median %.2f%%", median(r1), median(r2)))
  ## F  SNV identity
  o <- order(match(x$kind, c("meso", "liver", "mixed", "other")), match(x$line, names(CC)), rownames(VM)); Mx <- VM[o, ]; par(mar = c(5.4, 4.8, 3.6, 4)); pal <- colorRampPalette(c("#f7f6f2", "#f0a58f", "#a8380b"))(100)
  image(seq_len(nrow(Mx)), seq_len(ncol(Mx)), pmin(Mx, 0.6), col = pal, zlim = c(0, 0.6), axes = FALSE, xlab = "samples (mesothelioma lines, liver lines, mixed)", ylab = ""); axis(2, at = seq_len(ncol(Mx)), labels = colnames(Mx), las = 1, tick = FALSE, cex.axis = 0.62); box()
  bnd <- which(diff(as.integer(factor(paste(x$kind[o], x$line[o]), levels = unique(paste(x$kind[o], x$line[o]))))) != 0) + 0.5; abline(v = bnd, col = "grey40", lwd = 0.6)
  ttl("F  private SNVs of each line (CCLE) confirm the identities", "fraction of reads with the alternative allele at the line's private sites (depth >= 8); colour capped at 0.6")
  mtext("Trace liver transcripts in some i5 index groups: carry-over between samples that share an i5 index", outer = TRUE, side = 3, adj = 0, line = 2.6, font = 2, cex = 1)
  mtext("i5 (index2) defines 5 groups of samples; within a group the samples differ only by i7. CCLE: expression of 1019 lines and private somatic SNVs of the 11 lines of this run", outer = TRUE, side = 3, adj = 0, line = 1.4, cex = 0.65, col = INK2)
})
