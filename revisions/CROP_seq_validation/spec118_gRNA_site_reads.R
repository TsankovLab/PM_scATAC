# Spec 118: read distribution at the gRNA target sites in the bulk RNA-seq BAMs -- is there evidence of editing at the site of each guide in the samples that received it?
# guides: "Oligos NR4A2+CROPseq.xlsx" (sheet gRNA; 3 guides per target, 3 non-targeting, 3 for the NR4A2 enhancer), mapped to hg38 by exact match on both strands (BSgenome), PAM (NGG) and
#         Cas9 cut position (3 bp 5' of the PAM) recorded; hits inside the target gene (GENCODE v46 gene body) are kept, otherwise all hits
# per sample x site, from the STAR alignments (primary alignments, window = cut +/- W bp, W = 10):
#   n_cover  = read pairs with aligned bases (or a deletion) in the window
#   n_indel  = pairs with an insertion or deletion in the window         n_clip = pairs soft-clipped with the clip boundary in the window
#   n_junc   = pairs with a splice junction (N) boundary in the window (large deletions align as N; annotated junctions next to the site also count -> read against the controls)
#   indel = n_indel / n_cover is the primary editing metric (soft clips are frequent away from any guide: adapter read-through, exon boundaries); edit = (indel or clip) / n_cover; edit_all also counts junction boundaries
#   per-base profiles (cut +/- FLANK bp): depth, deleted bases, read ends -- the "read distribution" around the cut
# comparison: samples whose KO is the guide's target vs all other samples of the same (genotyped) cell line; per guide Fisher test on pooled pairs
# NR4A2 enhancer guides (CRISPRa, no cut expected): the same counts plus coverage (reads per million uniquely mapped) in +/- 1 kb, NR4A2_1-3 vs SCRAMBLE_1-3
# needs: BAMs in BAM_DIR (STAR, coordinate sorted, indexed), spec71 audit (genotyped line and KO of each sample)
# outputs: spec118_gRNA_sites.csv, spec118_site_counts.csv, spec118_guide_summary.csv, spec118_profiles.rds, fig33_gRNA_site_editing.{png,pdf}, fig34_gRNA_site_read_profiles.{png,pdf}, fig35_NR4A2_enhancer_guides.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages({ library(GenomicAlignments); library(Rsamtools); library(Biostrings); library(BSgenome.Hsapiens.UCSC.hg38) }); options(width = 220)
XLSX <- file.path(here, "Oligos NR4A2+CROPseq.xlsx"); BAM_DIR <- Sys.getenv("BAM_DIR", unset = "/sc/arion/scratch/giottb01/CROPseq_validation_bulkRNA/work/bam")
GTF <- "/sc/arion/projects/Tsankov_Normal_Lung/genome_assemblies/human/gencode.v46.annotation.gtf"; W <- 10L; FLANK <- 60L; ENH_FLANK <- 1000L
fp <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))); SH <- function(l) sub("NCI-|-211H", "", l)

## ---------------------------------------------------------------- guides -> genome
g <- as.data.frame(readxl::read_excel(XLSX, sheet = "gRNA")); colnames(g) <- c("target", "n", "seq"); g <- g[!is.na(g$seq), ]; g$target <- trimws(g$target); g$seq <- toupper(trimws(g$seq))
g$target[g$target == "NR4A2 Enhancer"] <- "NR4A2enh"; g$target[g$target == "Non Targeting Guide"] <- "NTC"; g$guide <- sprintf("%s_g%d", g$target, as.integer(g$n))
f_sites <- file.path(OUT, "spec118_gRNA_sites.csv")
if (!file.exists(f_sites) || Sys.getenv("REMAP", "0") == "1") {
  fw <- DNAStringSet(setNames(g$seq, g$guide)); rc <- reverseComplement(fw); chrs <- paste0("chr", c(1:22, "X", "Y")); H <- list()
  for (ch in chrs) { s <- BSgenome.Hsapiens.UCSC.hg38[[ch]]
    for (sd in c("+", "-")) { m <- matchPDict(PDict(if (sd == "+") fw else rc), s); n <- elementNROWS(m); if (!sum(n)) next
      H[[length(H) + 1]] <- data.frame(guide = rep(names(fw), n), chr = ch, start = unlist(start(m)), end = unlist(end(m)), strand = sd) } }
  S <- do.call(rbind, H); gs <- BSgenome.Hsapiens.UCSC.hg38
  S$PAM <- mapply(function(ch, s, e, sd) if (sd == "+") as.character(getSeq(gs, ch, e + 1, e + 3)) else as.character(reverseComplement(getSeq(gs, ch, s - 3, s - 1))), S$chr, S$start, S$end, S$strand)
  S$cut <- ifelse(S$strand == "+", S$end - 3L, S$start + 2L)                                    # Cas9 cuts between cut and cut + 1 (3 bp 5' of the PAM)
  S$target <- g$target[match(S$guide, g$guide)]; tg <- setdiff(unique(g$target), c("NTC", "NR4A2enh"))
  gl <- read.delim(pipe(sprintf("awk -F'\\t' '$3==\"gene\"' %s | grep -E 'gene_name \"(%s)\";' | cut -f1,4,5,7,9", GTF, paste(tg, collapse = "|"))), header = FALSE, quote = "", col.names = c("chr", "gs", "ge", "gstrand", "attr"))
  gl$gene <- sub('.*gene_name "([^"]+)".*', "\\1", gl$attr); i <- match(S$target, gl$gene); S$in_target_gene <- !is.na(i) & S$chr == gl$chr[i] & S$start >= gl$gs[i] & S$end <= gl$ge[i]; S$gene_strand <- gl$gstrand[i]
  S$n_hits <- as.integer(table(S$guide)[S$guide]); S <- S[S$in_target_gene | !S$guide %in% S$guide[S$in_target_gene], ]; S <- merge(g[, c("guide", "seq")], S, by = "guide", all.x = TRUE); write.csv(S, f_sites, row.names = FALSE) }
S <- read.csv(f_sites, stringsAsFactors = FALSE); cat("guides:", nrow(g), "| mapped exactly:", sum(!is.na(S$chr[!duplicated(S$guide)])), "| not found:", paste(S$guide[is.na(S$chr)], collapse = " "), "\n")
print(S[!is.na(S$chr), c("guide", "chr", "start", "end", "strand", "PAM", "cut", "in_target_gene", "n_hits")], row.names = FALSE); S <- S[!is.na(S$chr), ]; S$site <- sprintf("%s:%d", S$chr, S$cut)

## ---------------------------------------------------------------- samples
ids <- list.dirs(BAM_DIR, full.names = FALSE, recursive = FALSE); bam <- file.path(BAM_DIR, ids, paste0(ids, "_Aligned.sortedByCoord.out.bam")); ok <- file.exists(bam) & file.exists(paste0(bam, ".bai")); ids <- ids[ok]; bam <- setNames(bam[ok], ids)
q <- load_audit(); PFX <- c(H2052 = "NCI-H2052", MSTO = "MSTO-211H", H2452 = "NCI-H2452", H28 = "NCI-H28"); pm <- regmatches(ids, regexec("^(H2052|MSTO|H2452|H28)_([A-Za-z0-9]+)_R([0-9]+)$", ids))
SM <- data.frame(sample = ids, line = sapply(pm, function(m) if (length(m) == 4) PFX[[m[2]]] else NA), KO = sapply(pm, function(m) if (length(m) == 4) m[3] else NA), stringsAsFactors = FALSE)
i <- match(SM$sample, q$sample); SM$line[!is.na(i)] <- q$line[i[!is.na(i)]]; SM$KO[!is.na(i)] <- q$KO[i[!is.na(i)]]                         # genotyped line / KO where audited (NTC reassignment included)
SM$line[grepl("^(NR4A2|SCRAMBLE)_", SM$sample)] <- "NR4A2 CRISPRa"; SM$KO[grepl("^NR4A2_", SM$sample)] <- "NR4A2enh"; SM$KO[grepl("^SCRAMBLE_", SM$sample)] <- "NTC"; SM <- SM[!is.na(SM$line) & !grepl("^MIXED", SM$sample), ]
SM$uniq_M <- sapply(SM$sample, function(s) { l <- readLines(file.path(BAM_DIR, s, paste0(s, "_Log.final.out"))); as.numeric(sub(".*\\|\\s*", "", grep("Uniquely mapped reads number", l, value = TRUE))) / 1e6 })
cat("\nsamples with BAM:", nrow(SM), "\n"); print(table(SM$line, SM$KO))

## ---------------------------------------------------------------- counts and profiles per sample x site
f_cnt <- file.path(OUT, "spec118_site_counts.csv"); f_prof <- file.path(OUT, "spec118_profiles.rds")
if (!file.exists(f_cnt) || !file.exists(f_prof) || Sys.getenv("RECOUNT", "0") == "1") {
  S$fl <- ifelse(S$target == "NR4A2enh", ENH_FLANK, FLANK); gr <- GRanges(S$chr, IRanges(S$cut - S$fl, S$cut + S$fl)); lab <- sprintf("%s:%d-%d", S$chr, S$cut - S$fl, S$cut + S$fl)
  prm <- ScanBamParam(which = gr, what = "qname", flag = scanBamFlag(isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE)); CNT <- list(); PROF <- list()
  in_win <- function(rl, lo, hi, ends = FALSE) { u <- unlist(rl, use.names = FALSE); if (!length(u)) return(integer(0)); gp <- togroup(PartitioningByWidth(elementNROWS(rl)))
    h <- if (ends) (start(u) >= lo[gp] & start(u) <= hi[gp]) | (end(u) >= lo[gp] & end(u) <= hi[gp]) else start(u) <= hi[gp] & end(u) >= lo[gp]; unique(gp[h]) }
  for (s in SM$sample) { ga <- readGAlignments(bam[[s]], param = prm, with.which_label = TRUE); wl <- as.character(mcols(ga)$which_label)
    for (j in seq_len(nrow(S))) { a <- ga[wl == lab[j]]; cut <- S$cut[j]; lo <- rep(cut - W + 1L, length(a)); hi <- rep(cut + W, length(a)); z <- data.frame(sample = s, guide = S$guide[j], site = S$site[j], n_cover = 0L, n_indel = 0L, n_clip = 0L, n_junc = 0L, n_edit = 0L, n_edit_all = 0L, n_region = length(unique(mcols(a)$qname)))
      if (length(a)) { cg <- cigar(a); st <- start(a); qn <- mcols(a)$qname; M <- cigarRangesAlongReferenceSpace(cg, pos = st, ops = c("M", "=", "X")); Dl <- cigarRangesAlongReferenceSpace(cg, pos = st, ops = "D"); In <- cigarRangesAlongReferenceSpace(cg, pos = st, ops = "I"); N <- cigarRangesAlongReferenceSpace(cg, pos = st, ops = "N")
        cov_i <- union(in_win(M, lo, hi), in_win(Dl, lo, hi)); ind_i <- union(in_win(Dl, lo, hi), { u <- unlist(In, use.names = FALSE); gp <- togroup(PartitioningByWidth(elementNROWS(In))); unique(gp[start(u) >= lo[gp] & start(u) <= hi[gp] + 1L]) })
        clipL <- grepl("^[0-9]+S", cg) & st >= lo & st <= hi + 1L; clipR <- grepl("S$", cg) & end(a) >= lo - 1L & end(a) <= hi; clp_i <- which(clipL | clipR)
        Nb <- N; start(Nb) <- start(N) - 1L; end(Nb) <- end(N) + 1L; jn_i <- in_win(Nb, lo, hi, ends = TRUE)                                    # last / first aligned base flanking the junction
        z$n_cover <- length(unique(qn[union(cov_i, clp_i)])); z$n_indel <- length(unique(qn[ind_i])); z$n_clip <- length(unique(qn[clp_i])); z$n_junc <- length(unique(qn[jn_i])); z$n_edit <- length(unique(qn[union(ind_i, clp_i)])); z$n_edit_all <- length(unique(qn[Reduce(union, list(ind_i, clp_i, jn_i))]))
        fl <- S$fl[j]; x <- (cut - fl):(cut + fl); tab <- function(v) tabulate(v - (cut - fl) + 1L, nbins = length(x)); cv <- function(rl) { u <- unlist(rl, use.names = FALSE); if (!length(u)) return(integer(length(x))); as.integer(coverage(u, shift = -(cut - fl) + 1L, width = length(x))) }
        PROF[[paste(s, S$guide[j])]] <- rbind(depth = cv(M), deleted = cv(Dl), spliced = cv(N), ends = tab(c(st[st >= min(x) & st <= max(x)], end(a)[end(a) >= min(x) & end(a) <= max(x)]))) }
      CNT[[length(CNT) + 1]] <- z }
    cat(".") }
  cat("\n"); C <- do.call(rbind, CNT); write.csv(C, f_cnt, row.names = FALSE); saveRDS(PROF, f_prof) }
C <- read.csv(f_cnt, stringsAsFactors = FALSE); PROF <- readRDS(f_prof); C <- C[paste(C$guide, C$site) %in% paste(S$guide, S$site), ]; C <- merge(C, SM, by = "sample"); C$target <- S$target[match(C$guide, S$guide)]; C$own <- C$KO == C$target
C$indel <- ifelse(C$n_cover > 0, C$n_indel / C$n_cover, NA); C$edit <- ifelse(C$n_cover > 0, C$n_edit / C$n_cover, NA); C$edit_all <- ifelse(C$n_cover > 0, C$n_edit_all / C$n_cover, NA)   # indel = primary metric (soft clips are common away from any guide)

## ---------------------------------------------------------------- per guide x line: own KO vs the other samples of the line
BULK <- C[C$line %in% LINES & C$target != "NTC" & C$target != "NR4A2enh", ]
G <- do.call(rbind, lapply(split(BULK, list(BULK$guide, BULK$line), drop = TRUE), function(d) { o <- d[d$own, ]; b <- d[!d$own, ]; if (!nrow(o)) return(NULL)
  ft <- function(k) if (sum(o$n_cover) > 0 && sum(b$n_cover) > 0) fisher.test(matrix(c(sum(o[[k]]), sum(o$n_cover) - sum(o[[k]]), sum(b[[k]]), sum(b$n_cover) - sum(b[[k]])), 2), alternative = "greater")$p.value else NA
  pc <- function(x, k) 100 * sum(x[[k]]) / max(1, sum(x$n_cover))
  data.frame(target = d$target[1], guide = d$guide[1], line = d$line[1], n_KO_samples = nrow(o), KO_pairs = sum(o$n_cover), KO_indel_pct = pc(o, "n_indel"), other_pairs = sum(b$n_cover), other_indel_pct = pc(b, "n_indel"), p = ft("n_indel"),
             KO_indel_per_sample = paste(ifelse(o$n_cover > 0, sprintf("%.0f", 100 * o$indel), "-"), collapse = "/"), KO_clip_pct = pc(o, "n_clip"), other_clip_pct = pc(b, "n_clip"), KO_edit_pct = pc(o, "n_edit"), other_edit_pct = pc(b, "n_edit"), p_edit = ft("n_edit"),
             KO_junc_pct = pc(o, "n_junc"), other_junc_pct = pc(b, "n_junc")) }))
G <- G[order(match(G$target, c(KOS, "RUNX1", "RUNX2", "SNAI2")), G$guide, match(G$line, LINES)), ]; write.csv(G, file.path(OUT, "spec118_guide_summary.csv"), row.names = FALSE)
cat("\n=== indels at the gRNA sites: own-KO samples vs the other samples of the same line (% of read pairs covering the cut +/-", W, "bp with an insertion or deletion there) ===\n")
print(G[, c("target", "guide", "line", "n_KO_samples", "KO_pairs", "KO_indel_pct", "other_pairs", "other_indel_pct", "p", "KO_indel_per_sample", "KO_clip_pct", "other_clip_pct")], digits = 2, row.names = FALSE)
MINP <- 10; Gk <- G[G$KO_pairs >= MINP, ]; best <- do.call(rbind, lapply(split(Gk, list(Gk$target, Gk$line), drop = TRUE), function(d) d[which.max(d$KO_indel_pct), ])); TGS <- intersect(c(KOS, "RUNX1", "RUNX2", "SNAI2"), G$target)
mat <- function(v) { M <- matrix(NA, length(TGS), 4, dimnames = list(TGS, LINES)); M[cbind(match(best$target, TGS), match(best$line, LINES))] <- v; M }; BM <- mat(best$KO_indel_pct)
cat(sprintf("\n-- best guide per KO x line (>= %d pairs at the site): %% of read pairs with an indel in the KO samples\n", MINP)); print(round(BM, 1)); cat("-- guide giving it\n"); print(noquote(mat(sub(".*_", "", best$guide))))
SUMG <- tapply(Gk$KO_indel_pct, list(Gk$target, Gk$line), sum)[TGS, LINES]; cat("-- sum over the guides of the target (upper bound if the guides were pooled)\n"); print(round(SUMG, 1))
bg <- C[C$line %in% LINES & !C$own & C$target %in% TGS & C$n_cover >= 20, ]; cat(sprintf("\nbackground (samples without the guide, >= 20 pairs): indel median %.2f%%, 99th percentile %.2f%%; indel or clip median %.2f%%, 95th percentile %.2f%%\n", 100 * median(bg$indel), 100 * quantile(bg$indel, 0.99), 100 * median(bg$edit), 100 * quantile(bg$edit, 0.95)))
cat("guides with significant indels (BH FDR < 0.05) per target:\n"); G$fdr <- p.adjust(G$p, "BH"); print(tapply(G$fdr < 0.05 & G$KO_indel_pct > 2, list(G$guide, G$line), any)[unique(G$guide), LINES] * 1)
# does the editing level relate to the earlier per KO x line results (spec102: agreement across lines, transcript knock-down)?
f102 <- file.path(OUT, "spec102_agreement_vs_baseline.csv"); if (file.exists(f102)) { a <- read.csv(f102); a$indel_best <- BM[cbind(match(a$KO, rownames(BM)), match(a$line, LINES))]; a$indel_sum <- SUMG[cbind(match(a$KO, rownames(SUMG)), match(a$line, LINES))]
  ct <- function(x, y) { r <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)); sprintf("rho %.2f, p %.2g (n=%d)", r$estimate, r$p.value, sum(complete.cases(x, y))) }
  cat("\nindel % (best guide) vs transcript knock-down (NTC - KO):", ct(a$indel_best, a$knockdown), "| vs agreement across lines:", ct(a$indel_best, a$agree_none), "| vs baseline expression:", ct(a$indel_best, a$baseline), "\n")
  cat("mean indel % (best guide) per line:", paste(SH(LINES), round(colMeans(BM[KOS, ], na.rm = TRUE), 1), collapse = ", "), "\n"); write.csv(a, file.path(OUT, "spec118_indel_vs_agreement.csv"), row.names = FALSE) }

## ---------------------------------------------------------------- figures
TG <- TGS; PCH <- setNames(c(21, 22, 24, 23), LINES)
save_fig("fig33_gRNA_site_editing", 16, 9.6, function() {
  layout(matrix(1:12, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.2, 0.4), family = "sans")
  for (t in TG) { d <- BULK[BULK$target == t & BULK$n_cover >= 5, ]; gs <- sort(unique(BULK$guide[BULK$target == t])); par(mar = c(3.6, 4.2, 3, 0.6))
    plot(NA, xlim = c(0.5, length(gs) + 0.5), ylim = c(0, max(100 * d$indel, 5, na.rm = TRUE) * 1.15), xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "% read pairs with an indel at the cut", cex.axis = 0.8)
    for (k in seq_along(gs)) { e <- d[d$guide == gs[k], ]; b <- e[!e$own, ]; o <- e[e$own, ]; if (nrow(b)) points(k - 0.2 + runif(nrow(b), -0.1, 0.1), 100 * b$indel, pch = 16, cex = 0.5, col = "grey70")
      if (nrow(o)) points(k + 0.15 + seq(-0.14, 0.14, length.out = nrow(o)), 100 * o$indel, pch = PCH[o$line], bg = LINE_COL[o$line], col = INK2, cex = 1.3)
      x <- S[S$guide == gs[k], ][1, ]; mtext(sprintf("g%s\n%s", sub(".*_g", "", gs[k]), ifelse(is.na(x$PAM), "", x$PAM)), side = 1, at = k, line = 1.6, cex = 0.6); mtext(sprintf("n=%.0f", median(e$n_cover)), side = 1, at = k, line = 2.5, cex = 0.5, col = MUTED) }
    mtext(t, side = 3, adj = 0, line = 1.2, font = 2, cex = 0.9, col = if (t %in% names(KO_COL)) KO_COL[[t]] else INK); mtext("colour = samples with this KO; grey = all other samples", side = 3, adj = 0, line = 0.2, cex = 0.5, col = INK2) }
  par(mar = c(4.6, 5.6, 3, 4.6)); zl <- max(BM, na.rm = TRUE); pal <- colorRampPalette(c("#f7f6f2", "#f3c4bf", "#c8322b"))(100); nr <- nrow(BM)
  image(1:4, 1:nr, t(BM[nr:1, , drop = FALSE]), col = pal, zlim = c(0, zl), axes = FALSE, xlab = "", ylab = ""); text(rep(1:4, each = nr), rep(nr:1, 4), ifelse(is.na(BM), "n/a", sprintf("%.0f", BM)), cex = 0.8, col = ifelse(is.na(BM), MUTED, INK))
  axis(1, at = 1:4, labels = SH(LINES), tick = FALSE, cex.axis = 0.8); axis(2, at = nr:1, labels = rownames(BM), tick = FALSE, las = 1, cex.axis = 0.8); mtext("best guide: % indel, KO samples", side = 3, adj = 0, line = 1, font = 2, cex = 0.72)
  mtext("Evidence of editing at the gRNA target sites in the bulk RNA-seq reads", outer = TRUE, side = 3, adj = 0, line = 1.7, font = 2, cex = 1)
  mtext(sprintf("read pairs covering the Cas9 cut site (+/- %d bp) that carry an insertion or deletion there; one point per sample (>= 5 pairs); n = median pairs at the site; symbols = cell line (H2052 circle, MSTO square, H2452 triangle, H28 diamond)", W), outer = TRUE, side = 3, adj = 0, line = 0.5, cex = 0.62, col = INK2)
})
mean_prof <- function(samples, guide, row, norm = TRUE) { P <- lapply(samples, function(s) PROF[[paste(s, guide)]]); P <- P[!sapply(P, is.null)]; if (!length(P)) return(NULL); M <- sapply(P, function(p) { v <- p[row, ]; if (!norm) return(v); d <- p["depth", ] + p["deleted", ]; if (row == "depth") v / max(1, max(d)) else v / pmax(1, d) }); rowMeans(M) }
save_fig("fig34_gRNA_site_read_profiles", 16, 8.4, function() {
  layout(matrix(seq_len(2 * length(TG)), nrow = 2)); par(oma = c(0, 0, 3.2, 0.4), family = "sans"); x <- -FLANK:FLANK
  for (t in TG) { gb <- names(which.max(tapply(Gk$KO_indel_pct[Gk$target == t], Gk$guide[Gk$target == t], mean))); d <- BULK[BULK$guide == gb, ]; o <- d$sample[d$own]; b <- d$sample[!d$own & d$KO == "NTC"]
    for (rw in c("depth", "deleted")) { yo <- mean_prof(o, gb, rw); yb <- mean_prof(b, gb, rw); par(mar = c(3.4, 3.8, 2.8, 0.5)); yl <- c(0, max(c(yo, yb, if (rw == "deleted") 0.05 else 1), na.rm = TRUE))
      plot(NA, xlim = range(x), ylim = yl, bty = "n", las = 1, xlab = "", ylab = if (rw == "depth") "relative depth" else "fraction of reads with the base deleted", cex.axis = 0.75, cex.lab = 0.8); abline(v = 0.5, lty = 2, col = MUTED)
      if (!is.null(yb)) lines(x, yb, col = "grey45", lwd = 1.6); if (!is.null(yo)) lines(x, yo, col = if (t %in% names(KO_COL)) KO_COL[[t]] else INK, lwd = 1.8); mtext("position relative to the cut (bp)", side = 1, line = 2.2, cex = 0.55)
      if (rw == "depth") { mtext(sprintf("%s (%s)", t, sub(".*_", "", gb)), side = 3, adj = 0, line = 1, font = 2, cex = 0.85, col = if (t %in% names(KO_COL)) KO_COL[[t]] else INK); mtext(sprintf("KO n=%d, NTC n=%d", length(o), length(b)), side = 3, adj = 0, line = 0.1, cex = 0.5, col = INK2) } } }
  mtext("Read distribution around the cut site of the best guide per target: samples with the KO (colour) vs NTC samples (grey), 4 cell lines pooled", outer = TRUE, side = 3, adj = 0, line = 1.7, font = 2, cex = 0.95)
  mtext("top: depth along the site, each sample scaled to its maximum; bottom: fraction of reads in which the base is deleted; dashed = Cas9 cut", outer = TRUE, side = 3, adj = 0, line = 0.5, cex = 0.62, col = INK2)
})
## NR4A2 enhancer guides
E <- C[C$target == "NR4A2enh" & C$line == "NR4A2 CRISPRa", ]; E$rpm_region <- E$n_region / E$uniq_M; cat("\n=== NR4A2 enhancer guides: NR4A2 CRISPRa vs SCRAMBLE samples ===\n")
print(E[order(E$guide, E$sample), c("guide", "site", "sample", "n_region", "rpm_region", "n_cover", "n_edit")], digits = 3, row.names = FALSE)
es <- do.call(rbind, lapply(split(E, E$guide), function(d) { a <- d$rpm_region[d$KO == "NR4A2enh"]; b <- d$rpm_region[d$KO == "NTC"]; data.frame(guide = d$guide[1], site = d$site[1], CRISPRa_rpm = mean(a), scramble_rpm = mean(b), p = if (sd(c(a, b)) > 0) t.test(a, b)$p.value else NA) })); print(es, digits = 3, row.names = FALSE)
eg <- S[S$target == "NR4A2enh", ]; if (nrow(eg)) cat(sprintf("enhancer guides span %s:%d-%d (%.1f kb)\n", eg$chr[1], min(eg$cut), max(eg$cut), diff(range(eg$cut)) / 1e3))
save_fig("fig35_NR4A2_enhancer_guides", 13, 4.6, function() { gs <- sort(unique(E$guide)); layout(matrix(seq_along(gs), nrow = 1)); par(oma = c(0, 0, 3, 0.4), family = "sans"); x <- -ENH_FLANK:ENH_FLANK; bin <- 25
  for (gd in gs) { sm <- function(v) as.numeric(stats::filter(v, rep(1 / bin, bin), sides = 2)); pr <- function(ko) { ss <- E$sample[E$guide == gd & E$KO == ko]; M <- sapply(ss, function(s) { p <- PROF[[paste(s, gd)]]; if (is.null(p)) rep(0, length(x)) else p["depth", ] / SM$uniq_M[SM$sample == s] }); M }
    A <- pr("NR4A2enh"); B <- pr("NTC"); yl <- c(0, max(c(sm(rowMeans(A)), sm(rowMeans(B)), 1e-3), na.rm = TRUE) * 1.2); par(mar = c(4.2, 4.4, 3, 0.8)); plot(NA, xlim = range(x), ylim = yl, bty = "n", las = 1, xlab = "position relative to the guide (bp)", ylab = "depth per million uniquely mapped reads", cex.axis = 0.8); abline(v = 0, lty = 2, col = MUTED)
    lines(x, sm(rowMeans(B)), col = "grey45", lwd = 1.8); lines(x, sm(rowMeans(A)), col = "#c8322b", lwd = 1.8); r <- es[es$guide == gd, ]; mtext(sprintf("%s  (%s)", gd, r$site), side = 3, adj = 0, line = 1.2, font = 2, cex = 0.8)
    mtext(sprintf("reads within 1 kb per million: CRISPRa %.2f vs scramble %.2f, p = %s", r$CRISPRa_rpm, r$scramble_rpm, fp(r$p)), side = 3, adj = 0, line = 0.2, cex = 0.58, col = INK2); if (gd == gs[1]) legend("topright", c("NR4A2 enhancer CRISPRa (n=3)", "scramble (n=3)"), col = c("#c8322b", "grey45"), lwd = 2, bty = "n", cex = 0.7) }
  mtext("RNA-seq reads at the NR4A2 enhancer guide sites: CRISPRa vs scramble", outer = TRUE, side = 3, adj = 0, line = 1.5, font = 2, cex = 1); mtext(sprintf("mean depth, %d-bp running mean; dashed = guide position", bin), outer = TRUE, side = 3, adj = 0, line = 0.4, cex = 0.62, col = INK2) })
