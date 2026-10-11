# Spec 119: which guides does each bulk sample carry? Reads that contain a gRNA protospacer followed by the start of the sgRNA scaffold (20-mer + GTTT, either orientation) are counted among the read pairs STAR left unmapped
# (kept in the BAM with --outSAMunmapped Within and stored at the end of the file, so `samtools view <bam> '*'` reads only those): transcripts of the guide vector contain the protospacer
# followed by the sgRNA scaffold and do not align to the genome. The protospacer alone also occurs in transcripts of the target gene (some of which stay unmapped), hence the scaffold anchor.
# needs: samtools on the PATH (module load samtools), BAMs in BAM_DIR, "Oligos NR4A2+CROPseq.xlsx"
# outputs: spec119_guide_reads.csv (sample x guide read counts, unmapped pairs scanned), spec119_guide_calls.csv (guides detected per sample, expected target, label check)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R")); options(width = 220)
BAM_DIR <- Sys.getenv("BAM_DIR", unset = "/sc/arion/scratch/giottb01/CROPseq_validation_bulkRNA/work/bam"); NCORE <- as.integer(Sys.getenv("NCORE", unset = "8")); MIN_READS <- 3
g <- as.data.frame(readxl::read_excel(file.path(here, "Oligos NR4A2+CROPseq.xlsx"), sheet = "gRNA")); colnames(g) <- c("target", "n", "seq"); g <- g[!is.na(g$seq), ]; g$target <- trimws(g$target); g$seq <- toupper(trimws(g$seq))
g$target[g$target == "NR4A2 Enhancer"] <- "NR4A2enh"; g$target[g$target == "Non Targeting Guide"] <- "NTC"; g$guide <- sprintf("%s_g%d", g$target, as.integer(g$n)); rc <- function(x) chartr("ACGT", "TGCA", sapply(lapply(strsplit(x, ""), rev), paste, collapse = ""))
SCAF <- "GTTT"   # first bases of the sgRNA scaffold (GTTTTAGAGCTA / GTTTAAGAGCTA): requiring them right after the protospacer excludes reads of the target gene itself, where an NGG PAM follows
pat <- c(setNames(g$guide, paste0(g$seq, SCAF)), setNames(g$guide, rc(paste0(g$seq, SCAF)))); pf <- tempfile(fileext = ".txt"); writeLines(names(pat), pf)
ids <- list.dirs(BAM_DIR, full.names = FALSE, recursive = FALSE); bam <- setNames(file.path(BAM_DIR, ids, paste0(ids, "_Aligned.sortedByCoord.out.bam")), ids); bam <- bam[file.exists(bam) & file.exists(paste0(bam, ".bai"))]
f_out <- file.path(OUT, "spec119_guide_reads.csv")
if (!file.exists(f_out) || Sys.getenv("RECOUNT", "0") == "1") {
  one <- function(s) { x <- system2("bash", c("-c", shQuote(sprintf("samtools view %s '*' | cut -f10 | grep -oFf %s | sort | uniq -c", shQuote(bam[[s]]), pf))), stdout = TRUE); n <- as.numeric(system2("samtools", c("view", "-c", shQuote(bam[[s]]), "'*'"), stdout = TRUE))
    v <- setNames(numeric(nrow(g)), g$guide); if (length(x)) { cn <- as.numeric(sub("^\\s*([0-9]+)\\s+.*", "\\1", x)); sq <- sub("^\\s*[0-9]+\\s+", "", x); a <- tapply(cn, pat[sq], sum); v[names(a)] <- a }
    data.frame(sample = s, unmapped_reads = n, guide = g$guide, reads = as.numeric(v)) }
  R <- do.call(rbind, parallel::mclapply(names(bam), one, mc.cores = NCORE)); write.csv(R, f_out, row.names = FALSE) }
R <- read.csv(f_out, stringsAsFactors = FALSE); R$target <- g$target[match(R$guide, g$guide)]; R$rpm <- 1e6 * R$reads / pmax(1, R$unmapped_reads)
q <- load_audit(); pm <- regmatches(R$sample, regexec("^(H2052|MSTO|H2452|H28)_([A-Za-z0-9]+)_R([0-9]+)$", R$sample)); R$KO <- sapply(pm, function(m) if (length(m) == 4) m[3] else NA); i <- match(R$sample, q$sample); R$KO[!is.na(i)] <- q$KO[i[!is.na(i)]]
R$KO[grepl("^NR4A2_", R$sample)] <- "NR4A2enh"; R$KO[grepl("^SCRAMBLE_", R$sample)] <- "NTC"
CALL <- do.call(rbind, lapply(split(R, R$sample), function(d) { h <- d[d$reads >= MIN_READS, ]; h <- h[order(-h$reads), ]; tt <- tapply(d$reads, d$target, sum); top <- names(which.max(tt))
  data.frame(sample = d$sample[1], label = d$KO[1], unmapped_M = d$unmapped_reads[1] / 1e6, guide_reads = sum(d$reads), top_target = if (sum(tt) >= MIN_READS) top else NA, top_share = if (sum(tt) > 0) max(tt) / sum(tt) else NA,
             own_reads = if (!is.na(d$KO[1]) && d$KO[1] %in% names(tt)) tt[[d$KO[1]]] else NA, guides = paste(sprintf("%s:%d", h$guide, h$reads), collapse = " ")) }))
CALL$match <- ifelse(is.na(CALL$label) | is.na(CALL$top_target), NA, CALL$label == CALL$top_target); write.csv(CALL, file.path(OUT, "spec119_guide_calls.csv"), row.names = FALSE)
print(CALL[, c("sample", "label", "unmapped_M", "guide_reads", "top_target", "top_share", "own_reads", "match", "guides")], digits = 2, row.names = FALSE)
cat("\nlabel vs guide detected (samples with >=", MIN_READS, "guide reads):\n"); print(table(label = CALL$label, detected = CALL$top_target, useNA = "ifany"))
