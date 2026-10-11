# Spec 125a: allele counts in every bulk BAM at the SNVs that are private to one of the 11 cell lines sequenced in this run, taken from the CCLE / DepMap somatic-mutation table.
# lines: the 4 mesothelioma lines (NCI-H2052, MSTO-211H, NCI-H2452, NCI-H28) and the 7 liver lines present as samples or mixtures (JHH-7, Li-7, SK-HEP-1, SNU-398, Hep G2, Hep 3B, HuH-7)
# sites: SNVs (hg38) reported in exactly one of the 11 models, in genes of the expressed-gene list (spec72); counts of the reference and alternative base per BAM (base quality >= 20, primary alignments)
# needs: CCLE OmicsSomaticMutations.csv, DepMap Model.csv, BAMs (STAR, indexed)
# outputs: spec125_private_SNVs.csv (sites), spec125_private_SNV_counts.csv (sample x site: ref, alt, depth)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R")); suppressPackageStartupMessages({ library(Rsamtools); library(GenomicRanges) })
CCLE <- Sys.getenv("CCLE_DIR", unset = "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Public_data/CCLE"); BAM_DIR <- Sys.getenv("BAM_DIR", unset = "/sc/arion/scratch/giottb01/CROPseq_validation_bulkRNA/work/bam"); NCORE <- as.integer(Sys.getenv("NCORE", unset = "10"))
MODELS <- c(`ACH-000153` = "NCI-H2052", `ACH-000335` = "MSTO-211H", `ACH-000092` = "NCI-H2452", `ACH-000648` = "NCI-H28", `ACH-000848` = "JHH-7", `ACH-000471` = "Li-7", `ACH-000361` = "SK-HEP-1", `ACH-000221` = "SNU-398", `ACH-000739` = "Hep G2", `ACH-000625` = "Hep 3B", `ACH-000480` = "HuH-7")
f_sites <- file.path(OUT, "spec125_private_SNVs.csv")
if (!file.exists(f_sites)) { cmd <- sprintf("awk -F, 'NR==1 || /%s/' %s | cut -d, -f1-9,11,15,68", paste(names(MODELS), collapse = "|"), shQuote(file.path(CCLE, "OmicsSomaticMutations.csv"))); mu <- read.csv(pipe(cmd), stringsAsFactors = FALSE)
  mu <- mu[mu$ModelID %in% names(MODELS) & mu$VariantType == "SNV" & nchar(mu$Ref) == 1 & nchar(mu$Alt) == 1, ]; mu$line <- MODELS[mu$ModelID]; mu$key <- paste(mu$Chrom, mu$Pos, mu$Alt); n <- table(mu$key); mu <- mu[n[mu$key] == 1, ]
  genes <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF)))$genes; mu <- mu[mu$HugoSymbol %in% genes & mu$Chrom %in% paste0("chr", c(1:22, "X")), ]; write.csv(mu[, c("Chrom", "Pos", "Ref", "Alt", "AF", "DP", "HugoSymbol", "line")], f_sites, row.names = FALSE) }
S <- read.csv(f_sites, stringsAsFactors = FALSE); cat("private SNVs in expressed genes per line:\n"); print(table(S$line)[MODELS])
ids <- list.dirs(BAM_DIR, full.names = FALSE, recursive = FALSE); bam <- setNames(file.path(BAM_DIR, ids, paste0(ids, "_Aligned.sortedByCoord.out.bam")), ids); bam <- bam[file.exists(bam) & file.exists(paste0(bam, ".bai"))]
gr <- GRanges(S$Chrom, IRanges(S$Pos, S$Pos)); pp <- PileupParam(max_depth = 1e6, min_base_quality = 20, min_mapq = 30, min_nucleotide_depth = 0, distinguish_strands = FALSE, include_deletions = FALSE, include_insertions = FALSE)
one <- function(s) { p <- pileup(bam[[s]], scanBamParam = ScanBamParam(which = gr, flag = scanBamFlag(isSecondaryAlignment = FALSE, isDuplicate = NA)), pileupParam = pp); k <- paste(p$seqnames, p$pos); sk <- paste(S$Chrom, S$Pos)
  cnt <- function(base) { v <- tapply(p$count[p$nucleotide == base], k[p$nucleotide == base], sum); v[is.na(v)] <- 0; v }; tot <- tapply(p$count, k, sum); g <- function(v, key) { x <- v[key]; x[is.na(x)] <- 0; as.numeric(x) }
  ref <- alt <- numeric(nrow(S)); for (b in c("A", "C", "G", "T")) { v <- cnt(b); ref[S$Ref == b] <- g(v, sk[S$Ref == b]); alt[S$Alt == b] <- g(v, sk[S$Alt == b]) }
  data.frame(sample = s, site = seq_len(nrow(S)), ref = ref, alt = alt, depth = g(tot, sk)) }
R <- do.call(rbind, parallel::mclapply(names(bam), one, mc.cores = NCORE)); R <- R[R$depth > 0, ]; write.csv(R, file.path(OUT, "spec125_private_SNV_counts.csv"), row.names = FALSE); cat("done:", length(unique(R$sample)), "samples,", nrow(R), "covered sample x site\n")
