###############################################################################
# R2_Q14 -- genetic BAP1 status for the 87 TCGA MESO tumours.
#
# Built from the public TCGA PanCancer Atlas calls fetched by
# fetch_TCGA_BAP1_cbioportal.sh (somatic mutations, GISTIC copy number, structural
# variants).  This replaces the BAP1 mRNA tertile proxy for TCGA with a label that
# is genetic and independent of the expression data it will be tested against.
#
# Classification, per tumour, in priority order:
#   ALTERED       truncating mutation (nonsense, frameshift, splice, nonstop), OR a
#                 BAP1 fusion/rearrangement, OR GISTIC deep deletion (-2).
#                 Each of these is expected to inactivate the allele it hits.
#   VUS           missense or in-frame indel only.  Kept separate: BAP1 missense can
#                 be damaging (C91 is the catalytic cysteine of the UCH domain) but
#                 cannot be assumed so without functional annotation.
#   HEMIZYGOUS    GISTIC shallow deletion (-1) and nothing else -- one copy lost, no
#                 second hit DETECTED.
#   WILD-TYPE     no mutation, no SV, GISTIC diploid or gain.
# Biallelic = ALTERED with a shallow deletion as well (mutation + loss of the other
# copy), reported as a sub-flag.
#
# Caveat that limits "WILD-TYPE": exome and SNP-array calls miss a known share of
# BAP1 inactivation in mesothelioma (small intragenic deletions below array
# resolution, deep intronic / splicing events, promoter silencing).  WILD-TYPE here
# means "no BAP1 lesion detected by these assays", not "BAP1 intact".  BAP1 mRNA is
# printed per class as a sanity check on the direction.
#
# BAP1 protein (RPPA) is used as an independent check on the calls -- it is the
# nearest TCGA equivalent of the IHC used clinically and in MESOMICS.
#
# Input : tcga_genetic/*.tsv, tcga_genetic/samplelist_*.txt,
#         bulkRNA_meso/bulk_RNA_studies(_metadata).rds
# Output: tcga_genetic/TCGA_BAP1_genetic_status.csv
###############################################################################
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUT  <- file.path(ROOT, "git_repo_claude", "R2_Q14", "tcga_genetic")
BULK <- file.path(ROOT, "bulkRNA_meso")
setwd(OUT)

rd <- function(f) read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
mut <- rd("BAP1_mutations.tsv"); cna <- rd("BAP1_gistic.tsv")
l2  <- rd("BAP1_log2CNA.tsv");   sv  <- rd("BAP1_structural_variants.tsv")
rp  <- rd("BAP1_rppa.tsv")
seqd <- readLines("samplelist_sequenced.txt"); svd <- readLines("samplelist_sv.txt")

## our 87 tumours, ids converted TCGA.XX.YYYY.01 -> TCGA-XX-YYYY-01
meta <- readRDS(file.path(BULK, "bulk_RNA_studies_metadata.rds"))$tcga
expr <- as.matrix(readRDS(file.path(BULK, "bulk_RNA_studies.rds"))$tcga)
ours <- colnames(expr); id <- gsub("\\.", "-", ours)

TRUNC <- c("Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation","Splice_Site",
           "Nonstop_Mutation","Translation_Start_Site")
VUS   <- c("Missense_Mutation","In_Frame_Del","In_Frame_Ins")
agg <- function(x) if (length(x)) paste(unique(x), collapse = "; ") else ""

S <- data.frame(sample = ours, tcga_id = id, stringsAsFactors = FALSE)
S$mutation_assayed <- S$tcga_id %in% seqd
S$sv_assayed       <- S$tcga_id %in% svd
S$mutations   <- vapply(S$tcga_id, function(s) agg(paste0(mut$proteinChange, " (", mut$mutationType, ", VAF ", mut$vaf, ")")[mut$sampleId == s]), "")
S$truncating  <- vapply(S$tcga_id, function(s) any(mut$mutationType[mut$sampleId == s] %in% TRUNC), NA)
S$missense_or_inframe <- vapply(S$tcga_id, function(s) any(mut$mutationType[mut$sampleId == s] %in% VUS), NA)
S$fusion      <- vapply(S$tcga_id, function(s) agg(sv$eventInfo[sv$sampleId == s]), "")
S$gistic      <- cna$alteration[match(S$tcga_id, cna$sampleId)]
S$log2CNA     <- l2$value[match(S$tcga_id, l2$sampleId)]
S$STATUS_3P   <- as.character(meta$STATUS_3P[match(ours, rownames(meta))])

S$class <- ifelse(S$truncating | nzchar(S$fusion) | S$gistic %in% -2, "ALTERED",
            ifelse(S$missense_or_inframe, "VUS",
             ifelse(S$gistic %in% -1, "HEMIZYGOUS", "WILD-TYPE")))
S$biallelic <- S$class == "ALTERED" & (S$gistic %in% -2 |
               ((S$truncating | nzchar(S$fusion)) & S$gistic %in% -1))
S$class[!S$mutation_assayed & S$class == "WILD-TYPE"] <- "WILD-TYPE (no mutation data)"
## NB a tumour without mutation data can still be classed from copy number or SV
S$BAP1_mRNA <- as.numeric(expr["BAP1", ours])
S$BAP1_protein_rppa <- rp$value[match(S$tcga_id, rp$sampleId)]

## the tertile proxy used so far, for comparison
q <- quantile(S$BAP1_mRNA, c(1/3, 2/3))
S$tertile_proxy <- ifelse(S$BAP1_mRNA <= q[1], "low", ifelse(S$BAP1_mRNA >= q[2], "high", "middle"))
write.csv(S, "TCGA_BAP1_genetic_status.csv", row.names = FALSE)

cat("=== genetic BAP1 status, TCGA MESO (n =", nrow(S), ") ===\n")
print(table(S$class))
cat("\nbiallelic among ALTERED:", sum(S$biallelic), "of", sum(S$class == "ALTERED"), "\n")
cat("\nevidence behind ALTERED:\n")
a <- S[S$class == "ALTERED", ]
cat("  truncating mutation:", sum(a$truncating), "| fusion:", sum(nzchar(a$fusion)),
    "| deep deletion:", sum(a$gistic %in% -2), "\n")
cat("\nGISTIC call, all tumours:\n"); print(table(factor(S$gistic, -2:2,
    c("deep del","shallow del","diploid","gain","amp"))))

cat("\n=== sanity check: BAP1 mRNA by genetic class ===\n")
print(do.call(rbind, lapply(split(S, S$class), function(d) data.frame(
  class = d$class[1], n = nrow(d), mean = round(mean(d$BAP1_mRNA), 2),
  median = round(median(d$BAP1_mRNA), 2)))), row.names = FALSE)
w <- wilcox.test(S$BAP1_mRNA[S$class == "ALTERED"], S$BAP1_mRNA[S$class == "WILD-TYPE"])
cat(sprintf("ALTERED vs WILD-TYPE BAP1 mRNA: Wilcoxon p = %.3g\n", w$p.value))

cat("\ntumours not in the mutation sample list:",
    paste(sprintf("%s (%s)", S$sample[!S$mutation_assayed], S$class[!S$mutation_assayed]),
          collapse = ", "), "\n")

## ---- protein: the readout closest to clinical IHC ------------------------------
## RPPA measures BAP1 protein.  If the genetic calls are right, ALTERED tumours must
## have less BAP1 protein -- the same check IHC gives in MESOMICS, on an independent
## assay from the sequencing that made the calls.
cat("\n=== BAP1 PROTEIN (RPPA) by genetic class ===\n")
cat("tumours with RPPA:", sum(!is.na(S$BAP1_protein_rppa)), "of", nrow(S), "\n")
P <- S[!is.na(S$BAP1_protein_rppa), ]
print(do.call(rbind, lapply(split(P, P$class), function(d) data.frame(
  class = d$class[1], n = nrow(d), mean = round(mean(d$BAP1_protein_rppa), 3),
  median = round(median(d$BAP1_protein_rppa), 3)))), row.names = FALSE)
wp <- suppressWarnings(wilcox.test(P$BAP1_protein_rppa[P$class == "ALTERED"],
                                   P$BAP1_protein_rppa[P$class == "WILD-TYPE"]))
cat(sprintf("ALTERED vs WILD-TYPE BAP1 protein: Wilcoxon p = %.3g\n", wp$p.value))
cp <- suppressWarnings(cor.test(P$BAP1_protein_rppa, P$BAP1_mRNA, method = "spearman"))
cat(sprintf("protein vs mRNA across tumours: Spearman rho = %.3f, p = %.3g\n",
            cp$estimate, cp$p.value))
## how well does each readout separate ALTERED from WILD-TYPE?  (AUC)
auc <- function(x, g) { r <- rank(x); n1 <- sum(g); n0 <- sum(!g)
  (sum(r[g]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
PW <- P[P$class %in% c("ALTERED","WILD-TYPE"), ]
cat(sprintf("AUC, low value predicts ALTERED (%d vs %d): protein %.3f | mRNA %.3f\n",
            sum(PW$class == "ALTERED"), sum(PW$class == "WILD-TYPE"),
            1 - auc(PW$BAP1_protein_rppa, PW$class == "ALTERED"),
            1 - auc(PW$BAP1_mRNA,         PW$class == "ALTERED")))

cat("\n=== does the mRNA tertile proxy recover the genetic call? ===\n")
print(table(genetic = S$class, tertile_proxy = S$tertile_proxy))
cat("\n=== does arm-level STATUS_3P recover it? ===\n")
print(table(genetic = S$class, STATUS_3P = S$STATUS_3P, useNA = "ifany"))

cat("\n=== ALTERED and VUS tumours, detail ===\n")
print(S[S$class %in% c("ALTERED","VUS"),
        c("sample","class","biallelic","mutations","fusion","gistic","BAP1_mRNA","BAP1_protein_rppa")],
      row.names = FALSE)
cat("\nDONE\n")
