###############################################################################
# R2_Q14 -- the TCGA MESO paper's BAP1 status (supplement tab 2A) vs our genetic calls.
#
# Supplement: 21598290cd180804-sup-205173_2_supp_5073240_pg14sl.xlsx
#   2A_BAP1 status  per-tumour integrated BAP1 inactivation call (copy number,
#                   mutations, mRNA, RPPA, hit details)
# Ours: tcga_genetic/TCGA_BAP1_genetic_status.csv (cBioPortal PanCancer calls).
#
# Output: paper_BAP1_status_vs_ours.csv
###############################################################################
suppressPackageStartupMessages(library(readxl))
setwd("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R2_Q14")
X <- "21598290cd180804-sup-205173_2_supp_5073240_pg14sl.xlsx"

A <- as.data.frame(read_excel(X, sheet = "2A_BAP1 status"), check.names = FALSE)
cat("tab 2A:", nrow(A), "tumours | columns:", paste(sprintf("[%s]", names(A)), collapse = " "), "\n")
for (k in names(A)) {
  v <- A[[k]]
  if (!is.numeric(suppressWarnings(as.numeric(v))) || length(unique(v)) <= 25) {
    if (length(unique(v)) <= 25) { cat("\n--", k, "--\n"); print(table(v, useNA = "ifany")) }
  }
}
cat("\nhit_details for tumours with at least one hit:\n")
hc <- grep("hit", names(A), value = TRUE)
print(A[suppressWarnings(as.numeric(A[["#hits"]])) > 0, c(names(A)[1:2], hc)], row.names = FALSE)

G <- read.csv("tcga_genetic/TCGA_BAP1_genetic_status.csv", stringsAsFactors = FALSE)
names(A)[1] <- "tcga_id"
cat(sprintf("\ntumours: paper %d | ours %d | in both %d\n", nrow(A), nrow(G), sum(A$tcga_id %in% G$tcga_id)))
cat("ours but NOT in the paper table (", sum(!G$tcga_id %in% A$tcga_id), "):",
    paste(G$tcga_id[!G$tcga_id %in% A$tcga_id], collapse = ", "), "\n")
cat("classes of those missing:\n"); print(table(G$class[!G$tcga_id %in% A$tcga_id]))

M <- merge(A, G[, c("tcga_id","class","mutations","fusion","gistic","BAP1_mRNA","BAP1_protein_rppa")],
           by = "tcga_id", all.x = TRUE)
cat("\n=== paper BAP1_status vs our genetic class ===\n")
print(table(paper = M$BAP1_status, ours = M$class, useNA = "ifany"))
cat("\n=== paper Copy_number vs our GISTIC call ===\n")
print(table(paper = M$Copy_number, gistic = M$gistic, useNA = "ifany"))

num <- function(x) suppressWarnings(as.numeric(x))
cat("\n=== BAP1 mRNA and RPPA by PAPER status (their columns) ===\n")
for (k in intersect(c("mRNA", "RPPA"), names(M)))
  print(aggregate(num(M[[k]]) ~ M$BAP1_status, FUN = function(x) c(n = length(x), mean = round(mean(x), 3))))
write.csv(M, "paper_BAP1_status_vs_ours.csv", row.names = FALSE)
cat("\nDONE\n")
