###############################################################################
# R2_Q14 -- the TCGA MESO paper's BAP1-associated INFERRED TF ACTIVITY (supplement
#           tab 2B) against every TF result we have.
#
# Tab 2B is a t-test of TF activities inferred with the Osmanbeyoglu et al. 2017
# (Nat Commun 8:14249) framework, which models each tumour's expression (and
# proteomic) profile from TF-target priors and so estimates what a TF is DOING, not
# how much of its mRNA there is.  That is the reason IRF8 and EGR2 are strong there
# (p ~ 3e-9) while their mRNA shows nothing in our TCGA test.
#
# Conceptually the closest of our readouts is scATAC motif activity (chromVAR), which
# is also inferred from targets (motif-containing peaks) rather than from the TF's
# own locus.  All four readouts are compared:
#   TCGA TF mRNA, genetic label   TCGA_genetic_BAP1_TFonly_unadjusted.csv (t, unadjusted)
#   scATAC motif activity         BAP1_TF_limma_sarcadjusted.csv (t_unadjusted, t_adj)
#   scATAC gene score             BAP1_TFgenescore_limma.csv (t_gs, sarc-adjusted)
#   scATAC gene-score-filtered    BAP1_TF_gsfiltered_results.csv (t_motif)
#
# Sign convention of estimate_tf is not stated in the tab.  It is ASSUMED that
# estimate > 0 means higher activity in BAP1-inactivated tumours; the check printed
# below (agreement with the direction of our TCGA IRF9/IRF7 mRNA effects, which are
# higher in altered tumours) is supportive, not proof.
#
# Output: paper_TFactivity_vs_ours.csv, Plots/paper_TFactivity_vs_ours.pdf
###############################################################################
suppressPackageStartupMessages({ library(readxl); library(ggplot2); library(ggrepel) })
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(ROOT, "git_repo_claude", "R2_Q14"))
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))
X <- "21598290cd180804-sup-205173_2_supp_5073240_pg14sl.xlsx"

B <- as.data.frame(read_excel(X, sheet = "2B_bap1_inferred_TF_activity"), check.names = FALSE)
names(B)[1] <- "TF"
for (k in c("p_tf","p_adj_tf","estimate_tf")) B[[k]] <- as.numeric(B[[k]])
B <- B[order(B$p_tf), ]
cat("tab 2B:", nrow(B), "TFs with inferred activity\n")
for (a in c(0.001, 0.01, 0.05)) cat(sprintf("p_adj < %.3f : %d TFs (estimate > 0: %d, < 0: %d)\n", a,
    sum(B$p_adj_tf < a), sum(B$p_adj_tf < a & B$estimate_tf > 0), sum(B$p_adj_tf < a & B$estimate_tf < 0)))
cat("\ntop 25 by p:\n"); print(head(B, 25), row.names = FALSE, digits = 3)

ours <- list(
  tcga_mRNA   = { d <- read.csv("TCGA_genetic_BAP1_TFonly_unadjusted.csv", stringsAsFactors = FALSE)
                  data.frame(TF = d$TF, t = d$t, P = d$P) },
  atac_motif_unadj = { d <- read.csv("BAP1_TF_limma_sarcadjusted.csv", stringsAsFactors = FALSE)
                  data.frame(TF = d$TF, t = d$t_unadjusted, P = d$P_unadjusted) },
  atac_motif_adj   = { d <- read.csv("BAP1_TF_limma_sarcadjusted.csv", stringsAsFactors = FALSE)
                  data.frame(TF = d$TF, t = d$t_adj, P = d$P_adj) },
  atac_genescore   = { d <- read.csv("BAP1_TFgenescore_limma.csv", stringsAsFactors = FALSE)
                  data.frame(TF = d$TF, t = d$t_gs, P = d$P_gs) },
  atac_gs_filtered = { d <- read.csv("BAP1_TF_gsfiltered_results.csv", stringsAsFactors = FALSE)
                  data.frame(TF = d$TF, t = d$t_motif, P = d$P_motif) })

cat("\n=== sign-convention check (assumed estimate > 0 = higher in BAP1-inactivated) ===\n")
chk <- merge(B[B$TF %in% c("IRF9","IRF7","IRF1","IRF8","EGR2"), c("TF","estimate_tf","p_adj_tf")],
             ours$tcga_mRNA, by = "TF", all.x = TRUE)
print(chk, row.names = FALSE, digits = 3)

sig <- B$TF[B$p_adj_tf < 0.05]
C <- do.call(rbind, lapply(names(ours), function(k){
  m <- merge(B, ours[[k]], by = "TF")
  ct <- suppressWarnings(cor.test(m$estimate_tf, m$t, method = "spearman"))
  ms <- m[m$TF %in% sig, ]
  data.frame(our_readout = k, n_shared = nrow(m), rho = unname(ct$estimate), p = ct$p.value,
             pct_same_sign_all = 100 * mean(sign(m$estimate_tf) == sign(m$t)),
             n_paper_sig_shared = nrow(ms),
             pct_same_sign_paper_sig = if (nrow(ms)) 100 * mean(sign(ms$estimate_tf) == sign(ms$t)) else NA,
             n_paper_sig_same_sign_and_P05 = sum(sign(ms$estimate_tf) == sign(ms$t) & ms$P < 0.05),
             stringsAsFactors = FALSE)
}))
cat("\n=== congruence: paper inferred TF activity vs our readouts ===\n")
print(transform(C, rho = round(rho, 3), p = signif(p, 3), pct_same_sign_all = round(pct_same_sign_all),
                pct_same_sign_paper_sig = round(pct_same_sign_paper_sig)), row.names = FALSE)

W <- B
for (k in names(ours)) {
  o <- ours[[k]]; names(o) <- c("TF", paste0("t_", k), paste0("P_", k)); W <- merge(W, o, by = "TF", all.x = TRUE)
}
W <- W[order(W$p_tf), ]
write.csv(W, "paper_TFactivity_vs_ours.csv", row.names = FALSE)
cat("\n=== paper-significant TFs (p_adj < 0.05) with our t values ===\n")
print(W[W$p_adj_tf < 0.05, c("TF","estimate_tf","p_adj_tf","t_tcga_mRNA","P_tcga_mRNA",
                             "t_atac_motif_unadj","P_atac_motif_unadj","t_atac_genescore","P_atac_genescore")],
      row.names = FALSE, digits = 3)

pl <- lapply(c("atac_motif_unadj","tcga_mRNA","atac_genescore"), function(k){
  d <- W[!is.na(W[[paste0("t_", k)]]), ]
  r <- cor(d$estimate_tf, d[[paste0("t_", k)]], method = "spearman")
  d$lab <- ifelse(d$p_adj_tf < 0.001, d$TF, NA)
  ggplot(d, aes(estimate_tf, .data[[paste0("t_", k)]])) +
    geom_hline(yintercept = 0, linewidth = .2, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
    geom_point(colour = "grey65", size = 1.2) +
    geom_point(data = subset(d, p_adj_tf < 0.05), colour = "grey15", size = 1.5) +
    geom_text_repel(aes(label = lab), size = 2.1, colour = "grey15", max.overlaps = 30,
                    seed = 1, min.segment.length = 0) +
    gtheme_no_rot + xlab("paper: BAP1 effect on inferred TF activity (estimate)") +
    ylab(paste("ours: t,", gsub("_", " ", k))) +
    ggtitle(sprintf("n = %d TFs, rho = %.2f (dark = paper p_adj < 0.05)", nrow(d), r))
})
pdf("Plots/paper_TFactivity_vs_ours.pdf", width = 5.2, height = 4.3)
for (p in pl) print(p); dev.off()
cat("\nDONE\n")
