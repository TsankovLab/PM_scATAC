###############################################################
# R2_Q14 : bulk-cohort differential EXPRESSION of TFs by BAP1 status.
#
# Reviewer asks: beyond RUNX, which TFs are differentially expressed in
# BAP1-lost vs BAP1-retained tumors in the bulk cohorts?
#
# Genetic (DNA-level) BAP1-loss calls per cohort (same as _bulk_genetic_BAP1.R):
#   Bueno    - FISH.chrom3  ("del 3p"/"del 3" = Loss)   [near-degenerate: ~44/45 lost]
#   TCGA     - STATUS_3P    ("Lost" = Loss)              [only ~5 lost]
#   MESOMICS - IHC.BAP1     (NO = Loss, verified by BAP1 mRNA)  [powered: 63 vs 43]
#
# For every TF (chromVAR motif universe) present in each cohort we run a
# Wilcoxon DE (Loss vs Retained), OVERALL and WITHIN EPITHELIOID (histology
# held constant), plus a sarc-score-adjusted linear model. We then build a
# cross-cohort consistency table so "any other TF" can be answered directly.
###############################################################
suppressMessages({ library(ggplot2); library(ggrepel); library(presto) })

SC     <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
BULKDIR<- file.path(SC, "bulkRNA_meso")
TCOMP  <- file.path(SC, "tumor_compartment")
OUTDIR <- file.path(SC, "git_repo_claude/R2_Q14")
setwd(OUTDIR); dir.create("Plots", showWarnings=FALSE)
runx <- c("RUNX1","RUNX2","RUNX3","CBFB")

bl   <- readRDS(file.path(BULKDIR,"bulk_RNA_studies.rds"))
meta <- readRDS(file.path(BULKDIR,"bulk_RNA_studies_metadata.rds"))
mt   <- read.csv("BAP1_TF_activity_master_table.csv", stringsAsFactors=FALSE)
tf_universe <- unique(toupper(mt$TF))
cat("TF universe (chromVAR motifs):", length(tf_universe), "TFs\n")

align_meta <- function(s){
  m <- bl[[s]]; md <- meta[[s]]
  if (all(colnames(m) %in% rownames(md))) md <- md[colnames(m),,drop=FALSE] else
    if ("Sample" %in% colnames(md)) md <- md[match(colnames(m), md$Sample),,drop=FALSE]
  md
}
get_sarc <- function(s, md, m){
  if ("sarc_score" %in% colnames(md)) return(as.numeric(md$sarc_score))
  nmf <- read.csv(file.path(TCOMP,"scrna","cnmf_genelist_25_nfeat_5000.csv"), stringsAsFactors=FALSE)
  g <- intersect(head(nmf[["cNMF20"]][!is.na(nmf[["cNMF20"]])],20), rownames(m))
  colMeans(m[g,,drop=FALSE], na.rm=TRUE)
}

## ================= build genetic BAP1 status per cohort =================
build <- function(s){
  m <- bl[[s]]; md <- align_meta(s)
  status <- rep(NA_character_, ncol(m)); src <- NA_character_
  if (s=="mesomics") {
    ihc <- as.character(md$IHC.BAP1)
    bap1 <- as.numeric(m["BAP1",])
    mn <- tapply(bap1, ihc, mean, na.rm=TRUE)
    loss_lab <- names(which.min(mn[c("YES","NO")]))
    status[ihc==loss_lab] <- "Loss"
    status[ihc==setdiff(c("YES","NO"), loss_lab)] <- "Retained"
    src <- paste0("IHC.BAP1 (", loss_lab, "=Loss)")
  }
  if (s=="bueno") {
    f <- as.character(md$FISH.chrom3)
    status[grepl("del *3p|del *3\\b|del3", f, ignore.case=TRUE)] <- "Loss"
    status[f=="3"] <- "Retained"
    src <- "FISH.chrom3 (del 3p/3=Loss)"
  }
  if (s=="tcga") {
    st <- as.character(md$STATUS_3P)
    status[st=="Lost"]   <- "Loss"
    status[st=="Gained" | st=="Neutral" | st=="Not Called"] <- "Retained"
    src <- "STATUS_3P (Lost=Loss)"
  }
  subtype <- as.character(md$subtype)
  list(m=m, status=status, subtype=subtype, sarc=get_sarc(s,md,m), src=src)
}

## ================= full TF DE scan (overall + epithelioid) =================
saf_w <- function(xx, gg){
  ok <- is.finite(xx) & !is.na(gg); xx <- xx[ok]; gg <- gg[ok]
  if (length(unique(gg)) < 2 || min(table(gg)) < 2) return(NA_real_)
  tryCatch(suppressWarnings(wilcox.test(xx ~ gg)$p.value), error=function(e) NA_real_)
}
scan_cohort <- function(s){
  b <- build(s); m <- b$m; status <- b$status; sub <- b$subtype; ss <- b$sarc
  keep <- !is.na(status)
  cat("\n[",s,"] src:", b$src, " Loss:", sum(status=="Loss",na.rm=TRUE),
      " Retained:", sum(status=="Retained",na.rm=TRUE),
      " (epithelioid Loss:", sum(status=="Loss" & sub=="Epithelioid",na.rm=TRUE),
      "Retained:", sum(status=="Retained" & sub=="Epithelioid",na.rm=TRUE), ")\n")
  tf_in <- intersect(tf_universe, toupper(rownames(m)))
  # map back to actual rownames (case)
  rn_up <- toupper(rownames(m)); idx <- match(tf_in, rn_up)
  genes <- rownames(m)[idx]
  out <- do.call(rbind, lapply(seq_along(genes), function(i){
    g <- genes[i]; x <- as.numeric(m[g,])
    grp <- status
    de_ov  <- mean(x[grp=="Loss"],na.rm=TRUE) - mean(x[grp=="Retained"],na.rm=TRUE)
    p_ov   <- saf_w(x, grp)
    epi <- sub=="Epithelioid" & !is.na(sub)
    de_epi <- mean(x[epi & grp=="Loss"],na.rm=TRUE) - mean(x[epi & grp=="Retained"],na.rm=TRUE)
    p_epi  <- saf_w(x[epi], grp[epi])
    bin <- as.integer(grp=="Loss")
    lm_p <- tryCatch({cf <- summary(lm(x ~ bin + ss))$coefficients
                      if ("bin" %in% rownames(cf)) cf["bin","Pr(>|t|)"] else NA_real_},
                     error=function(e) NA_real_)
    data.frame(TF=toupper(g), delta_overall=de_ov, p_overall=p_ov,
               delta_epi=de_epi, p_epi=p_epi, lm_p_adjSarc=lm_p, stringsAsFactors=FALSE)
  }))
  out$padj_overall <- p.adjust(out$p_overall, "BH")
  out$padj_epi     <- p.adjust(out$p_epi, "BH")
  out$cohort <- s
  out
}
cohorts <- c("bueno","tcga","mesomics")
scan <- do.call(rbind, lapply(cohorts, scan_cohort))
write.csv(scan, "bulk_TF_DE_by_BAP1_percohort.csv", row.names=FALSE)

## per-cohort top hits
cat("\n================ TOP differentially-expressed TFs (overall) per cohort ================\n")
for (s in cohorts){
  d <- scan[scan$cohort==s & is.finite(scan$p_overall), ]
  d <- d[order(d$p_overall), ]
  cat("\n---", s, "--- (delta = Loss - Retained; +up in loss) top 12 by p_overall\n")
  print(head(d[,c("TF","delta_overall","p_overall","padj_overall","delta_epi","p_epi","lm_p_adjSarc")],12),
        row.names=FALSE, digits=3)
  cat("  # padj_overall<0.05:", sum(d$padj_overall<0.05,na.rm=TRUE),
      " | # padj_epi<0.05:", sum(d$padj_epi<0.05,na.rm=TRUE), "\n")
}

## RUNX specifically
cat("\n================ RUNX family across cohorts ================\n")
print(scan[scan$TF %in% runx, c("cohort","TF","delta_overall","p_overall","padj_overall","delta_epi","p_epi","lm_p_adjSarc")],
      row.names=FALSE, digits=3)

## ================= cross-cohort consistency (the "any other TF" answer) ==============
# powered cohorts only for the overall contrast: mesomics is the only balanced one,
# but we still ask which TFs move the SAME direction & are nominally sig in >=2 cohorts.
cons <- do.call(rbind, lapply(split(scan, scan$TF), function(x){
  x <- x[is.finite(x$delta_overall) & is.finite(x$p_overall), ]
  if (nrow(x) < 2) return(NULL)
  same_dir <- length(unique(sign(x$delta_overall)))==1
  data.frame(TF=x$TF[1], n_cohorts=nrow(x), consistent_dir=same_dir,
             mean_delta=mean(x$delta_overall), min_p=min(x$p_overall),
             n_nomsig=sum(x$p_overall<0.05,na.rm=TRUE),
             cohorts=paste(x$cohort,collapse=","), stringsAsFactors=FALSE)
}))
cons_sig <- cons[cons$consistent_dir & cons$n_nomsig>=2, ]
cons_sig <- cons_sig[order(-abs(cons_sig$mean_delta)), ]
write.csv(cons[order(cons$min_p),], "bulk_TF_DE_by_BAP1_crosscohort.csv", row.names=FALSE)
cat("\n================ TFs consistent-direction & nominally sig in >=2 cohorts ================\n")
if (nrow(cons_sig)) print(head(cons_sig,30), row.names=FALSE, digits=3) else cat("  (none)\n")

## ================= plots =================
# (1) MESOMICS is the powered cohort: volcano of all TFs, overall contrast
dm <- scan[scan$cohort=="mesomics" & is.finite(scan$p_overall), ]
dm$sig <- dm$padj_overall<0.05
dm$lab <- ifelse(dm$TF %in% runx | (dm$sig & abs(dm$delta_overall)>quantile(abs(dm$delta_overall),.97,na.rm=TRUE)),
                 dm$TF, NA)
pV <- ggplot(dm, aes(delta_overall, -log10(p_overall))) +
  geom_point(aes(color=sig), size=.9) +
  scale_color_manual(values=c(`FALSE`="grey78",`TRUE`="#c0392b"),
                     labels=c("ns","BH<0.05"), name=NULL) +
  geom_vline(xintercept=0, color="grey50") +
  geom_text_repel(aes(label=lab), size=2.3, max.overlaps=30) +
  theme_bw(base_size=10) +
  xlab("TF expression: BAP1-Loss - Retained  (>0 = up in loss)") + ylab("-log10 p") +
  ggtitle("MESOMICS: TF differential expression by genetic BAP1 status",
          "powered cohort (63 loss / 43 retained); RUNX labelled")
ggsave("Plots/bulk_TF_DE_by_BAP1_mesomics_volcano.pdf", pV, width=7, height=5)

# (2) RUNX delta across cohorts, overall vs epithelioid
rr <- scan[scan$TF %in% runx, c("cohort","TF","delta_overall","delta_epi")]
rl <- reshape(rr, varying=c("delta_overall","delta_epi"), v.names="delta",
              timevar="contrast", times=c("overall","epithelioid"), direction="long")
rl$cohort <- factor(rl$cohort, levels=cohorts); rl$TF <- factor(rl$TF, levels=runx)
pR <- ggplot(rl, aes(TF, delta, fill=contrast)) +
  geom_col(position="dodge") + geom_hline(yintercept=0,color="grey40") +
  facet_wrap(~cohort, nrow=1) +
  scale_fill_manual(values=c(overall="#7f8c8d", epithelioid="#e67e22"), name=NULL) +
  theme_bw(base_size=10) + ylab("RUNX expr: Loss - Retained") + xlab("") +
  ggtitle("RUNX expression by genetic BAP1 status (overall vs within-epithelioid)")
ggsave("Plots/bulk_TF_DE_RUNX_overall_vs_epi.pdf", pR, width=8, height=3.4)

cat("\nDONE\n")
