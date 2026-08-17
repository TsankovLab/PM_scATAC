###############################################################
# R2_Q14 : GENETIC BAP1 loss in the bulk cohorts vs TF activity (expression).
#
# Genetic (DNA-level) BAP1-loss calls available:
#   Bueno    – FISH.chrom3  (chr3p deletion; "del 3p"/"del 3" = loss)
#   TCGA     – STATUS_3P    (chr3p arm SNP-array CNV; "Lost")  [+ per-sample
#              BAP1 copy-number segment from TCGA_CNV_hg38.rds if matchable]
#   MESOMICS – IHC.BAP1     (protein IHC; direction VERIFIED here via BAP1 mRNA)
#
# For each cohort: (a) does genetic BAP1 loss track epithelioid histology?
# (b) RUNX + genome-wide TF expression by genetic BAP1 status, OVERALL and
# WITHIN EPITHELIOID (histology-controlled), + sarc-score-adjusted lm.
###############################################################
suppressMessages({ library(ggplot2); library(presto); library(GenomicRanges) })

SC     <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
BULKDIR<- file.path(SC, "bulkRNA_meso")
TCOMP  <- file.path(SC, "tumor_compartment")
OUTDIR <- file.path(SC, "git_repo_claude/R2_Q14")
setwd(OUTDIR); dir.create("Plots", showWarnings=FALSE)
runx <- c("RUNX1","RUNX2","RUNX3","CBFB")

bl   <- readRDS(file.path(BULKDIR,"bulk_RNA_studies.rds"))
meta <- readRDS(file.path(BULKDIR,"bulk_RNA_studies_metadata.rds"))
# TF universe from the scATAC chromVAR motifs (to annotate "any other TF")
mt   <- read.csv("BAP1_TF_activity_master_table.csv", stringsAsFactors=FALSE)
tf_universe <- unique(toupper(mt$TF))

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
    # verify direction with BAP1 mRNA
    bap1 <- as.numeric(m["BAP1",])
    mn <- tapply(bap1, ihc, mean, na.rm=TRUE)
    cat("\n[MESOMICS] BAP1 mRNA by IHC.BAP1 value:\n"); print(round(mn,3))
    loss_lab <- names(which.min(mn[c("YES","NO")]))   # lower mRNA = loss
    cat("[MESOMICS] IHC value with lower BAP1 mRNA (=Loss):", loss_lab, "\n")
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
    status[st=="Gained" | st=="Neutral" | st=="Not Called"] <- "Retained"  # non-loss
    src <- "STATUS_3P (Lost=Loss)"
  }
  subtype <- as.character(md$subtype)
  data.frame(sample=colnames(m), cohort=s, BAP1_genetic=status, subtype=subtype,
             sarc_score=get_sarc(s,md,m), BAP1_mrna=as.numeric(m["BAP1",]),
             src=src, stringsAsFactors=FALSE)
}
tabs <- lapply(c("bueno","tcga","mesomics"), build); names(tabs) <- c("bueno","tcga","mesomics")

## ================= (a) does genetic loss track histology? =================
cat("\n=== genetic BAP1 status counts & subtype cross-tab ===\n")
hist_rows <- list()
for (s in names(tabs)) {
  d <- tabs[[s]]; d <- d[!is.na(d$BAP1_genetic), ]
  cat("\n[",s,"] source:", unique(d$src), " n:", nrow(d),
      " Loss:", sum(d$BAP1_genetic=="Loss"), " Retained:", sum(d$BAP1_genetic=="Retained"), "\n")
  tb <- table(d$BAP1_genetic, d$subtype); print(tb)
  fp <- tryCatch(suppressWarnings(fisher.test(tb, simulate.p.value=TRUE)$p.value), error=function(e) NA)
  epi_frac_loss <- mean(d$subtype[d$BAP1_genetic=="Loss"]=="Epithelioid", na.rm=TRUE)
  epi_frac_ret  <- mean(d$subtype[d$BAP1_genetic=="Retained"]=="Epithelioid", na.rm=TRUE)
  cat("Fisher p(status x subtype):", signif(fp,3),
      " | epithelioid frac: Loss=", round(epi_frac_loss,2), " Retained=", round(epi_frac_ret,2), "\n")
  hist_rows[[s]] <- data.frame(cohort=s, n=nrow(d), fisher_p=fp,
                               epi_frac_loss=epi_frac_loss, epi_frac_ret=epi_frac_ret)
}
write.csv(do.call(rbind,hist_rows), "bulk_genetic_BAP1_histology.csv", row.names=FALSE)

## ================= (b) RUNX by genetic BAP1: overall / epithelioid / adj =================
saf_w <- function(xx, gg){                    # robust Wilcoxon: NA if not testable
  ok <- is.finite(xx) & !is.na(gg); xx <- xx[ok]; gg <- gg[ok]
  if (length(unique(gg)) < 2 || min(table(gg)) < 2) return(NA_real_)
  tryCatch(suppressWarnings(wilcox.test(xx ~ gg)$p.value), error=function(e) NA_real_)
}
test_gene <- function(d, g){
  x <- d[[g]]
  if (is.null(x) || all(is.na(x)))
    return(data.frame(gene=g, delta_overall=NA, p_overall=NA, n_epi=NA,
                      delta_epi=NA, p_epi=NA, lm_p_adjSarc=NA))
  grp <- d$BAP1_genetic; ss <- d$sarc_score; sub <- d$subtype
  de  <- mean(x[grp=="Loss"],na.rm=TRUE) - mean(x[grp=="Retained"],na.rm=TRUE)
  ov  <- saf_w(x, grp)
  epi <- sub=="Epithelioid" & !is.na(sub)
  ep  <- saf_w(x[epi], grp[epi])
  de_epi <- mean(x[epi & grp=="Loss"],na.rm=TRUE) - mean(x[epi & grp=="Retained"],na.rm=TRUE)
  bin <- as.integer(grp=="Loss")
  lm_p <- tryCatch({cf <- summary(lm(x ~ bin + ss))$coefficients
                    if ("bin" %in% rownames(cf)) cf["bin","Pr(>|t|)"] else NA_real_},
                   error=function(e) NA_real_)
  data.frame(gene=g, delta_overall=de, p_overall=ov, n_epi=sum(epi),
             delta_epi=de_epi, p_epi=ep, lm_p_adjSarc=lm_p)
}
runx_rows <- list()
for (s in names(tabs)) {
  d <- tabs[[s]]; d <- d[!is.na(d$BAP1_genetic),]
  m <- bl[[s]]
  for (g in runx) d[[g]] <- if (g %in% rownames(m)) as.numeric(m[g, d$sample]) else NA
  rr <- do.call(rbind, lapply(runx, function(g) test_gene(d,g)))
  rr$cohort <- s; runx_rows[[s]] <- rr
}
runx_tab <- do.call(rbind, runx_rows)
write.csv(runx_tab, "bulk_genetic_RUNX_by_BAP1.csv", row.names=FALSE)
cat("\n=== RUNX by GENETIC BAP1 (delta = Loss - Retained; <0 = RUNX down in loss) ===\n")
print(runx_tab[,c("cohort","gene","delta_overall","p_overall","delta_epi","p_epi","lm_p_adjSarc")],
      row.names=FALSE, digits=3)

## ================= (c) genome-wide TF scan WITHIN EPITHELIOID =================
cat("\n=== within-epithelioid genome-wide diff by genetic BAP1 (TFs highlighted) ===\n")
de_list <- list()
for (s in names(tabs)) {
  d <- tabs[[s]]; d <- d[!is.na(d$BAP1_genetic) & !is.na(d$subtype) & d$subtype=="Epithelioid",]
  if (sum(d$BAP1_genetic=="Loss",na.rm=TRUE)<5 | sum(d$BAP1_genetic=="Retained",na.rm=TRUE)<5) { cat("[",s,"] too few epithelioid genetic calls, skip\n"); next }
  m <- bl[[s]][, d$sample, drop=FALSE]
  res <- wilcoxauc(as.matrix(m), d$BAP1_genetic)
  res <- res[res$group=="Loss", c("feature","auc","logFC","pval","padj")]
  res$cohort <- s; res$is_TF <- toupper(res$feature) %in% tf_universe
  de_list[[s]] <- res
}
if (length(de_list)) {
  de_all <- do.call(rbind, de_list)
  write.csv(de_all, "bulk_genetic_epithelioid_DE_allgenes.csv", row.names=FALSE)
  # cross-cohort consistent TFs: same direction in >=2 cohorts, TF, nominally sig
  tfd <- de_all[de_all$is_TF, ]
  cons <- do.call(rbind, lapply(split(tfd, tfd$feature), function(x){
    if (nrow(x)<2) return(NULL)
    dir_consistent <- length(unique(x$auc>0.5))==1
    data.frame(TF=x$feature[1], n_cohorts=nrow(x), consistent_dir=dir_consistent,
               mean_auc=mean(x$auc), min_p=min(x$pval), n_sig=sum(x$pval<0.05))
  }))
  cons <- cons[cons$consistent_dir & cons$n_sig>=2, ]
  cons <- cons[order(-abs(cons$mean_auc-0.5)),]
  write.csv(cons, "bulk_genetic_epithelioid_consistentTFs.csv", row.names=FALSE)
  cat("\nTFs differing by genetic BAP1 within epithelioid, consistent & sig in >=2 cohorts:\n")
  print(head(cons,30), row.names=FALSE, digits=3)
}

## ================= plots =================
# RUNX delta within epithelioid by cohort
rt <- runx_tab; rt$cohort <- factor(rt$cohort, levels=c("bueno","tcga","mesomics"))
pR <- ggplot(rt, aes(gene, delta_epi, fill=cohort)) +
  geom_col(position="dodge") + geom_hline(yintercept=0,color="grey40") +
  theme_bw(base_size=10) + ylab("RUNX expr: Loss - Retained (epithelioid only)") + xlab("") +
  ggtitle("RUNX vs GENETIC BAP1 loss, within epithelioid tumors")
ggsave("Plots/bulk_genetic_RUNX_epithelioid.pdf", pR, width=6, height=3.6)
cat("\nDONE\n")
