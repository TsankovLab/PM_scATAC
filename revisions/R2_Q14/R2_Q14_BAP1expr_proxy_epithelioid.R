###############################################################
# R2_Q14 : BAP1 EXPRESSION as a proxy for BAP1 status, WITHIN EPITHELIOID.
#
# Bypass the (degenerate/underpowered) genetic calls: use continuous BAP1
# expression to split BAP1-high vs BAP1-low, restricted to EPITHELIOID tumors
# so histology is held ~constant, and ask which TFs differ.
#
#  PART 1  bulk RNA : within Epithelioid samples of each cohort ->
#          (a) TF expression vs BAP1 mRNA (Spearman, continuous proxy)
#          (b) TF DE, BAP1-low tertile vs BAP1-high tertile (wilcoxauc)
#  PART 2  scATAC   : within EPITHELIOID (low scS-score) samples ->
#          per-cell BAP1 gene score high vs low TF activity (chromVAR),
#          meta-analysed across the epithelioid samples.
###############################################################
suppressMessages({ library(ggplot2); library(ggrepel); library(presto) })

SC     <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
BULKDIR<- file.path(SC, "bulkRNA_meso")
TCOMP  <- file.path(SC, "tumor_compartment")
OUTDIR <- file.path(SC, "git_repo_claude/R2_Q14")
setwd(OUTDIR); dir.create("Plots", showWarnings=FALSE)
runx <- c("RUNX1","RUNX2","RUNX3","CBFB")

mt   <- read.csv("BAP1_TF_activity_master_table.csv", stringsAsFactors=FALSE)
tf_universe <- unique(toupper(mt$TF))

############################  PART 1 : BULK RNA  ###############################
bl   <- readRDS(file.path(BULKDIR,"bulk_RNA_studies.rds"))
meta <- readRDS(file.path(BULKDIR,"bulk_RNA_studies_metadata.rds"))
align_meta <- function(s){
  m <- bl[[s]]; md <- meta[[s]]
  if (all(colnames(m) %in% rownames(md))) md <- md[colnames(m),,drop=FALSE] else
    if ("Sample" %in% colnames(md)) md <- md[match(colnames(m), md$Sample),,drop=FALSE]
  md
}
bulk_cohort <- function(s){
  m <- bl[[s]]; md <- align_meta(s)
  if (!("BAP1" %in% rownames(m)) || is.null(md$subtype)) return(NULL)
  epi <- which(as.character(md$subtype)=="Epithelioid")
  if (length(epi) < 12) { cat("[",s,"] only",length(epi),"epithelioid, skip\n"); return(NULL) }
  m <- m[, epi, drop=FALSE]
  bap1 <- as.numeric(m["BAP1",])
  genes <- rownames(m)[toupper(rownames(m)) %in% tf_universe]
  cat("[",s,"] epithelioid n=",length(epi)," TFs=",length(genes),
      " BAP1 mRNA range",round(min(bap1),2),"-",round(max(bap1),2),"\n")
  ## (a) continuous: TF expr vs BAP1 mRNA (Spearman)
  cor_rows <- do.call(rbind, lapply(genes, function(g){
    x <- as.numeric(m[g,]); ok <- is.finite(x) & is.finite(bap1)
    if (sum(ok) < 4 || sd(x[ok])==0 || sd(bap1[ok])==0)
      return(data.frame(TF=toupper(g), rho=NA_real_, p=NA_real_, stringsAsFactors=FALSE))
    ct <- suppressWarnings(cor.test(x[ok], bap1[ok], method="spearman"))
    data.frame(TF=toupper(g), rho=unname(ct$estimate), p=ct$p.value, stringsAsFactors=FALSE)
  }))
  cor_rows$padj <- p.adjust(cor_rows$p, "BH"); cor_rows$cohort <- s
  ## (b) tertile split BAP1-low vs BAP1-high
  q <- quantile(bap1, c(1/3, 2/3), na.rm=TRUE)
  grp <- ifelse(bap1<=q[1], "BAP1low", ifelse(bap1>=q[2], "BAP1high", NA))
  keep <- !is.na(grp)
  de <- NULL
  if (sum(grp=="BAP1low",na.rm=TRUE)>=4 && sum(grp=="BAP1high",na.rm=TRUE)>=4){
    de <- wilcoxauc(as.matrix(m[genes, keep, drop=FALSE]), grp[keep])
    de <- de[de$group=="BAP1low", c("feature","auc","logFC","pval","padj")]  # auc>0.5 => higher in BAP1-low
    names(de) <- c("TF","auc_low","logFC_low","p_de","padj_de"); de$TF <- toupper(de$TF); de$cohort <- s
  }
  list(cor=cor_rows, de=de,
       n_epi=length(epi), n_low=sum(grp=="BAP1low",na.rm=TRUE), n_high=sum(grp=="BAP1high",na.rm=TRUE))
}
cohorts <- names(bl)
bres <- lapply(cohorts, bulk_cohort); names(bres) <- cohorts
bres <- bres[!sapply(bres, is.null)]

cor_all <- do.call(rbind, lapply(bres, `[[`, "cor"))
de_all  <- do.call(rbind, lapply(bres, function(x) x$de))
write.csv(cor_all, "epi_proxy_bulk_TFvsBAP1_spearman.csv", row.names=FALSE)
if (!is.null(de_all)) write.csv(de_all, "epi_proxy_bulk_TF_DE_lowVShigh.csv", row.names=FALSE)

cat("\n=========== BULK: TF expr vs BAP1 mRNA within EPITHELIOID (Spearman) ===========\n")
for (s in names(bres)){
  d <- cor_all[cor_all$cohort==s,]; d <- d[order(d$p),]
  cat("\n---",s,"--- (rho<0 = TF UP as BAP1 goes DOWN, i.e. up in BAP1-low) top 12\n")
  print(head(d[,c("TF","rho","p","padj")],12), row.names=FALSE, digits=3)
  cat("  #padj<0.05:", sum(d$padj<0.05,na.rm=TRUE),"\n")
}
cat("\n=== BULK RUNX within epithelioid (vs BAP1 mRNA) ===\n")
print(cor_all[cor_all$TF %in% runx, c("cohort","TF","rho","p","padj")], row.names=FALSE, digits=3)

## cross-cohort consistent (continuous): same-sign rho & nominally sig in >=2 cohorts
cons_b <- do.call(rbind, lapply(split(cor_all, cor_all$TF), function(x){
  x <- x[is.finite(x$rho) & is.finite(x$p), ]
  if (nrow(x)<2) return(NULL)
  data.frame(TF=x$TF[1], n=nrow(x), same_sign=length(unique(sign(x$rho)))==1,
             mean_rho=mean(x$rho), min_p=min(x$p), n_sig=sum(x$p<0.05), stringsAsFactors=FALSE)
}))
cons_b <- cons_b[cons_b$same_sign & cons_b$n_sig>=2,]; cons_b <- cons_b[order(-abs(cons_b$mean_rho)),]
write.csv(cons_b, "epi_proxy_bulk_consistentTFs.csv", row.names=FALSE)
cat("\n=== BULK cross-cohort consistent TFs (same-sign rho, sig in >=2 epithelioid cohorts) ===\n")
if (nrow(cons_b)) print(head(cons_b,25), row.names=FALSE, digits=3) else cat("  (none)\n")

############################  PART 2 : scATAC  ################################
# epithelioid = low scS-score scATAC samples; reuse per-cell within-sample AUCs.
ps  <- read.csv("scATAC_per_sample_RUNX_BAP1_sarc.csv", stringsAsFactors=FALSE)
ws  <- read.csv("within_sample_TFactivity_BAP1_high_vs_low.csv", stringsAsFactors=FALSE)
scS_med <- median(ps$sarc_score, na.rm=TRUE)
epi_samples <- ps$sample[is.finite(ps$sarc_score) & ps$sarc_score < scS_med]
epi_samples <- intersect(epi_samples, unique(ws$sample))
cat("\n=========== scATAC: epithelioid (low scS) samples used ===========\n")
cat("scS median:", round(scS_med,3), " | epithelioid samples:", paste(epi_samples, collapse=", "), "\n")
print(ps[ps$sample %in% epi_samples, c("sample","BAP1_status","BAP1_gs","sarc_score")],
      row.names=FALSE, digits=3)

we <- ws[ws$sample %in% epi_samples, ]
# meta across epithelioid samples: effect = mean(auc-0.5); consistency of direction
meta_sc <- do.call(rbind, lapply(split(we, we$TF), function(x){
  eff <- x$auc - 0.5
  data.frame(TF=x$TF[1], mean_eff=mean(eff, na.rm=TRUE),
             n_samples=nrow(x), n_pos=sum(eff>0), n_neg=sum(eff<0),
             n_sig=sum(x$padj<0.05, na.rm=TRUE), stringsAsFactors=FALSE)
}))
meta_sc$consistent <- pmax(meta_sc$n_pos, meta_sc$n_neg) == meta_sc$n_samples
meta_sc <- meta_sc[order(-abs(meta_sc$mean_eff)), ]
write.csv(meta_sc, "epi_proxy_scATAC_TFactivity_meta.csv", row.names=FALSE)
cat("\n=== scATAC within-epithelioid TF activity, BAP1-high vs low (meta over samples) ===\n")
cat("(mean_eff>0 = higher activity in BAP1-accessible/high cells; top 20 by |effect|)\n")
print(head(meta_sc[,c("TF","mean_eff","n_samples","n_pos","n_neg","n_sig","consistent")],20),
      row.names=FALSE, digits=3)
cat("\n=== scATAC RUNX within epithelioid ===\n")
print(meta_sc[meta_sc$TF %in% runx, c("TF","mean_eff","n_pos","n_neg","n_sig","consistent")],
      row.names=FALSE, digits=3)
cat("consistent-direction TFs (all epithelioid samples agree):", sum(meta_sc$consistent),"\n")

###############################  PLOTS  ######################################
## (1) bulk MESOMICS Spearman volcano (largest epithelioid cohort)
big <- names(which.max(sapply(bres, function(x) x$n_epi)))
db <- cor_all[cor_all$cohort==big,]; db$sig <- db$padj<0.05
db$lab <- ifelse(db$TF %in% runx | (db$sig & abs(db$rho)>quantile(abs(db$rho),.98,na.rm=TRUE)), db$TF, NA)
p1 <- ggplot(db, aes(rho, -log10(p))) +
  geom_point(aes(color=sig), size=.9) +
  scale_color_manual(values=c(`FALSE`="grey78",`TRUE`="#c0392b"), labels=c("ns","BH<0.05"), name=NULL) +
  geom_vline(xintercept=0,color="grey50") + geom_text_repel(aes(label=lab), size=2.3, max.overlaps=30) +
  theme_bw(base_size=10) +
  xlab("Spearman rho: TF expr vs BAP1 mRNA (rho<0 = up in BAP1-low)") + ylab("-log10 p") +
  ggtitle(sprintf("BULK %s: TF expression vs BAP1 mRNA, WITHIN EPITHELIOID", toupper(big)),
          sprintf("n=%d epithelioid tumors; BAP1 expression proxy; RUNX labelled", bres[[big]]$n_epi))
ggsave("Plots/epi_proxy_bulk_spearman_volcano.pdf", p1, width=7, height=5)

## (2) scATAC meta effect: top TFs + RUNX highlighted
ms <- meta_sc; ms$is_runx <- ms$TF %in% runx
top <- head(ms[order(-abs(ms$mean_eff)),], 25)
top <- rbind(top, ms[ms$is_runx & !(ms$TF %in% top$TF),])
top$TF <- factor(top$TF, levels=top$TF[order(top$mean_eff)])
p2 <- ggplot(top, aes(mean_eff, TF, fill=is_runx)) +
  geom_col() + geom_vline(xintercept=0,color="grey40") +
  scale_fill_manual(values=c(`FALSE`="#34495e",`TRUE`="#e74c3c"),
                    labels=c("other TF","RUNX"), name=NULL) +
  theme_bw(base_size=9) +
  xlab("mean(AUC-0.5) across epithelioid scATAC samples\n(>0 = higher in BAP1-high cells)") + ylab("") +
  ggtitle("scATAC: TF activity, BAP1-high vs low WITHIN EPITHELIOID samples",
          sprintf("meta over %d low-scS samples: %s", length(epi_samples), paste(epi_samples,collapse=", ")))
ggsave("Plots/epi_proxy_scATAC_meta_bars.pdf", p2, width=7, height=6)

## (3) cross-modal: does bulk (rho) agree with scATAC (mean_eff)?
xm <- merge(
  aggregate(rho~TF, cor_all, mean),
  meta_sc[,c("TF","mean_eff")], by="TF")
xmr <- suppressWarnings(cor.test(xm$rho, xm$mean_eff, method="spearman"))
xm$lab <- ifelse(xm$TF %in% runx, xm$TF, NA)
p3 <- ggplot(xm, aes(rho, mean_eff)) +
  geom_hline(yintercept=0,color="grey80") + geom_vline(xintercept=0,color="grey80") +
  geom_point(color="grey55", size=.8) +
  geom_point(data=subset(xm,!is.na(lab)), color="#e74c3c", size=2) +
  geom_text_repel(aes(label=lab), size=2.6) +
  theme_bw(base_size=10) +
  xlab("BULK epithelioid: mean rho (TF expr vs BAP1 mRNA)") +
  ylab("scATAC epithelioid: mean(AUC-0.5) BAP1-hi vs lo") +
  ggtitle("Do bulk-RNA and scATAC BAP1-proxy TF signals agree (within epithelioid)?",
          sprintf("Spearman rho=%.2f, p=%.2g", xmr$estimate, xmr$p.value))
ggsave("Plots/epi_proxy_crossmodal_bulk_vs_scATAC.pdf", p3, width=6, height=5)

cat("\nDONE\n")
