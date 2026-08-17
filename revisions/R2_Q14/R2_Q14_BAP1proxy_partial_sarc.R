###############################################################
# R2_Q14 : Which BAP1-proxy TFs survive after partialling out the
#          RESIDUAL sarcomatoid gradient WITHIN epithelioid?
#
# The within-epithelioid BAP1-expression proxy (epi_proxy_bulk_*) surfaced a
# TF signature (SMARCC1/HOXA up-with-BAP1; SOX11/BARX2/POU3F3 up-in-BAP1-low).
# Caveat: BAP1 mRNA may still track a residual sarcomatoid-within-epithelioid
# gradient. Here we compute the PARTIAL Spearman correlation of each TF vs
# BAP1 mRNA controlling for the cNMF20 sarcomatoid score, restricted to
# epithelioid tumors. TFs whose |partial rho| stays large & significant are
# genuinely sarc-INDEPENDENT BAP1 effects; TFs that collapse toward 0 were
# the residual-histology artefact.
###############################################################
suppressMessages({ library(ggplot2); library(ggrepel) })

SC     <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
BULKDIR<- file.path(SC, "bulkRNA_meso")
TCOMP  <- file.path(SC, "tumor_compartment")
OUTDIR <- file.path(SC, "git_repo_claude/R2_Q14")
setwd(OUTDIR); dir.create("Plots", showWarnings=FALSE)
runx <- c("RUNX1","RUNX2","RUNX3","CBFB")

mt   <- read.csv("BAP1_TF_activity_master_table.csv", stringsAsFactors=FALSE)
tf_universe <- unique(toupper(mt$TF))
cat("TF universe (chromVAR motifs):", length(tf_universe), "TFs\n")

bl   <- readRDS(file.path(BULKDIR,"bulk_RNA_studies.rds"))
meta <- readRDS(file.path(BULKDIR,"bulk_RNA_studies_metadata.rds"))

align_meta <- function(s){
  m <- bl[[s]]; md <- meta[[s]]
  if (all(colnames(m) %in% rownames(md))) md <- md[colnames(m),,drop=FALSE] else
    if ("Sample" %in% colnames(md)) md <- md[match(colnames(m), md$Sample),,drop=FALSE]
  md
}
# cNMF20 sarcomatoid program score (same definition as the genetic-scan script)
get_sarc <- function(md, m){
  if ("sarc_score" %in% colnames(md)) return(as.numeric(md$sarc_score))
  nmf <- read.csv(file.path(TCOMP,"scrna","cnmf_genelist_25_nfeat_5000.csv"), stringsAsFactors=FALSE)
  g <- intersect(head(nmf[["cNMF20"]][!is.na(nmf[["cNMF20"]])],20), rownames(m))
  colMeans(m[g,,drop=FALSE], na.rm=TRUE)
}

# ---- partial Spearman: cor(x, y | z), rank-based; returns rho + p (t, df=n-3) ----
pspearman <- function(x, y, z){
  ok <- is.finite(x) & is.finite(y) & is.finite(z)
  n <- sum(ok); if (n < 6) return(c(rho=NA_real_, p=NA_real_, n=n))
  rx <- rank(x[ok]); ry <- rank(y[ok]); rz <- rank(z[ok])
  if (sd(rx)==0 || sd(ry)==0 || sd(rz)==0) return(c(rho=NA_real_, p=NA_real_, n=n))
  rxy <- cor(rx,ry); rxz <- cor(rx,rz); ryz <- cor(ry,rz)
  den <- sqrt((1-rxz^2)*(1-ryz^2))
  if (!is.finite(den) || den==0) return(c(rho=NA_real_, p=NA_real_, n=n))
  pr <- (rxy - rxz*ryz)/den
  pr <- max(min(pr, 1), -1)
  tval <- pr*sqrt((n-3)/(1-pr^2))
  p <- 2*pt(-abs(tval), df=n-3)
  c(rho=pr, p=p, n=n)
}
# marginal Spearman (for the before/after comparison)
mspearman <- function(x, y){
  ok <- is.finite(x) & is.finite(y); n <- sum(ok)
  if (n < 4 || sd(x[ok])==0 || sd(y[ok])==0) return(c(rho=NA_real_, p=NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman"))
  c(rho=unname(ct$estimate), p=ct$p.value)
}

partial_cohort <- function(s){
  m <- bl[[s]]; md <- align_meta(s)
  if (!("BAP1" %in% rownames(m)) || is.null(md$subtype)) return(NULL)
  epi <- which(as.character(md$subtype)=="Epithelioid")
  if (length(epi) < 12) { cat("[",s,"] only",length(epi),"epithelioid, skip\n"); return(NULL) }
  m <- m[, epi, drop=FALSE]; md <- md[epi,,drop=FALSE]
  bap1 <- as.numeric(m["BAP1",])
  sarc <- get_sarc(md, m)
  bs   <- mspearman(bap1, sarc)   # BAP1-vs-sarc confound magnitude within epithelioid
  cat(sprintf("[ %s ] epithelioid n=%d | BAP1~sarc rho=%.2f (p=%.2g)\n",
              s, length(epi), bs["rho"], bs["p"]))
  genes <- rownames(m)[toupper(rownames(m)) %in% tf_universe]
  out <- do.call(rbind, lapply(genes, function(g){
    x <- as.numeric(m[g,])
    mg <- mspearman(x, bap1)               # marginal (already computed before, recomputed here)
    pg <- pspearman(x, bap1, sarc)         # partial, controlling sarc
    data.frame(TF=toupper(g),
               rho_marg=unname(mg["rho"]), p_marg=unname(mg["p"]),
               rho_part=unname(pg["rho"]), p_part=unname(pg["p"]),
               stringsAsFactors=FALSE)
  }))
  out$padj_marg <- p.adjust(out$p_marg, "BH")
  out$padj_part <- p.adjust(out$p_part, "BH")
  out$cohort <- s
  out$bap1_sarc_rho <- unname(bs["rho"])
  list(out=out, n_epi=length(epi), bap1_sarc=bs)
}

cohorts <- names(bl)
res <- lapply(cohorts, partial_cohort); names(res) <- cohorts
res <- res[!sapply(res, is.null)]
all <- do.call(rbind, lapply(res, `[[`, "out"))
write.csv(all, "epi_proxy_partial_sarc_percohort.csv", row.names=FALSE)

## ---- per-cohort: how much does partialling shrink the signal? ----
cat("\n=========== marginal vs partial (sarc-adjusted) Spearman, per cohort ===========\n")
for (s in names(res)){
  d <- all[all$cohort==s,]
  cat(sprintf("\n--- %s --- (n_epi=%d)  #padj_marg<0.05: %d  ->  #padj_part<0.05: %d\n",
              s, res[[s]]$n_epi, sum(d$padj_marg<0.05,na.rm=TRUE), sum(d$padj_part<0.05,na.rm=TRUE)))
  dd <- d[is.finite(d$p_part),]; dd <- dd[order(dd$p_part),]
  dd$survives <- dd$padj_part<0.05
  cat("  top 12 by partial p (rho_part<0 = up in BAP1-low even after sarc adj):\n")
  print(head(dd[,c("TF","rho_marg","padj_marg","rho_part","padj_part","survives")],12),
        row.names=FALSE, digits=3)
}

## ---- RUNX + the reported signature TFs, before/after ----
sig_tfs <- c("SMARCC1","HOXA1","HOXA2","HOXA3","HOXA4","HOXA5","ZNF589","MITF","TCF4","E2F8",
             "SOX11","BARX2","POU3F3","SIX4","ZNF219","FOXD2","PRDM16","ID1")
cat("\n=========== RUNX family: marginal vs partial ===========\n")
print(all[all$TF %in% runx, c("cohort","TF","rho_marg","p_marg","rho_part","p_part")],
      row.names=FALSE, digits=3)
cat("\n=========== reported BAP1-proxy signature TFs: marginal vs partial ===========\n")
st <- all[all$TF %in% sig_tfs, c("cohort","TF","rho_marg","padj_marg","rho_part","padj_part")]
st <- st[order(st$TF, st$cohort),]
print(st, row.names=FALSE, digits=3)

## ---- cross-cohort: which TFs stay same-sign & significant AFTER partialling ----
cons <- do.call(rbind, lapply(split(all, all$TF), function(x){
  x <- x[is.finite(x$rho_part) & is.finite(x$p_part),]
  if (nrow(x) < 2) return(NULL)
  data.frame(TF=x$TF[1], n=nrow(x),
             same_sign_part = length(unique(sign(x$rho_part)))==1,
             mean_rho_marg = mean(x$rho_marg, na.rm=TRUE),
             mean_rho_part = mean(x$rho_part, na.rm=TRUE),
             n_sig_marg = sum(x$padj_marg<0.05, na.rm=TRUE),
             n_sig_part = sum(x$padj_part<0.05, na.rm=TRUE),
             min_p_part = min(x$p_part, na.rm=TRUE),
             cohorts = paste(x$cohort, collapse=","), stringsAsFactors=FALSE)
}))
surv <- cons[cons$same_sign_part & cons$n_sig_part>=2,]
surv <- surv[order(-abs(surv$mean_rho_part)),]
write.csv(cons[order(cons$min_p_part),], "epi_proxy_partial_sarc_crosscohort.csv", row.names=FALSE)
cat("\n=========== TFs SURVIVING sarc-partialling (same-sign partial rho, padj_part<0.05 in >=2 cohorts) ===========\n")
if (nrow(surv)) print(head(surv,30), row.names=FALSE, digits=3) else cat("  (none)\n")

## ---- how many marginal hits collapse? per cohort attrition table ----
cat("\n=========== attrition: marginal-significant TFs that survive partialling (per cohort) ===========\n")
for (s in names(res)){
  d <- all[all$cohort==s,]
  msig <- d[which(d$padj_marg<0.05),]
  if (!nrow(msig)) { cat(sprintf("--- %s --- 0 marginal hits\n", s)); next }
  kept <- sum(msig$padj_part<0.05, na.rm=TRUE)
  cat(sprintf("--- %s --- marginal-sig: %d  ->  still sig after sarc adj: %d (%.0f%%)\n",
              s, nrow(msig), kept, 100*kept/nrow(msig)))
}

###############################  PLOTS  ######################################
big <- names(which.max(sapply(res, function(x) x$n_epi)))
db <- all[all$cohort==big,]

## (1) marginal rho vs partial rho — points off the diagonal are sarc-driven; on the
##     diagonal = sarc-independent BAP1 effect.
db$lab <- ifelse(db$TF %in% c(sig_tfs,runx) | (db$padj_part<0.05 & abs(db$rho_part)>quantile(abs(db$rho_part),.985,na.rm=TRUE)),
                 db$TF, NA)
db$class <- ifelse(db$padj_part<0.05, "survives sarc-adj",
                   ifelse(db$padj_marg<0.05, "lost after sarc-adj", "ns"))
p1 <- ggplot(db, aes(rho_marg, rho_part)) +
  geom_abline(slope=1, intercept=0, color="grey70", linetype=2) +
  geom_hline(yintercept=0, color="grey85") + geom_vline(xintercept=0, color="grey85") +
  geom_point(aes(color=class), size=1) +
  scale_color_manual(values=c(`survives sarc-adj`="#c0392b",`lost after sarc-adj`="#2980b9",`ns`="grey80"), name=NULL) +
  geom_text_repel(aes(label=lab), size=2.3, max.overlaps=40) +
  theme_bw(base_size=10) +
  xlab("marginal Spearman rho (TF vs BAP1 mRNA)") +
  ylab("partial rho | cNMF20 sarc score") +
  ggtitle(sprintf("%s epithelioid: which BAP1-proxy TFs survive sarc-partialling?", toupper(big)),
          sprintf("on the diagonal = sarc-independent BAP1 effect; collapse toward 0 = residual histology (n=%d)", res[[big]]$n_epi))
ggsave("Plots/epi_proxy_partial_marg_vs_partial.pdf", p1, width=7.5, height=6)

## (2) signature TFs: marginal vs partial rho, faceted by cohort
sd <- all[all$TF %in% c(sig_tfs, runx),]
sd$TFf <- factor(sd$TF, levels=c(sig_tfs, runx))
sl <- reshape(sd[,c("cohort","TF","rho_marg","rho_part")],
              varying=c("rho_marg","rho_part"), v.names="rho",
              timevar="type", times=c("marginal","partial(sarc-adj)"), direction="long")
sl$TF <- factor(sl$TF, levels=rev(c(sig_tfs, runx)))
p2 <- ggplot(sl, aes(rho, TF, color=type)) +
  geom_vline(xintercept=0, color="grey60") +
  geom_point(size=2, alpha=.85) +
  facet_wrap(~cohort, nrow=1) +
  scale_color_manual(values=c(marginal="#95a5a6",`partial(sarc-adj)`="#c0392b"), name=NULL) +
  theme_bw(base_size=9) +
  xlab("Spearman rho vs BAP1 mRNA (within epithelioid)") + ylab("") +
  ggtitle("BAP1-proxy signature TFs: does the correlation survive sarc-partialling?")
ggsave("Plots/epi_proxy_partial_signature_dumbbell.pdf", p2, width=9, height=6)

cat("\nDONE\n")
