#!/usr/bin/env Rscript
# R1_Q3_SOX9_myeloid_IL33_spatial.R
# Do SOX9-high malignant spots that are surrounded by Myeloid cells also express
# IL33? Reuses the proximity framework (spa_all_magic.rds; dominant cell type from
# the deconvolution `predictions`; Epithelioid*->Malignant; SOX9_high = top 50%
# SOX9 within malignant). For each malignant spot we compute distance to the
# nearest Myeloid spot (spot-spacing units) -> "myeloid-near" if <= 1.5 spacings.
# IL33 read from a RAW assay if available (MAGIC smooths signal across neighbours
# and would bleed myeloid IL33 into adjacent spots), with MAGIC IL33 as secondary.
# Tests: 2x2 SOX9(high/low) x myeloid(near/far) on IL33; within SOX9-high, IL33
# near-vs-far; and IL33 ~ myeloid proximity correlation.

suppressPackageStartupMessages({library(Seurat); library(ggplot2)})
projdir <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3'
spa <- readRDS('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/meso_spatial/spa_all_magic.rds')
epi_types <- c('Epithelioid-1','Epithelioid-2','Ephitelioid-p','Epithelioid-HBB')
NEAR_THR <- 1.5    # spot-spacing units => immediate myeloid neighbourhood

get_spot_coords <- function(obj) {
  img <- obj@images[[Images(obj)[1]]]; img_h <- nrow(img@image)
  co <- GetTissueCoordinates(obj, image=Images(obj)[1])
  if ('imagerow' %in% colnames(co)) cbind(x=as.numeric(co$imagecol), y=img_h-as.numeric(co$imagerow))
  else cbind(x=as.numeric(co[,2]), y=img_h-as.numeric(co[,1]))
}
min_dist <- function(q, t) if(nrow(t)==0) rep(NA_real_,nrow(q)) else
  apply(q,1,function(p) min(sqrt((t[,1]-p[1])^2+(t[,2]-p[2])^2)))
# IL33 from a raw assay if possible
il33_vec <- function(obj){
  pref <- intersect(c('Spatial','SCT','RNA'), Assays(obj))
  for (a in c(pref, setdiff(Assays(obj),'predictions'))){
    d <- tryCatch(GetAssayData(obj,assay=a,slot='data'), error=function(e) NULL)
    if (!is.null(d) && 'IL33' %in% rownames(d)) return(list(v=as.numeric(d['IL33',]), assay=a))
  }
  list(v=rep(NA,ncol(obj)), assay=NA)
}

rows <- list(); spotrows <- list()
for (sl in names(spa)) {
  obj <- spa[[sl]]
  pr  <- tryCatch(GetAssayData(obj,assay='predictions',slot='data'), error=function(e) NULL)
  if (is.null(pr)) { message(sl,": no predictions, skip"); next }
  pr  <- pr[rownames(pr)!='max',,drop=FALSE]
  dom <- ifelse(rownames(pr)[apply(pr,2,which.max)] %in% epi_types,'Malignant',
                rownames(pr)[apply(pr,2,which.max)])
  if (!('Myeloid' %in% dom) || sum(dom=='Malignant')<20) { message(sl,": no myeloid/malig, skip"); next }
  if (!('MAGIC_to_impute' %in% Assays(obj))) { message(sl,": no MAGIC assay, skip"); next }
  magic <- GetAssayData(obj,assay='MAGIC_to_impute',slot='data')
  if (!('SOX9' %in% rownames(magic))) { message(sl,": SOX9 not in MAGIC, skip"); next }
  sox9  <- as.numeric(magic['SOX9',])
  il    <- il33_vec(obj); il33 <- il$v
  il33_magic <- if ('IL33' %in% rownames(magic)) as.numeric(magic['IL33',]) else rep(NA,ncol(obj))
  coords <- get_spot_coords(obj)
  is_malig <- dom=='Malignant'
  thr <- quantile(sox9[is_malig],0.50)
  grp <- ifelse(sox9>=thr,'SOX9_high','SOX9_low'); grp[!is_malig]<-NA
  # spacing
  ns <- min(300,nrow(coords)); sd <- as.matrix(dist(coords[sample(nrow(coords),ns),])); diag(sd)<-Inf
  spacing <- median(apply(sd,1,min))
  myelo <- coords[dom=='Myeloid',,drop=FALSE]
  mdist <- rep(NA_real_,ncol(obj)); mdist[is_malig] <- min_dist(coords[is_malig,,drop=FALSE],myelo)/spacing
  near  <- ifelse(mdist<=NEAR_THR,'myeloid_near','myeloid_far')
  cat(sprintf("\n%s: malig=%d (SOX9hi=%d) myeloid=%d  IL33 assay=%s\n",
              sl,sum(is_malig),sum(grp=='SOX9_high',na.rm=TRUE),nrow(myelo),il$assay))
  m <- is_malig & !is.na(near)
  spotrows[[sl]] <- data.frame(slide=sl, SOX9=grp[m], myeloid=near[m], mdist=mdist[m],
                               IL33=il33[m], IL33_magic=il33_magic[m], sox9_expr=sox9[m])
}
spots <- do.call(rbind, spotrows)
write.csv(spots, file.path(projdir,'SOX9_myeloid_IL33_spots.csv'), row.names=FALSE)

## ---- group summary (raw IL33) ----
useIL <- if (all(is.na(spots$IL33))) 'IL33_magic' else 'IL33'
cat(sprintf("\n=== Using %s for IL33 ===\n", useIL))
spots$IL <- spots[[useIL]]
agg <- aggregate(IL ~ slide+SOX9+myeloid, spots, function(x) c(mean=mean(x),median=median(x),n=length(x)))
agg <- do.call(data.frame, agg); names(agg)[4:6]<-c('mean','median','n')
cat("\n=== IL33 by SOX9 x myeloid (per slide) ===\n"); print(agg)
write.csv(agg, file.path(projdir,'SOX9_myeloid_IL33_groupmeans.csv'), row.names=FALSE)

## ---- key tests ----
cat("\n=== tests (pooled across slides) ===\n")
hi <- spots[spots$SOX9=='SOX9_high',]
t1 <- wilcox.test(IL~myeloid, hi)
cat(sprintf("Within SOX9-high: IL33 near vs far  -> near median=%.3f far median=%.3f  p=%.3g\n",
            median(hi$IL[hi$myeloid=='myeloid_near']), median(hi$IL[hi$myeloid=='myeloid_far']), t1$p.value))
hi_near <- spots$SOX9=='SOX9_high' & spots$myeloid=='myeloid_near'
t2 <- wilcox.test(spots$IL[hi_near], spots$IL[!hi_near])
cat(sprintf("SOX9-high & myeloid-near vs all others -> %.3f vs %.3f  p=%.3g\n",
            median(spots$IL[hi_near]), median(spots$IL[!hi_near]), t2$p.value))
cs <- cor.test(hi$IL, hi$mdist, method='spearman')
cat(sprintf("Within SOX9-high: IL33 ~ myeloid distance (Spearman) rho=%.3f p=%.3g (neg=near->highIL33)\n",
            cs$estimate, cs$p.value))

## ---- plot: IL33 across the 4 groups, per slide ----
spots$grp4 <- factor(paste(spots$SOX9, spots$myeloid), levels=c(
  'SOX9_high myeloid_near','SOX9_high myeloid_far','SOX9_low myeloid_near','SOX9_low myeloid_far'))
p <- ggplot(spots, aes(grp4, IL, fill=grp4)) +
  geom_violin(scale='width', alpha=.8, linewidth=.3) + geom_boxplot(width=.15, outlier.shape=NA, fill='white', linewidth=.3) +
  facet_wrap(~slide, scales='free_y') +
  scale_fill_manual(values=c('#D62728','#F4A6A6','#5B8FF9','#C7D7F5'), guide='none') +
  labs(title=sprintf('IL33 (%s) in malignant spots by SOX9 status x Myeloid proximity', useIL),
       subtitle='myeloid-near = within 1.5 spot-spacings of nearest Myeloid spot',
       x=NULL, y='IL33 expression') +
  theme_bw(base_size=10) + theme(axis.text.x=element_text(angle=25,hjust=1,size=7))
ggsave(file.path(projdir,'SOX9_myeloid_IL33.pdf'), p, width=9, height=4.5)
cat("\nWrote spots csv + groupmeans + plot\n")
