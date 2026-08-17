###############################################################
# R2_Q3 : map the P4 chr8q-amplified clone SPATIALLY (Visium, slides C1/D1).
# Within malignant (epithelioid-deconvolved) spots, score chr8q coverage
# (mean log-norm expression of chr8q genes, relative to normal reference spots),
# split into chr8q-amp / chr8q-low clones, plot on the slide, and test spatial
# clustering with Moran's I.
###############################################################
suppressMessages({ library(Seurat); library(EnsDb.Hsapiens.v86); library(ggplot2); library(ggforce); library(patchwork) })
set.seed(1)
SC  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUT <- file.path(SC,"git_repo_claude","R2_Q3"); setwd(OUT); dir.create("Plots", showWarnings=FALSE)
spa_path <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/meso_spatial/spa_all_magic.rds"
slides <- c("C1","D1"); CEN8 <- 44905425        # hg38 chr8 centromere, bp

## ---- STEP 1: define the chr8q gene set -------------------------------------
## Everything with its midpoint distal to the centromere counts as chr8q. This is the
## expression proxy for the arm-level CNV score used in the single-cell analyses.
edb <- EnsDb.Hsapiens.v86
ge <- genes(edb, columns=c("gene_name","seq_name","gene_seq_start","gene_seq_end"),
            filter=SeqNameFilter("8"), return.type="data.frame")
ge$mid <- (ge$gene_seq_start+ge$gene_seq_end)/2
chr8q_genes <- unique(ge$gene_name[ge$mid > CEN8 & ge$gene_name!=""])
cat("chr8q genes:", length(chr8q_genes), "\n")

## ---- STEP 2: plotting + statistics helpers ---------------------------------
## get_spot_coords : Visium spot centres in image pixel space. The y axis is flipped
##   (ih - imagerow) because image rows count downward while ggplot counts upward.
get_spot_coords <- function(srt){ img<-Images(srt)[1]; io<-srt@images[[img]]; ih<-nrow(io@image)
  co<-GetTissueCoordinates(srt,image=img)
  if("imagerow"%in%colnames(co)) data.frame(x=as.numeric(co$imagecol), y=ih-as.numeric(co$imagerow))
  else data.frame(x=as.numeric(co[,2]), y=ih-as.numeric(co[,1])) }
## spatial_base : the H&E image as a raster background, axes locked to pixel extent
spatial_base <- function(img,ih,iw) ggplot()+annotation_raster(as.raster(img),0,iw,0,ih)+
  scale_x_continuous(limits=c(0,iw),expand=c(0,0))+scale_y_continuous(limits=c(0,ih),expand=c(0,0))+
  theme_void()+theme(aspect.ratio=ih/iw)
## feat_gg : continuous value per spot, drawn as circles at the real spot radius;
##   sorted so high-value spots are drawn last and stay visible
feat_gg <- function(srt,vals,title,cols,psf=0.5){ io<-srt@images[[Images(srt)[1]]]; img<-io@image; ih<-nrow(img); iw<-ncol(img)
  r<-io@spot.radius*max(ih,iw)*psf; xy<-get_spot_coords(srt); df<-data.frame(x0=xy$x,y0=xy$y,value=as.numeric(vals))
  df<-df[order(df$value),]; spatial_base(img,ih,iw)+geom_circle(data=df,aes(x0=x0,y0=y0,r=r,fill=value),color=NA,alpha=.9)+
  scale_fill_gradientn(colours=cols,name="chr8q")+ggtitle(title)+theme(plot.title=element_text(size=10,hjust=.5,face="bold"),legend.key.height=unit(.4,"cm")) }
## class_gg : same, for a categorical spot label (clone / non-malignant)
class_gg <- function(srt,cls,title,pal,psf=0.5){ io<-srt@images[[Images(srt)[1]]]; img<-io@image; ih<-nrow(img); iw<-ncol(img)
  r<-io@spot.radius*max(ih,iw)*psf; xy<-get_spot_coords(srt); lev<-rev(names(pal))
  df<-data.frame(x0=xy$x,y0=xy$y,class=factor(cls,levels=lev)); df<-df[order(match(df$class,lev)),]
  spatial_base(img,ih,iw)+geom_circle(data=df,aes(x0=x0,y0=y0,r=r,fill=class),color=NA,alpha=.9)+
  scale_fill_manual(values=pal,name=NULL)+ggtitle(title)+theme(plot.title=element_text(size=10,hjust=.5,face="bold"),legend.text=element_text(size=8)) }
## moran : Moran's I of a spot-level value on a 6-nearest-neighbour, row-normalised
##   spatial weight matrix. Significance by 299 permutations of the values over fixed
##   positions, so the null keeps the tissue geometry and only breaks the value-location
##   link. I ~ 0 = salt-and-pepper, I -> 1 = contiguous territories.
moran <- function(x, xy){ d<-as.matrix(dist(xy)); diag(d)<-Inf
  W<-t(apply(d,1,function(r){ k<-order(r)[1:6]; w<-numeric(length(r)); w[k]<-1; w })); W<-W/pmax(rowSums(W),1)
  xc<-x-mean(x); I<-(length(x)/sum(W))*sum(W*outer(xc,xc))/sum(xc^2)
  perm<-replicate(299,{ xp<-sample(xc); (length(x)/sum(W))*sum(W*outer(xp,xp))/sum(xp^2) })
  c(I=I, p=mean(abs(perm)>=abs(I))) }

## ---- STEP 3: per slide, score chr8q and split the malignant spots ----------
spa <- readRDS(spa_path)[slides]
tab <- list(); plots <- list()
for (sl in slides){
  obj <- spa[[sl]]
  ## 3a. call each spot malignant or not from the deconvolution weights: malignant if
  ##     the epithelioid tumour fractions outweigh the summed normal fractions.
  ##     ("Ephitelioid-p" is misspelt in the object; keep the typo or the lookup fails.)
  pred <- GetAssayData(obj, assay="predictions", slot="data")
  mal <- as.numeric(pred["Epithelioid-1",])+as.numeric(pred["Epithelioid-2",])+as.numeric(pred["Ephitelioid-p",])
  nor <- as.numeric(pred["Fibroblasts",])+as.numeric(pred["Tcells",])+as.numeric(pred["Myeloid",])+as.numeric(pred["Endothelial",])+as.numeric(pred["Bcells",])
  is_mal <- mal>nor; is_ref <- !is_mal
  ## 3b. CP10K log-normalise the raw counts (Seurat 5 renamed slot= to layer=, so try both)
  cnt <- tryCatch(GetAssayData(obj,assay="Spatial",layer="counts"), error=function(e) GetAssayData(obj,assay="Spatial",slot="counts"))
  ln <- log1p(sweep(cnt,2,colSums(cnt),"/")*1e4)
  ## 3c. chr8q score = mean expression over the chr8q genes, then subtract the mean over
  ##     NORMAL spots so the score is relative to non-tumour tissue on the same slide --
  ##     the spatial analogue of the normal-cell reference used in the single-cell CNV.
  g8 <- intersect(chr8q_genes, rownames(ln))
  s8 <- colMeans(as.matrix(ln[g8,,drop=FALSE]))
  s8 <- s8 - mean(s8[is_ref])
  ## 3d. split the MALIGNANT spots in two on the chr8q score alone (k-means on the scaled
  ##     score). Targeted on purpose: genome-wide clustering dilutes a single-arm event.
  mi <- which(is_mal); km<-kmeans(scale(s8[mi]),2,nstart=25)$cluster
  amp <- which.max(tapply(s8[mi],km,mean)); clone_mal<-ifelse(km==amp,"chr8q-amp","chr8q-low")
  cls <- rep("non-malignant", length(s8)); cls[mi]<-clone_mal
  ## 3e. is the split spatially organised? Moran's I on the chr8q score, malignant spots
  ##     only -- a clone that is real tissue structure gives contiguous territories,
  ##     a clustering artefact gives I ~ 0.
  xy <- get_spot_coords(obj); mI <- moran(s8[mi], xy[mi,])
  tab[[sl]] <- data.frame(slide=sl, n_spots=length(s8), n_malignant=length(mi),
    frac_amp=round(mean(clone_mal=="chr8q-amp"),2),
    chr8q_amp=round(mean(s8[mi][clone_mal=="chr8q-amp"]),3), chr8q_low=round(mean(s8[mi][clone_mal=="chr8q-low"]),3),
    moranI=round(mI["I"],3), moran_p=signif(mI["p"],2))
  ## 3f. figures: blank the normal spots so the score panel shows tumour only
  s8p<-s8; s8p[is_ref]<-NA
  plots[[paste0(sl,"_score")]] <- feat_gg(obj, s8p, paste0(sl," chr8q score (malignant spots)"),
                                          c("#2166ac","#f7f7f7","#b2182b"))
  plots[[paste0(sl,"_clone")]] <- class_gg(obj, cls, paste0(sl," chr8q clone (Moran I=",round(mI["I"],2),", p=",signif(mI["p"],2),")"),
    c(`chr8q-amp`="#d62728",`chr8q-low`="#1f77b4",`non-malignant`="#e5e5e5"))
}
## ---- STEP 4: summary table + 2x2 figure (score | clone, one row per slide) --
TAB <- do.call(rbind, tab); rownames(TAB)<-NULL
write.csv(TAB, "P4_visium_chr8q_spatial.csv", row.names=FALSE)
cat("=== P4 Visium chr8q clone spatial summary ===\n"); print(TAB, row.names=FALSE)

fig <- (plots[["C1_score"]] | plots[["C1_clone"]]) / (plots[["D1_score"]] | plots[["D1_clone"]]) +
  plot_annotation(title="P4 Visium: chr8q-amplified clone spatial map (slides C1, D1)", theme=theme(plot.title=element_text(face="bold",size=13)))
ggsave("Plots/R2_Q3_P4_visium_chr8q.pdf", fig, width=13, height=11)
cat("\nDONE\n")
