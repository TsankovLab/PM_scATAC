###############################################################################
# R2_Q3 -- shared settings and helpers for the epiAneufinder subclone analysis.
#
# Every numbered script starts with  source("00_common.R").  Nothing here runs an
# analysis; it only fixes the paths, the parameters and the two pieces of logic that
# more than one script needs (the clone split and the ARI), so that the clone
# definition cannot drift between the figures and the statistics.
###############################################################################

## ---- paths ------------------------------------------------------------------
SC     <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
ROOT   <- file.path(SC, "git_repo_claude", "R2_Q3")     # this folder; all outputs land here
ARCHR  <- file.path(SC, "tumor_compartment", "scatac_ArchR")   # malignant-cell scATAC project
SPATIAL<- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/meso_spatial/spa_all_magic.rds"
BLACKLIST <- paste0("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Public_data/",
                    "blacklisted_regions/ENCODE_blacklist/hg38-blacklist.v2.bed")
setwd(ROOT)
dir.create("Plots", showWarnings = FALSE)

## ---- epiAneufinder parameters (identical for every sample) -------------------
WINDOW   <- 5e6    # window size. 10 Mb leaves only ~100 windows genome-wide and the
                   # caller errors out; 100 kb is too sparse (99.4% of calls "normal").
MINFRAGS <- 5000   # per-cell fragment floor. The package default (20000) would keep
                   # 154 of P4's 3082 malignant cells.
NCORES   <- 8
KVAL     <- 4      # segmentation depth: epiAneufinder finds 2^k segments per chromosome
                   # (k = 4 -> up to 16), by recursive Anderson-Darling breakpoint search
GENOME   <- "BSgenome.Hsapiens.UCSC.hg38"
EXCLUDE  <- c("chrX", "chrY", "chrM")

## ---- clone calling ----------------------------------------------------------
K_CLONES <- 2      # epiAneufinder does NOT estimate the number of subclones --
                   # split_subclones(tree_depth = k) is literally cutree(k = k). We cut
                   # every sample's own dendrogram at 2, the depth that isolates P4's
                   # chr8q clone, and apply the same depth everywhere for comparability.
MIN_CLONE <- function(n) max(20, 0.05 * n)   # a clone must hold >=20 cells and >=5%

## the 9 tumours with enough malignant scATAC cells to run (P3 n=21 and P13 n=18 are
## extracted but not analysed)
SAMPLES <- c("P1", "P4", "P5", "P8", "P10", "P11", "P12", "P14", "P23")

## ---- hg38 centromeres, Mb (for naming an arm from a bin coordinate) ----------
CENmb <- c(chr1=123.6, chr2=93.1, chr3=91.1, chr4=50.4, chr5=48.3, chr6=59.2, chr7=59.5,
           chr8=44.9, chr9=44.4, chr10=40.6, chr11=52.7, chr12=35.9, chr13=17.0,
           chr14=17.1, chr15=19.1, chr16=36.3, chr17=23.8, chr18=17.7, chr19=25.8,
           chr20=28.0, chr21=11.9, chr22=14.5)
armof <- function(chr, mb) paste0(sub("chr", "", chr), ifelse(mb < CENmb[chr], "p", "q"))
CEN8 <- 44.9e6                                  # chr8 centromere, bp -- chr8q = start > CEN8

## ---- helpers ----------------------------------------------------------------
## Path to one sample's epiAneufinder result table.
epi_table <- function(S) sprintf("out_5Mb/%s/epiAneufinder_results/results_table.tsv", S)

## Read one sample's calls.  Returns
##   M   : bins x cells integer matrix, 0 = loss, 1 = normal, 2 = gain
##   bin : "chr8:44900000" style key per row, in genome order
##   d   : the raw table (seq/start/end kept for coordinates)
read_epi <- function(S){
  d <- data.table::fread(epi_table(S)); data.table::setnames(d, 1, "idx")
  cc <- setdiff(names(d), c("idx", "seq", "start", "end"))
  M <- as.matrix(d[, ..cc]); colnames(M) <- sub("^cell-", "", cc)
  list(M = M, bin = paste0(d$seq, ":", d$start), d = d)
}

## The clone split, defined in exactly one place.
## Hierarchical clustering (ward.D) on the Euclidean distance between cells over their
## 0/1/2 call vectors -- this is what epiAneufinder's own split_subclones() does -- cut
## at K_CLONES, then drop clusters below the size floor.
## Returns a named character vector cell -> "c1"/"c2" over the retained cells only.
epi_clone_split <- function(S){
  e <- read_epi(S)
  ct <- t(e$M)                                  # cells x bins
  dm <- dist(ct); dm[is.na(dm)] <- 0
  cl <- cutree(hclust(dm, "ward.D"), k = K_CLONES)
  names(cl) <- colnames(e$M)
  keep <- names(table(cl))[table(cl) >= MIN_CLONE(ncol(e$M))]
  cl <- cl[as.character(cl) %in% keep]
  setNames(paste0("c", cl), names(cl))
}

## Adjusted Rand index between two labellings of the same cells.
ARI <- function(a, b){
  tb <- table(a, b); n <- sum(tb)
  ci <- sum(choose(rowSums(tb), 2)); cj <- sum(choose(colSums(tb), 2))
  ex <- ci * cj / choose(n, 2)
  (sum(choose(tb, 2)) - ex) / (0.5 * (ci + cj) - ex)
}

## ---- shared palette ---------------------------------------------------------
SAMPCOL <- setNames(colorRampPalette(c("#4c78a8","#f58518","#54a24b","#b279a2","#e45756",
                                       "#72b7b2","#ff9da6","#9d755d","#bab0ac"))(length(SAMPLES)),
                    SAMPLES)
COL_HI <- "#1baf7a"    # the P4 chr8q clone, everywhere it is highlighted

set.seed(1)
