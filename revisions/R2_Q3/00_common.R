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
## epiAneufinder does NOT estimate the number of subclones: split_subclones(tree_depth)
## is literally cutree(k = tree_depth), so asking for 2 gives 2 in every tumour whether
## or not there is anything there.  Instead we apply ONE cohort-wide criterion:
##
##   a tumour is divided in two only if its primary branch point separates cells that
##   differ by at least ARM_MIN in mean copy-number profile on at least one chromosome
##   arm; otherwise the tumour is reported as a single clone.
##
## The criterion is arm-level on purpose.  Genome-wide distance does NOT discriminate --
## the RMS profile difference of the primary split is 0.12-0.19 in all nine tumours and
## P4, the one subclone with orthogonal validation, sits at the BOTTOM of that range
## (0.117).  What sets P4 apart is that its difference is concentrated on one arm rather
## than spread thinly over the genome, and averaging within an arm is also what removes
## the per-bin sampling noise that shallow cells generate (see README, "why arm level").
LINKAGE      <- "ward.D2"   # ward.D2 is Ward's criterion on unsquared distances, which
                            # is what dist() returns; ward.D assumes squared input
ARM_MIN      <- 0.55        # calibrated on P4: its chr8q split scores 0.571, the weakest
                            # separation we accept, so nothing is called that is less
                            # clear-cut than the clone the Visium data validate (step 9).
                            # The next tumour below it is P8 at 0.465.
MIN_BINS_ARM <- 5           # arms with fewer informative 5 Mb bins are not tested: their
                            # means are dominated by one or two bins (19p, 22q, ...)
MIN_CLONE <- function(n) max(20, 0.05 * n)   # a clone must hold >=20 cells and >=5%
N_NULL    <- 20             # random splits used to report the noise floor of ARM_SEP

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

## Arm-level separation between two sets of cells: the mean difference of their CNV
## profiles on the arm where that difference is largest.  Profile = (fraction of cells
## called GAIN) - (fraction called LOSS), so the value is on a -1..1 scale and 0.55 means
## "on this arm, 55% more of one clone's cells carry the gain than of the other's".
## Returns the signed value, named after the arm (e.g. `8q`).
arm_separation <- function(X, arm, a, b, ok){
  dl <- rowMeans(X[, a, drop = FALSE]) - rowMeans(X[, b, drop = FALSE])
  ad <- vapply(split(dl, arm), mean, numeric(1))[ok]
  ad[which.max(abs(ad))]
}

## The clone split, defined in exactly one place.
##
##   1. hierarchical clustering (LINKAGE) of the cells on the Euclidean distance between
##      their 0/1/2 call vectors, normalised to a per-bin RMS so the scale is the same in
##      every tumour -- this is the tree epiAneufinder's own split_subclones() builds
##   2. take the primary branch point (k = 2)
##   3. accept it as a clone boundary only if both branches clear MIN_CLONE and the two
##      branches differ by >= ARM_MIN on some chromosome arm; otherwise the tumour is one
##      clone.  The same threshold is applied to all nine tumours.
##
## Returns a named character vector cell -> "c1"/"c2" covering EVERY cell of the sample
## (c1 = the major clone, c2 = the minor one, so "c2 - c1" always reads as "minor clone
## relative to the bulk of the tumour").  attr(, "stats") carries one row describing the
## decision: the separation, the arm it is on, the mean separation of N_NULL random
## splits of the same sizes (the noise floor), and whether the split was accepted.
epi_clone_split <- function(S){
  e   <- read_epi(S)
  X   <- (e$M == 2) - (e$M == 0); storage.mode(X) <- "numeric"   # profile scale, -1..1
  arm <- armof(e$d$seq, e$d$start / 1e6)
  ok  <- names(which(table(arm) >= MIN_BINS_ARM))
  cn  <- colnames(e$M); n <- length(cn)

  dm <- dist(t(X)); dm[is.na(dm)] <- 0
  dm <- dm / sqrt(nrow(X))                      # RMS per-bin call difference, 0..2
  g  <- cutree(hclust(dm, LINKAGE), k = 2)
  A  <- cn[g == 1]; B <- cn[g == 2]
  if (length(A) < length(B)){ tmp <- A; A <- B; B <- tmp }   # A = major, B = minor

  sep <- arm_separation(X, arm, B, A, ok)                    # signed, minor - major
  set.seed(1)
  nul <- mean(replicate(N_NULL, {
           r <- sample(rep(1:2, c(length(A), length(B))))
           abs(arm_separation(X, arm, cn[r == 2], cn[r == 1], ok)) }))
  ok_size  <- length(B) >= MIN_CLONE(n)
  accepted <- ok_size && abs(sep) >= ARM_MIN

  cl <- if (accepted) setNames(rep(c("c1", "c2"), c(length(A), length(B))), c(A, B))[cn]
        else          setNames(rep("c1", n), cn)
  attr(cl, "stats") <- data.frame(sample = S, n_cells = n,
      n_clones = if (accepted) 2L else 1L, n_major = length(A), n_minor = length(B),
      arm = names(sep), arm_sep = unname(sep), arm_sep_null = nul,
      size_ok = ok_size, accepted = accepted, stringsAsFactors = FALSE)
  cl
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
