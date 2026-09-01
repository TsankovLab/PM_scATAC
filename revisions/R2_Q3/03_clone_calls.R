###############################################################################
# STEP 3 -- turn the 9 epiAneufinder result tables into ONE set of clone definitions.
#
# This is the only script that decides what a clone is; every figure and every
# statistic downstream reads its output, so the clone labels cannot drift.
#
#   per sample : cut epiAneufinder's own dendrogram at k = 2 (epi_clone_split in
#                00_common.R), keep clusters of >= max(20, 5%) cells
#   per clone  : a CNV profile over the shared 5 Mb bins,
#                profile = (fraction of cells called GAIN) - (fraction called LOSS)
#   per sample : the between-clone difference (clone2 - clone1) and the arm carrying
#                the largest one -- this is what says whether a split is one focal
#                event (P4 chr8q) or diffuse noise
#
# Outputs
#   epi_clone_labels.csv     cell -> clone, cells named "<S>#<barcode>" to match ArchR
#   epi_clone_profiles.rds   Z (leaf x bin), DL (sample x bin delta), meta, coord
#   epi_subclone_summary.csv one row per tumour: clone sizes, top arm, delta
#   epi_leaves.csv           one row per clone: size, chr8q gain, driver arm
###############################################################################
suppressMessages({ library(data.table) })
source("00_common.R")

prof <- list(); meta <- list(); delta <- list(); lab <- list(); COORD <- NULL
for (S in SAMPLES){
  e  <- read_epi(S); M <- e$M; d <- e$d
  cl <- epi_clone_split(S)
  gs <- sort(unique(cl))
  cat(sprintf("%-5s %d clones (%s cells of %d)\n", S, length(gs),
              paste(as.integer(table(cl)[gs]), collapse = "/"), ncol(M)))
  lab[[S]] <- data.table(sample = S, cell = paste0(S, "#", names(cl)), clone = unname(cl))

  is8q <- d$seq == "chr8" & d$start > CEN8
  pl <- list()
  for (g in gs){
    cells <- names(cl)[cl == g]
    p <- rowMeans(M[, cells, drop = FALSE] == 2) - rowMeans(M[, cells, drop = FALSE] == 0)
    names(p) <- e$bin
    leaf <- sprintf("%s %s", S, g)
    prof[[leaf]] <- p; pl[[g]] <- p
    meta[[leaf]] <- data.table(leaf = leaf, sample = S, clone = g, n_cells = length(cells),
                               chr8q_gain = mean(M[is8q, cells, drop = FALSE] == 2))
  }
  ## every sample splits into exactly two retained clones; guard in case that changes
  stopifnot(identical(gs, c("c1", "c2")))
  delta[[S]] <- pl[["c2"]] - pl[["c1"]]
  if (is.null(COORD)) COORD <- data.table(bin = e$bin, chr = d$seq, mb = d$start / 1e6)
}
LAB <- rbindlist(lab); fwrite(LAB, "epi_clone_labels.csv")
cat("\nlabelled cells:", nrow(LAB), "\n")

meta <- rbindlist(meta); setDF(meta, rownames = meta$leaf)
bins <- Reduce(intersect, lapply(prof, names))          # identical grid across samples
Z  <- t(sapply(prof[meta$leaf], function(v) v[bins])); colnames(Z) <- bins
DL <- t(sapply(delta[SAMPLES],  function(v) v[bins])); colnames(DL) <- bins
co <- COORD[match(bins, COORD$bin)]
arm <- paste0(sub("chr", "", co$chr), ifelse(co$mb < CENmb[co$chr], "p", "q"))
cat("leaves:", nrow(Z), " shared", WINDOW/1e6, "Mb bins:", ncol(Z), "\n")

## ---- driver arm per clone ---------------------------------------------------
## largest departure of this clone's profile from its sibling's
meta$driver <- NA_character_
for (l in meta$leaf){
  sib <- meta$leaf[meta$sample == meta[l, "sample"] & meta$leaf != l]
  ref <- if (length(sib)) colMeans(Z[sib, , drop = FALSE]) else colMeans(Z)
  dz  <- Z[l, ] - ref; b <- names(which.max(abs(dz)))
  meta[l, "driver"] <- sprintf("%s %s", armof(sub(":.*", "", b),
                                              as.numeric(sub(".*:", "", b)) / 1e6),
                               ifelse(dz[b] < 0, "loss", "gain"))
}
## the P4 clone carrying the chr8q amplification -- the one the Visium data validate
p4 <- meta[meta$sample == "P4", ]
CHR8Q <- p4$leaf[which.max(p4$chr8q_gain)]
meta$class <- ifelse(meta$leaf == CHR8Q, "P4 chr8q-amplified clone", "subclone")
cat(sprintf("P4 chr8q clone: %s (%d cells, chr8q gain %.3f vs %.3f in its sibling)\n",
            CHR8Q, meta[CHR8Q, "n_cells"], meta[CHR8Q, "chr8q_gain"],
            min(p4$chr8q_gain)))
write.csv(meta, "epi_leaves.csv", row.names = FALSE)

## ---- what separates the two clones of each tumour ---------------------------
amax <- t(sapply(SAMPLES, function(S){
  v <- DL[S, ]; i <- which.max(abs(v))
  c(arm = unname(arm[i]), delta = unname(round(v[i], 3))) }))
strength <- setNames(apply(abs(DL), 1, max), SAMPLES)
SUMM <- data.frame(sample = SAMPLES,
                   clone1_n = meta$n_cells[match(paste(SAMPLES, "c1"), meta$leaf)],
                   clone2_n = meta$n_cells[match(paste(SAMPLES, "c2"), meta$leaf)],
                   top_arm  = amax[, "arm"], delta = as.numeric(amax[, "delta"]),
                   max_abs_delta = round(strength, 3), row.names = NULL)
write.csv(SUMM, "epi_subclone_summary.csv", row.names = FALSE)
cat("\n=== largest between-clone difference per tumour ===\n"); print(SUMM, row.names = FALSE)

saveRDS(list(Z = Z, DL = DL, meta = meta, coord = co, arm = arm,
             amax = amax, strength = strength, chr8q_leaf = CHR8Q),
        "epi_clone_profiles.rds")
cat("\nDONE -> epi_clone_profiles.rds\n")
