###############################################################################
# STEP 3 -- turn the 9 epiAneufinder result tables into ONE set of clone definitions.
#
# This is the only script that decides what a clone is; every figure and every
# statistic downstream reads its output, so the clone labels cannot drift.
#
#   per sample : take the primary branch point of the cell dendrogram and accept it as a
#                clone boundary only if the two branches differ by >= ARM_MIN on some
#                chromosome arm and both clear the size floor (epi_clone_split in
#                00_common.R).  A tumour that fails the test is ONE clone -- no cell is
#                discarded either way.  ARM_MIN is calibrated on P4 (see 00_common.R).
#   per clone  : a CNV profile over the shared 5 Mb bins,
#                profile = (fraction of cells called GAIN) - (fraction called LOSS)
#   per sample : for the tumours that split, the between-clone difference (minor clone
#                c2 minus major clone c1) and the arm carrying the largest one
#
# Outputs
#   epi_clone_labels.csv     cell -> clone, cells named "<S>#<barcode>" to match ArchR
#   epi_clone_decision.csv   one row per tumour: the split test and its result
#   epi_clone_profiles.rds   Z (leaf x bin), DL (split sample x bin delta), meta, coord
#   epi_subclone_summary.csv one row per tumour: clone sizes, separation, top arm
#   epi_leaves.csv           one row per clone: size, chr8q gain, driver arm
###############################################################################
suppressMessages({ library(data.table) })
source("00_common.R")

prof <- list(); meta <- list(); delta <- list(); lab <- list(); dec <- list(); COORD <- NULL
for (S in SAMPLES){
  e  <- read_epi(S); M <- e$M; d <- e$d
  cl <- epi_clone_split(S)
  st <- attr(cl, "stats"); dec[[S]] <- st
  gs <- sort(unique(cl))
  cat(sprintf("%-5s %5d cells | %s %+.3f (noise floor %.3f) -> %s\n", S, ncol(M),
              st$arm, st$arm_sep, st$arm_sep_null,
              if (st$accepted) sprintf("2 clones, %d / %d", st$n_major, st$n_minor)
              else "1 clone (below threshold)"))
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
  if (length(gs) == 2) delta[[S]] <- pl[["c2"]] - pl[["c1"]]
  if (is.null(COORD)) COORD <- data.table(bin = e$bin, chr = d$seq, mb = d$start / 1e6)
}
LAB <- rbindlist(lab); fwrite(LAB, "epi_clone_labels.csv")
SPLIT <- names(delta)                                   # tumours that passed the test
cat(sprintf("\nlabelled cells: %d | %d of %d tumours split (%s)\n",
            nrow(LAB), length(SPLIT), length(SAMPLES), paste(SPLIT, collapse = ", ")))

## ---- the split test, tumour by tumour ---------------------------------------
DEC <- rbindlist(dec)
write.csv(DEC, "epi_clone_decision.csv", row.names = FALSE)
cat("\n=== primary split test (ARM_MIN =", ARM_MIN, ") ===\n")
print(transform(as.data.frame(DEC), arm_sep = round(arm_sep, 3),
                arm_sep_null = round(arm_sep_null, 3)), row.names = FALSE)

meta <- rbindlist(meta); setDF(meta, rownames = meta$leaf)
bins <- Reduce(intersect, lapply(prof, names))          # identical grid across samples
Z  <- t(sapply(prof[meta$leaf], function(v) v[bins])); colnames(Z) <- bins
DL <- t(sapply(delta[SPLIT],   function(v) v[bins])); colnames(DL) <- bins
co <- COORD[match(bins, COORD$bin)]
arm <- paste0(sub("chr", "", co$chr), ifelse(co$mb < CENmb[co$chr], "p", "q"))
cat("leaves:", nrow(Z), " shared", WINDOW/1e6, "Mb bins:", ncol(Z), "\n")

## ---- driver arm per clone ---------------------------------------------------
## largest departure of this clone's profile from its sibling's; a tumour that was not
## split has no sibling and no driver
meta$driver <- "unsplit"
for (l in meta$leaf[meta$sample %in% SPLIT]){
  sib <- meta$leaf[meta$sample == meta[l, "sample"] & meta$leaf != l]
  dz  <- Z[l, ] - Z[sib, ]; b <- names(which.max(abs(dz)))
  meta[l, "driver"] <- sprintf("%s %s", armof(sub(":.*", "", b),
                                              as.numeric(sub(".*:", "", b)) / 1e6),
                               ifelse(dz[b] < 0, "loss", "gain"))
}
## the P4 clone carrying the chr8q amplification -- the one the Visium data validate
p4 <- meta[meta$sample == "P4", ]
stopifnot(nrow(p4) == 2)          # the calibration sample must split; see 00_common.R
CHR8Q <- p4$leaf[which.max(p4$chr8q_gain)]
meta$class <- ifelse(meta$leaf == CHR8Q, "P4 chr8q-amplified clone",
              ifelse(meta$driver == "unsplit", "whole tumour (not split)", "subclone"))
cat(sprintf("P4 chr8q clone: %s (%d cells, chr8q gain %.3f vs %.3f in its sibling)\n",
            CHR8Q, meta[CHR8Q, "n_cells"], meta[CHR8Q, "chr8q_gain"],
            min(p4$chr8q_gain)))
write.csv(meta, "epi_leaves.csv", row.names = FALSE)

## ---- what separates the two clones of each tumour that split ----------------
amax <- t(sapply(SPLIT, function(S){
  v <- DL[S, ]; i <- which.max(abs(v))
  c(arm = unname(arm[i]), delta = unname(round(v[i], 3))) }))
strength <- setNames(apply(abs(DL), 1, max), SPLIT)
SUMM <- data.frame(sample = SAMPLES,
                   n_cells  = DEC$n_cells,
                   n_clones = DEC$n_clones,
                   clone1_n = meta$n_cells[match(paste(SAMPLES, "c1"), meta$leaf)],
                   clone2_n = meta$n_cells[match(paste(SAMPLES, "c2"), meta$leaf)],
                   sep_arm  = DEC$arm, sep = round(DEC$arm_sep, 3),
                   top_bin_arm   = amax[match(SAMPLES, SPLIT), "arm"],
                   top_bin_delta = as.numeric(amax[match(SAMPLES, SPLIT), "delta"]),
                   max_abs_delta = round(unname(strength[match(SAMPLES, SPLIT)]), 3),
                   row.names = NULL)
write.csv(SUMM, "epi_subclone_summary.csv", row.names = FALSE)
cat("\n=== clones per tumour ===\n"); print(SUMM, row.names = FALSE)

saveRDS(list(Z = Z, DL = DL, meta = meta, coord = co, arm = arm, split = SPLIT,
             decision = as.data.frame(DEC), amax = amax, strength = strength,
             chr8q_leaf = CHR8Q),
        "epi_clone_profiles.rds")
cat("\nDONE -> epi_clone_profiles.rds\n")
