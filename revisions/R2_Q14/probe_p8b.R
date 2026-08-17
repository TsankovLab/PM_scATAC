###############################################################
# R2_Q14 probe B: per-sample chr3p / BAP1 CNV bimodality (is ANY epithelioid
# BAP1-lost sample subclonal for chr3p, i.e. has a retained subclone?),
# normalized to each cell's genome-wide median. Plus fixed scATAC chr3p score.
###############################################################
suppressMessages({ library(ArchR); library(Matrix) })
addArchRThreads(threads = 1); addArchRGenome("hg38")

BASE   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment"
ICNV   <- file.path(BASE, "scrna", "infercnv")
OUTDIR <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R2_Q14"
setwd(OUTDIR)
CHR3_CENT <- 90e6; BAP1_S <- 52401003; BAP1_E <- 52410357

## ---- inferCNV per-sample chr3p analysis ----
obj <- readRDS(file.path(ICNV, "run.final.infercnv_obj"))
ed  <- obj@expr.data
go  <- obj@gene_order
chrcol <- go[[1]]; if (!grepl("^chr", chrcol[1])) chrcol <- paste0("chr", chrcol)
startcol <- go[[2]]; stopcol <- go[[3]]
chr3p    <- which(chrcol == "chr3" & startcol < CHR3_CENT)
bap1_win <- which(chrcol == "chr3" & startcol >= (BAP1_S-5e6) & stopcol <= (BAP1_E+5e6))

annot <- read.table(file.path(ICNV,"annot_df.txt"), header=FALSE, sep="\t", quote="\"",
                    stringsAsFactors=FALSE); colnames(annot) <- c("cell","sample")
annot <- annot[annot$cell %in% colnames(ed), ]

# per-cell genome median (neutral baseline) and normalized chr3p / bap1 window
gmed <- apply(ed, 2, median)
chr3p_raw <- colMeans(ed[chr3p, , drop=FALSE])
bap1_raw  <- colMeans(ed[bap1_win, , drop=FALSE])
chr3p_norm <- chr3p_raw / gmed
bap1_norm  <- bap1_raw  / gmed

samp <- setNames(annot$sample, annot$cell)
df <- data.frame(cell=colnames(ed), sample=samp[colnames(ed)],
                 chr3p_raw, chr3p_norm, bap1_raw, bap1_norm)
df <- df[!is.na(df$sample), ]
write.csv(df, file.path(OUTDIR,"allsample_scRNA_chr3p_cnv.csv"), row.names=FALSE)

cat("=== per-sample chr3p (normalized to genome median) ===\n")
tum <- setdiff(unique(df$sample), "reference")
summ <- do.call(rbind, lapply(sort(tum), function(s){
  v <- df$chr3p_norm[df$sample==s]; b <- df$bap1_norm[df$sample==s]
  data.frame(sample=s, n=length(v),
    chr3p_med=median(v), chr3p_iqr=IQR(v),
    frac_retained_chr3p=mean(v > 0.98), frac_lost_chr3p=mean(v < 0.92),
    bap1_med=median(b), frac_retained_bap1=mean(b > 0.98), frac_lost_bap1=mean(b < 0.90))
}))
ref_chr3p <- median(df$chr3p_norm[df$sample=="reference"])
cat("reference chr3p_norm median:", round(ref_chr3p,3), " (neutral ~1.0)\n\n")
print(summ, digits=3, row.names=FALSE)
write.csv(summ, file.path(OUTDIR,"persample_chr3p_bimodality.csv"), row.names=FALSE)

cat("\n>> Candidate subclonal samples = both frac_retained and frac_lost sizeable\n")

## ---- scATAC: BAP1 gene score + FIXED chr3p accessibility, all tumor samples ----
cat("\n=== scATAC gene score ===\n")
archp <- loadArchRProject(file.path(BASE,"scatac_ArchR"), showLogo=FALSE)
cd <- getCellColData(archp)
gsm <- getMatrixFromProject(archp, useMatrix="GeneScoreMatrix")
rd  <- rowData(gsm)
cat("rowData cols:", paste(colnames(rd), collapse=", "), "\n")
gs  <- assay(gsm); rownames(gs) <- rd$name
seqn <- as.character(rd$seqnames); st <- rd$start
chr3p_atac_genes <- rd$name[seqn=="chr3" & st < CHR3_CENT]
chr3p_atac_genes <- intersect(chr3p_atac_genes, rownames(gs))
cat("chr3p genes in scATAC:", length(chr3p_atac_genes), "\n")

cellsamp <- cd[colnames(gs), "Sample"]
tum_atac <- c("P1","P3","P4","P5","P8","P10","P11","P12","P13","P14","P23")
keep <- cellsamp %in% tum_atac
atac <- data.frame(cell=colnames(gs)[keep], sample=cellsamp[keep],
                   BAP1_gs = gs["BAP1", keep],
                   chr3p_access = colMeans(gs[chr3p_atac_genes, keep, drop=FALSE]))
write.csv(atac, file.path(OUTDIR,"allsample_scATAC_bap1_chr3p.csv"), row.names=FALSE)
cat("\nper-sample scATAC BAP1 gene score (median) & chr3p access:\n")
sa <- do.call(rbind, lapply(sort(unique(atac$sample)), function(s){
  d <- atac[atac$sample==s,]
  data.frame(sample=s, n=nrow(d), BAP1_gs_med=median(d$BAP1_gs),
             BAP1_gs_iqr=IQR(d$BAP1_gs), chr3p_access_med=median(d$chr3p_access))
}))
print(sa, digits=3, row.names=FALSE)
cat("\nDONE\n")
