###############################################################
# R2_Q14 P8 probe: is P8 subclonal for chr3p/BAP1, and can we split it in
# both scRNA (inferCNV) and scATAC (gene score)?  Reconnaissance only.
###############################################################
suppressMessages({ library(ArchR); library(Matrix) })
addArchRThreads(threads = 1); addArchRGenome("hg38")

BASE   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment"
ICNV   <- file.path(BASE, "scrna", "infercnv")
OUTDIR <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R2_Q14"
setwd(OUTDIR)

# BAP1 hg38: chr3:52,401,003-52,410,357  (chr3p21.1)
BAP1_CHR <- "chr3"; BAP1_S <- 52401003; BAP1_E <- 52410357
CHR3_CENT <- 90e6   # ~chr3 centromere; p-arm = start < this

## ---------------- (1) inferCNV : P8 chr3p / BAP1 CNV ----------------
cat("=== loading inferCNV object ===\n")
obj <- readRDS(file.path(ICNV, "run.final.infercnv_obj"))
cat("class:", class(obj), "\nslots:", paste(slotNames(obj), collapse=", "), "\n")
ed <- obj@expr.data
go <- obj@gene_order
cat("expr.data dim:", paste(dim(ed), collapse=" x "), "\n")
cat("gene_order cols:", paste(colnames(go), collapse=", "), "\n")
print(head(go, 3))
# normalize chr labels
chrcol <- go[[1]]; if (!grepl("^chr", chrcol[1])) chrcol <- paste0("chr", chrcol)
startcol <- go[[2]]; stopcol <- go[[3]]
cat("BAP1 in gene_order:", "BAP1" %in% rownames(go), "\n")
if ("BAP1" %in% rownames(go)) print(go["BAP1", ])

# P8 cells from annot_df
annot <- read.table(file.path(ICNV, "annot_df.txt"), header=FALSE, sep="\t",
                    stringsAsFactors=FALSE, quote="\"")
colnames(annot) <- c("cell","sample")
p8_cells <- annot$cell[annot$sample == "P8"]
p8_cells <- intersect(p8_cells, colnames(ed))
cat("P8 cells in inferCNV:", length(p8_cells), "\n")

chr3p <- which(chrcol == "chr3" & startcol < CHR3_CENT)
bap1_win <- which(chrcol == "chr3" & startcol >= (BAP1_S - 5e6) & stopcol <= (BAP1_E + 5e6))
cat("chr3p genes:", length(chr3p), " | BAP1 +/-5Mb window genes:", length(bap1_win), "\n")

p8ed <- ed[, p8_cells, drop=FALSE]
chr3p_score <- colMeans(p8ed[chr3p, , drop=FALSE])
bap1_score  <- colMeans(p8ed[bap1_win, , drop=FALSE])
bap1_gene   <- if ("BAP1" %in% rownames(p8ed)) p8ed["BAP1", ] else rep(NA, length(p8_cells))

# subclusters
og <- tryCatch(read.table(file.path(ICNV, "infercnv_subclusters.observation_groupings.txt"),
               header=TRUE, sep=" ", stringsAsFactors=FALSE, row.names=1, check.names=FALSE),
               error=function(e) NULL)
subcl <- if (!is.null(og)) og[p8_cells, 1] else NA

p8tab <- data.frame(cell=p8_cells, subcluster=subcl,
                    chr3p_cnv=chr3p_score, bap1win_cnv=bap1_score, BAP1_cnv=bap1_gene)
write.csv(p8tab, file.path(OUTDIR, "P8_scRNA_infercnv_probe.csv"), row.names=FALSE)

cat("\n--- P8 chr3p CNV score distribution (1 = neutral, <1 = loss) ---\n")
print(round(quantile(chr3p_score, c(0,.05,.1,.25,.5,.75,.9,.95,1)), 4))
cat("BAP1-window CNV quantiles:\n"); print(round(quantile(bap1_score, c(0,.1,.25,.5,.75,.9,1)),4))
cat("\nmean chr3p CNV by subcluster (top/bottom):\n")
sm <- sort(tapply(chr3p_score, subcl, mean))
print(round(head(sm,8),4)); cat("...\n"); print(round(tail(sm,8),4))
cat("n subclusters in P8:", length(unique(subcl)), "\n")

## ---------------- (2) scATAC : P8 BAP1 gene score + chr3p access ----------------
cat("\n=== loading ArchR (tumor_compartment/scatac_ArchR) ===\n")
archp <- loadArchRProject(file.path(BASE, "scatac_ArchR"), showLogo=FALSE)
cd <- getCellColData(archp)
cat("scATAC samples:\n"); print(table(cd$Sample))
p8_atac <- rownames(cd)[cd$Sample == "P8"]
cat("P8 scATAC cells:", length(p8_atac), "\n")

gsm <- getMatrixFromProject(archp, useMatrix="GeneScoreMatrix")
grn <- rowData(gsm)
gs  <- assay(gsm)
rownames(gs) <- grn$name
gs_p8 <- gs[, colnames(gs) %in% p8_atac, drop=FALSE]
cat("BAP1 in GeneScoreMatrix:", "BAP1" %in% rownames(gs), "\n")
bap1_gs <- if ("BAP1" %in% rownames(gs)) gs_p8["BAP1", ] else rep(NA, ncol(gs_p8))
# chr3p accessibility score: mean gene score of chr3p genes
chr3p_genes_atac <- grn$name[as.character(seqnames(grn)) == "chr3" & start(grn) < CHR3_CENT]
chr3p_genes_atac <- intersect(chr3p_genes_atac, rownames(gs_p8))
chr3p_atac <- colMeans(gs_p8[chr3p_genes_atac, , drop=FALSE])
cat("chr3p genes in scATAC gene score:", length(chr3p_genes_atac), "\n")

atacTab <- data.frame(cell=colnames(gs_p8), BAP1_genescore=bap1_gs, chr3p_access=chr3p_atac)
write.csv(atacTab, file.path(OUTDIR, "P8_scATAC_genescore_probe.csv"), row.names=FALSE)
cat("\n--- P8 scATAC BAP1 gene score quantiles ---\n")
print(round(quantile(bap1_gs, c(0,.1,.25,.5,.75,.9,1)),4))
cat("BAP1 gs vs chr3p access spearman rho:",
    round(cor(bap1_gs, chr3p_atac, method="spearman"),3), "\n")
cat("\nDONE\n")
