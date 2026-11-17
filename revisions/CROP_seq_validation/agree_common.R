# Shared helpers for the bulk KO analyses (spec118-121, spec195, final6KO_normalise.R). Sourced after repro_common.R and profiles.R.
source(local({a <- grep("^--file=", commandArgs(FALSE), value = TRUE); for (d in c(if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))), getwd())) { while (!file.exists(file.path(d, "config.R")) && d != dirname(d)) d <- dirname(d); if (file.exists(file.path(d, "config.R"))) break }; Sys.getenv("MESO_CONFIG", file.path(d, "config.R"))}))   # data / results / reference roots: see config.R at the repo root
suppressPackageStartupMessages({ library(edgeR); library(limma) })
KOS <- ORDER; CM_RDS <- paste0(REPO_DIR, "/files/scRNA_PM_modules.rds")
cm_top50 <- function() { m <- readRDS(CM_RDS); m <- lapply(m, head, 50); m[order(as.integer(sub("Cm", "", names(m))))] }

# normalised expression of one line: list(y = genes x samples, d = DGEList (TMM) for voom). counts = symbol-level counts of the line's samples, already gene-filtered
normalise <- function(counts, norm) { d <- calcNormFactors(DGEList(counts))
  y <- switch(norm, cp10k = log2(sweep(counts, 2, colSums(counts), "/") * 1e4 + 1), tmm = cpm(d, log = TRUE, prior.count = 2), voomw = cpm(d, log = TRUE, prior.count = 2),
              vst = { suppressPackageStartupMessages(library(DESeq2)); x <- round(counts); storage.mode(x) <- "integer"; SummarizedExperiment::assay(vst(DESeqDataSetFromMatrix(x, data.frame(a = rep(1, ncol(x))), ~ 1), blind = TRUE)) })
  list(y = y, d = d) }
