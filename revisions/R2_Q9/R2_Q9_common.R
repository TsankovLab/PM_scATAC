# R2_Q9 -- shared setup: packages, aesthetics of git_repo/utils, and the improved sample palette.
# Every R2_Q9 script regenerates figures of git_repo that use palette_sample; the plotting code is copied from
# the original scripts and only the palette changes. Projects are loaded read-only (never saved) and all output
# goes to R2_Q9/Plots.
source(local({a <- grep("^--file=", commandArgs(FALSE), value = TRUE); for (d in c(if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))), getwd())) { while (!file.exists(file.path(d, "config.R")) && d != dirname(d)) d <- dirname(d); if (file.exists(file.path(d, "config.R"))) break }; Sys.getenv("MESO_CONFIG", file.path(d, "config.R"))}))   # data / results / reference roots: see config.R at the repo root
suppressPackageStartupMessages({
  library(ArchR); library(SummarizedExperiment); library(ComplexHeatmap); library(circlize)
  library(ggplot2); library(patchwork); library(paletteer); library(dplyr); library(RColorBrewer); library(viridis); library(scales)
})
root  = DATA_DIR
utils = file.path (root, 'git_repo', 'utils')
source (file.path (utils, 'ggplot_aestetics.R'))
source (file.path (utils, 'palettes.R'))                 # other palettes; palette_sample is replaced below
palette_sample_git = palette_sample
source (file.path (root, 'git_repo_claude', 'R2_Q9', 'R2_Q9_palette_sample.R'))
outdir = file.path (root, 'git_repo_claude', 'R2_Q9', 'Plots')
dir.create (outdir, showWarnings = FALSE)

addArchRThreads (threads = 1)
addArchRGenome ('hg38')
load_project = function (dir) suppressMessages (loadArchRProject (file.path (root, dir), showLogo = FALSE))
