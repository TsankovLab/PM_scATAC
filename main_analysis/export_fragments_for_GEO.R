### Export fragment files of annotated cells for GEO submission ####
# Get fragments from each sample
# library (rtracklayer)
# fragments_l = getFragmentsFromProject (archp)
# tmpdir = '$SCRATCH_DIR/meso_scatac_fragments'
# for (sam in names (fragments_l))
#   {
#     if (!file.exists(file.path(tmpdir,paste0(sam,'_fragments.bed'))))
#       export(fragments_l[[sam]], file.path(tmpdir, paste0(sam,'_fragments.bed')))
#   }
source(local({a <- grep("^--file=", commandArgs(FALSE), value = TRUE); for (d in c(if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))), getwd())) { while (!file.exists(file.path(d, "config.R")) && d != dirname(d)) d <- dirname(d); if (file.exists(file.path(d, "config.R"))) break }; Sys.getenv("MESO_CONFIG", file.path(d, "config.R"))}))   # data / results / reference roots: see config.R at the repo root
library(ArchR)
library(GenomicRanges)
library(rtracklayer)

fragments_l <- getFragmentsFromProject(archp)

# Export normal lung meso
fragments_l = fragments_l[c('RPL_280_neg_1','RPL_280_neg_2','RPL_Epi_1','RPL_Epi_2')]

tmpdir <- paste0(SCRATCH_DIR, "/meso_scatac_fragments")
dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)

# adjust if you want
sort_mem <- "20G"
n_threads <- 4
sort_tmp <- SCRATCH_DIR
force = T
for (sam in names(fragments_l) ) {
  
  bedfile <- file.path(tmpdir, paste0(sam, "_fragments.bed"))
  gzfile  <- file.path(tmpdir, paste0(sam, "_fragments.sorted.bed.gz"))
  
  if (!file.exists(gzfile)| force) {
    
    gr <- fragments_l[[sam]]
    
    # keep only barcode (remove "P1#")
    mcols(gr)$name <- sub("^.*#", "", as.character(mcols(gr)$RG))
    
    # required BED fields
    mcols(gr)$score <- 1
    strand(gr) <- "*"
    
    # export unsorted BED
    export(gr, bedfile, format = "BED")
    
    # sort -> bgzip
    cmd <- sprintf(
      "export LC_ALL=C; sort -T %s -S %s -k1,1V -k2,2n %s | bgzip -@ %d -c > %s",
      shQuote(sort_tmp),
      shQuote(sort_mem),
      shQuote(bedfile),
      n_threads,
      shQuote(gzfile)
    )
    
    status <- system(cmd)
    if (status != 0) {
      stop(sprintf("Sorting/bgzip failed for sample %s", sam))
    }
    
    # tabix index
    status <- system2("tabix", args = c("-f", "-p", "bed", gzfile))
    if (status != 0) {
      stop(sprintf("tabix indexing failed for sample %s", sam))
    }
  }
}

# # Export barcode annotation and fragments for GEO submission ####
# if (!exists ('fragments_l')) fragments_l = getFragmentsFromProject (archp)
# library(rtracklayer)

# # Example GRanges object
# # Export to a BED file
# lapply (names(fragments_l), function(x) export (fragments_l[[x]], paste0(x,"_fragments.bed")))




