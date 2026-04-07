# Mesothelioma scATAC-seq & scRNA-seq Analysis

Computational analysis accompanying the pleural mesothelioma single-cell chromatin accessibility and transcriptomics study. The pipeline covers joint scATAC/scRNA-seq processing, cell-type-specific subclustering, transcription factor (TF) footprinting, regulatory hub detection, and tumor biology analyses.

## Data Availability

Fragment files and Seurat objects are available on GEO under series **GSE311510**.

## Acknowledgements

Several utility scripts in `utils/` were adapted from the [GreenleafLab/HDMA](https://github.com/GreenleafLab/HDMA) repository.

---

## Repository Structure

```
git_repo/
├── main_analysis/                  # Whole-cohort scATAC-seq analysis
│   ├── scatac_main_ArchR.R         # ArchR project setup, QC, clustering, peak calling
│   ├── scatac_main_hub_analysis.R  # Regulatory hub detection across all cell types
│   ├── scatac_main_chrombpnet.R    # ChromBPNet TF footprinting [UNDER CONSTRUCTION]
│   └── export_fragments_for_GEO.R # Fragment file export for GEO submission
│
├── tumor_analysis/                 # Malignant cell compartment
│   ├── scatac_tumor_SOX9.R         # SOX9-driven chromatin programs
│   ├── scatac_tumor_SOX9_chromBPnet.R  # [UNDER CONSTRUCTION]
│   ├── scatac_tumor_P23.R          # Patient P23-specific differential peaks
│   ├── discover_oncogenic_TFs.R    # chromVAR + gene score TF ranking
│   ├── enrichment_cnmfs.R          # Pathway enrichment of cNMF tumor modules
│   ├── epigenomic_features_correlated_to_scS_score.R  # Sarcomatoid score correlations
│   ├── bulkRNA_analysis.R          # Bulk RNA-seq survival and subtype analysis
│   ├── CCLE_cell_lines_RNAseq.R    # DepMap/CCLE cell line RNA + CNV analysis
│   └── CRISPR_cropseq.R            # CRISPRi CROPseq analysis
│
├── tnk_analysis/                   # T cell / NK cell compartment
│   ├── scatac_NKTcells_ArchR.R     # NKT scATAC subclustering
│   ├── scatac_NKT_cells_hub_analysis.R  # Hub analysis in NKT cells
│   ├── scatac_NKT_cells_chromBPnet.R    # [UNDER CONSTRUCTION]
│   ├── scrna_tnk.R                 # T/NK scRNA-seq analysis + hdWGCNA
│   └── HiChIP_CD8_Ex_NR4A2_enhancer.R  # HiChIP validation of NR4A2 enhancer
│
├── myeloid_analysis/               # Myeloid compartment
│   ├── scatac_myeloid_ArchR.R      # Myeloid scATAC subclustering
│   ├── scatac_myeloid_hub_analysis.R  # Hub analysis in myeloid cells
│   └── scatac_myeloid_chromBPnet.R    # [UNDER CONSTRUCTION]
│
├── stroma_analysis/                # Stromal compartment (fibroblasts, endothelial, smooth muscle)
│   ├── scatac_stroma_ArchR.R       # Stroma scATAC subclustering
│   ├── scatac_stroma_hub_analysis.R
│   ├── scatac_stroma_lineage_analysis.R
│   ├── scatac_endothelial_ArchR.R
│   ├── scatac_endothelial_hub_analysis.R
│   ├── scatac_endothelial_chromBPnet.R  # [UNDER CONSTRUCTION]
│   ├── scrna_stroma.R
│   └── scrna_endothelial.R
│
├── B_cells_analysis/               # B cell / plasma cell compartment
│   ├── scatac_Bcells_ArchR.R
│   └── scrna_BCells.R
│
├── utils/                          # Shared utility scripts
│   ├── load_packages.R             # Load all standard R packages
│   ├── useful_functions.R          # General helper functions
│   ├── scATAC_functions.R          # ArchR/scATAC helper functions
│   ├── ggplot_aestetics.R          # ggplot2 themes
│   ├── palettes.R                  # Color palettes
│   ├── DAG.R                       # Differential accessibility (DAG wrapper)
│   ├── DEG_standard.R              # Differential expression wrapper
│   ├── chromVAR.R                  # chromVAR deviation scoring
│   ├── activeTFs.R                 # Active TF identification
│   ├── fGSEA_enrichment.R          # Gene set enrichment wrapper
│   ├── Hubs_finder.R               # Regulatory hub detection algorithm
│   ├── hubs_track.R                # Hub visualization tracks
│   ├── knnGen.R                    # KNN graph generation
│   ├── addCoax.R                   # Co-accessibility utilities
│   ├── callPeaks.R                 # Peak calling helpers
│   ├── DEG_standard.R              # Differential expression wrapper
│   ├── HDMA_pipeline.py            # HDMA pipeline
│   ├── HIC_processing.sh / HiC_Pro.sh  # HiC data processing
│   ├── chromBPnet_*.sh / *.py / *.R    # ChromBPNet pipeline [UNDER CONSTRUCTION]
│   ├── TFmodisco_*.sh              # TF-MoDISco motif discovery [UNDER CONSTRUCTION]
│   └── finemo_*.sh                 # FiNeMo motif calling [UNDER CONSTRUCTION]
│
├── files/                          # Reference files and precomputed objects
│   ├── barcode_annotation.csv      # Cell barcode → cell type mapping (all cells)
│   ├── normal_lung_mesothelium_barcodes.csv     # Normal mesothelium barcodes (combined)
│   ├── normal_lung_mesothelium_scrna_barcodes.csv  # Normal mesothelium scRNA barcodes
│   ├── normal_lung_mesothelium_scatac_barcodes.csv # Normal mesothelium scATAC barcodes
│   ├── normal_lung_peaks.rds       # Peak set from normal lung reference
│   ├── SCENIC_regulons.rds         # pySCENIC output regulons
│   ├── PCA_regulons.rds            # PCA of regulon activity
│   ├── scRNA_PM_modules.rds        # cNMF gene modules (scRNA)
│   ├── cnmf_myeloid_per_sample.rds # Myeloid cNMF modules
│   ├── malignant_cnmf_genelist_25_nfeat_5000.rds  # Tumor cNMF gene lists
│   ├── SE_regions_SE_database.rds  # Super-enhancer database regions
│   ├── Yang_Sox9_genelist.rds      # SOX9 target gene list (Yang et al.)
│   ├── Blum_et_al_SE_score.csv     # Super-enhancer scores from Blum et al.
│   ├── Riegel_CD8ex_DAP.csv        # CD8 exhaustion DAPs from Riegel et al.
│   ├── Tcell_exhaustion_genes_PMID37091230.csv
│   ├── hg38-blacklist.v2.bed       # ENCODE hg38 blacklist regions
│   ├── h.all.v7.4.symbols.gmt      # MSigDB Hallmark gene sets
│   ├── c5.bp.v7.1.symbol.gmt       # MSigDB GO biological process gene sets
│   └── ENCODE_query.txt            # ENCODE API query used for data retrieval
│
├── environment.yaml                # Conda environment spec (meso_scatac)
├── software_used.csv               # Full software inventory with versions and environments
└── readme.md                       # This file
```

---

## Setup

### 1. Create the main conda environment

```bash
conda env create -f environment.yaml
conda activate meso_scatac
```

Then install GitHub-only R packages from within R:

```r
# ArchR (primary scATAC framework)
devtools::install_github("GreenleafLab/ArchR", ref="master", repos=BiocManager::repositories())

# chromVAR motif set
devtools::install_github("GreenleafLab/chromVARmotifs")

# scATAC utility functions
devtools::install_github("crazyhottommy/scATACutils")

# Copy number from scATAC
devtools::install_github("colomemaria/epiAneufinder")

# Fast Wilcoxon for marker detection
devtools::install_github("immunogenomics/presto")

# hdWGCNA (co-expression networks from single-cell data)
devtools::install_github("smorabit/hdWGCNA", ref="dev")

# Additional Bioconductor packages
BiocManager::install(c("JASPAR2020", "EnhancedVolcano", "karyoploteR",
                       "DropletUtils", "SCENT", "scCustomize"))
```

### 2. Additional conda environments

| Environment   | Purpose                                                        | Python |
|---------------|----------------------------------------------------------------|--------|
| `chrombpnet`  | ChromBPNet deep learning footprinting (**under construction**) | 3.8    |
| `finemo`      | FiNeMo motif calling (**under construction**)                  | 3.10   |

> cNMF and pySCENIC outputs are provided as precomputed `.rds` files in `files/` and do not need to be re-run.

See `software_used.csv` for full version details of every package in every environment.

---

## Running the analysis

All analysis scripts expect to be run from the project directory (one level above `git_repo/`). Each script sources shared utilities from `git_repo/utils/` using relative paths.

### scATAC-seq (ArchR)

```r
conda activate meso_scatac
R
# Main cohort
source('git_repo/main_analysis/scatac_main_ArchR.R')

# Cell-type-specific subclustering (run independently)
source('git_repo/tumor_analysis/scatac_tumor_SOX9.R')
source('git_repo/tnk_analysis/scatac_NKTcells_ArchR.R')
source('git_repo/myeloid_analysis/scatac_myeloid_ArchR.R')
source('git_repo/stroma_analysis/scatac_stroma_ArchR.R')
source('git_repo/B_cells_analysis/scatac_Bcells_ArchR.R')
```

### scRNA-seq

```r
conda activate meso_scatac
R
source('git_repo/tnk_analysis/scrna_tnk.R')
source('git_repo/stroma_analysis/scrna_stroma.R')
source('git_repo/B_cells_analysis/scrna_BCells.R')
```

---

## Key R packages

The table below lists the main R packages. For complete version pinning see `software_used.csv`; for the conda spec see `environment.yaml`.

| Category              | Package(s)                                                                 |
|-----------------------|----------------------------------------------------------------------------|
| scATAC-seq            | ArchR, Signac, chromVAR, chromVARmotifs, epiAneufinder                    |
| Genome annotation     | BSgenome.Hsapiens.UCSC.hg38/hg19, EnsDb.Hsapiens.v86, TxDb.*.hg38/hg19, TFBSTools, JASPAR2020 |
| Genomic ranges        | GenomicRanges, GenomicFeatures, regioneR, rtracklayer, rGREAT              |
| scRNA-seq             | Seurat, Signac, scran, scuttle, SingleCellExperiment, sctransform          |
| Batch correction      | harmony                                                                    |
| Differential analysis | limma, edgeR, DESeq2, presto, fgsea, clusterProfiler, DOSE, topGO, GSA    |
| Gene signatures       | UCell, SCENT, hdWGCNA, WGCNA, destiny                                     |
| CNV inference         | inferCNV, epiAneufinder                                                    |
| Survival              | survival, survminer, RTCGA                                                 |
| Visualization         | ggplot2, ComplexHeatmap, patchwork, cowplot, ggpubr, ggrepel, ggridges, circlize, viridis, RColorBrewer, paletteer, ggtree, corrplot, EnhancedVolcano, scCustomize |
| Annotation            | org.Hs.eg.db, org.Mm.eg.db, biomaRt                                       |

---

## Notes

- Scripts with `[UNDER CONSTRUCTION]` in the tree above (ChromBPNet, TF-MoDISco, FiNeMo pipeline) are incomplete and should not be run as-is.
- HPC job submission uses LSF (`bsub`). All batch scripts assume the Minerva cluster.
