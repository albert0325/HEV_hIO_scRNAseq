# =============================================================================
# 01_quality_control.R
# scRNA-seq quality control, sample integration, and dimensionality reduction
#
# Two 10x Genomics libraries are loaded and integrated via SCTransform:
#   - HEV_norm : cells exposed to HEV under standard conditions (a.k.a. HEV-1 in the paper)
#   - HEV_pp   : cells exposed to HEV under alternative passage conditions (a.k.a. HEV-2 in the paper)
#
# Infected cells are flagged by detectable ORF1 or ORF2 viral transcripts.
#
# Input:  data/HEV_norm/filtered_feature_bc_matrix/
#         data/HEV_pp/filtered_feature_bc_matrix/
# Output: outputs/S9A.tiff   — QC metrics projected onto UMAP
#         hev                — Seurat object carried forward to 02_annotation.R
#
# =============================================================================

source("R/utils.R")
library(Seurat)
library(tidyverse)
library(future)

out_dir <- "outputs/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Step 1: Load raw count matrices
# -----------------------------------------------------------------------------

hev.data  <- Read10X("data/HEV_norm/filtered_feature_bc_matrix/")
dat1      <- CreateSeuratObject(counts = hev.data,  project = "hev",   min.cells = 3, min.features = 200)

hevpp.data <- Read10X("data/HEV_pp/filtered_feature_bc_matrix/")
dat2       <- CreateSeuratObject(counts = hevpp.data, project = "hevpp", min.cells = 3, min.features = 200)

hev <- merge(dat1, y = dat2, add.cell.ids = c("norm", "pp"), project = "HEV")

# -----------------------------------------------------------------------------
# Step 2: Quality control
# Thresholds chosen by inspection of per-sample distributions.
# nFeature_RNA < 9,000 and nCount_RNA < 100,000 remove likely doublets.
# percent.mt  < 50.
# -----------------------------------------------------------------------------

hev[["percent.mt"]] <- PercentageFeatureSet(hev, pattern = "^MT-")

hev <- subset(
  hev,
  subset = nFeature_RNA > 200   &
           nFeature_RNA < 9000  &
           percent.mt   < 50    &
           nCount_RNA   < 100000
)

cat(sprintf("Cells retained after QC: %d\n", ncol(hev)))

# -----------------------------------------------------------------------------
# Step 3: Flag HEV-infected cells
# A cell is called infected if it carries ≥1 UMI from ORF1 or ORF2,
# the two principal HEV open reading frames present in the custom genome index.
# -----------------------------------------------------------------------------

inf.counts      <- FetchData(hev, vars = c("ORF1", "ORF2"), layer = "counts")
hev$infected    <- (inf.counts$ORF1 + inf.counts$ORF2) > 0

# Infection summary per library
cat("\nInfection status per library:\n")
print(table(hev$orig.ident, hev$infected))

hev$infected_status <- ifelse(hev$infected, "Positive", "Negative")

# -----------------------------------------------------------------------------
# Step 4: Integration via SCTransform
# Each library is normalised independently with SCTransform before
# anchor-based integration to remove batch effects between conditions.
# -----------------------------------------------------------------------------

options(future.globals.maxSize = Inf)

hev.list <- SplitObject(hev, split.by = "orig.ident")
hev.list <- lapply(hev.list, SCTransform, verbose = FALSE)

features <- SelectIntegrationFeatures(object.list = hev.list, nfeatures = 3000)
hev.list <- PrepSCTIntegration(object.list = hev.list, anchor.features = features)
anchors  <- FindIntegrationAnchors(
  object.list          = hev.list,
  normalization.method = "SCT",
  anchor.features      = features
)
hev <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# -----------------------------------------------------------------------------
# Step 5: Dimensionality reduction and clustering
# -----------------------------------------------------------------------------

hev <- RunPCA(hev, verbose = FALSE)
hev <- RunUMAP(hev, dims = 1:30)
hev <- FindNeighbors(hev, dims = 1:30)
hev <- FindClusters(hev, resolution = 0.8)

# -----------------------------------------------------------------------------
# Step 6: Visualise QC metrics on UMAP embedding (Suppl. Fig. S9A)
# Projecting QC metrics onto the UMAP confirms that low-quality cells
# do not form spatially coherent clusters post-filtering.
# -----------------------------------------------------------------------------

p <- FeaturePlot(
  hev,
  features  = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  reduction = "umap",
  ncol      = 3
) & NoAxes()

ggsave(
  filename = file.path(out_dir, "S9A.tiff"),
  plot     = p,
  device   = "tiff",
  width    = 12,
  height   = 6,
  units    = "in",
  dpi      = 300
)
