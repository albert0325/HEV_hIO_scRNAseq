# =============================================================================
# 04_ISG_response.R
# =============================================================================

source("R/utils.R")
library(Seurat)
library(tidyverse)

out_dir <- "outputs/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Step 1: Load the mock + HEV-norm annotated object
# -----------------------------------------------------------------------------

seurat_obj <- readRDS("data/hev_mock_norm_annotated.rds")

# Assign treatment labels from library identity
seurat_obj$treatments <- case_when(
  seurat_obj$orig.ident == "hev"     ~ "hev",
  seurat_obj$orig.ident == "IO_mock" ~ "mock",
  TRUE                               ~ NA_character_
)

# -----------------------------------------------------------------------------
# Step 2: Compute ISG module score (Schoggins et al.)
# -----------------------------------------------------------------------------

sig_schoggins <- read.csv("data/isg_Schoggins.csv", stringsAsFactors = FALSE) %>%
  pull(gene)

seurat_obj <- AddModuleScore(
  object   = seurat_obj,
  features = list(sig_schoggins),
  ctrl     = 5,
  assay    = "SCT",
  slot     = "data",
  name     = "sig_schoggins"
)

# AddModuleScore appends a numeric suffix; strip it for cleaner column naming
score_cols <- grep("^sig_", colnames(seurat_obj@meta.data), value = TRUE)
colnames(seurat_obj@meta.data)[match(score_cols, colnames(seurat_obj@meta.data))] <-
  sub("[0-9]+$", "", score_cols)

# -----------------------------------------------------------------------------
# Step 3: Visualise and test ISG scores per cell type
# -----------------------------------------------------------------------------

plot_signature_score(
  seurat_obj      = seurat_obj,
  score_colname   = "sig_schoggins",
  score_label     = "Schoggins ISG signature score",
  celltype_colname= "SingleR_label",
  out_dir         = out_dir
)
