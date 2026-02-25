# =============================================================================
# 03_stemness_pseudotime.R
# Developmental potency scoring with CytoTRACE2 and pseudotime trajectory
# inference with Monocle3.
#
# CytoTRACE2 estimates each cell's differentiation status from transcriptional
# diversity; higher scores indicate greater stemness / developmental potency.
# Pseudotime is computed independently per library (hev, hevpp) using their
# respective UMAP embeddings, then merged back into the joint Seurat object.
# Root cells are defined as the top-10 stem cells ranked by CytoTRACE2 score.
#
# Input:  hev              — Seurat object from 02_annotation.R
# Output: outputs/S11A.tiff — UMAP coloured by CytoTRACE2 score
#         outputs/S11B.tiff — UMAP coloured by CytoTRACE2 potency category
#         outputs/S11C.tiff — UMAP coloured by Monocle3 pseudotime
# =============================================================================

source("R/utils.R")
library(Seurat)
library(CytoTRACE2)
library(monocle3)
library(tidyverse)
library(viridis)

out_dir <- "outputs/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Step 1: CytoTRACE2 stemness scoring
# SCT counts are used as input; CytoTRACE2 is run in gene-expression matrix
# mode (is_seurat = FALSE) to ensure compatibility with the SCT assay layer.
# -----------------------------------------------------------------------------

Idents(hev) <- "SingleR_label"

sct_data <- LayerData(hev, assay = "SCT", layer = "counts")
cytotrace2_result <- cytotrace2(sct_data, species = "human", is_seurat = FALSE)

hev$CytoTRACE2_Score    <- cytotrace2_result$CytoTRACE2_Score
hev$CytoTRACE2_Potency  <- cytotrace2_result$CytoTRACE2_Potency
hev$CytoTRACE2_Relative <- cytotrace2_result$CytoTRACE2_Relative

# Suppl. Fig. S11A — Continuous stemness score
p1 <- FeaturePlot(hev, features = "CytoTRACE2_Score", reduction = "umap") +
  ggtitle("CytoTRACE2 stemness score") &
  NoAxes()

ggsave(file.path(out_dir, "S11A.tiff"), plot = p1, device = "tiff",
       width = 8, height = 6, units = "in", dpi = 300)

# Suppl. Fig. S11B — Discrete potency categories
p2 <- DimPlot(hev, group.by = "CytoTRACE2_Potency", reduction = "umap") &
  NoAxes()

ggsave(file.path(out_dir, "S11B.tiff"), plot = p2, device = "tiff",
       width = 8, height = 6, units = "in", dpi = 300)

# -----------------------------------------------------------------------------
# Step 2: Monocle3 pseudotime trajectory
#
# Pseudotime is computed separately for each library to avoid cross-condition
# trajectory distortions. UMAP coordinates are transferred from the integrated
# Seurat object so that trajectories are anchored to the same embedding.
# Root cells: top 10 stem cells by CytoTRACE2 score, providing a biologically
# grounded and data-driven starting point for the trajectory.
# -----------------------------------------------------------------------------

run_pseudotime <- function(seurat_obj, counts_layer) {

  Idents(seurat_obj) <- "SingleR_label"

  counts_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = counts_layer)
  cell_metadata <- seurat_obj@meta.data[colnames(counts_matrix), , drop = FALSE]
  gene_metadata <- data.frame(gene_short_name = rownames(counts_matrix),
                              row.names       = rownames(counts_matrix))

  cds <- new_cell_data_set(
    counts_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )

  # Transfer the integrated UMAP embedding to ensure trajectory consistency
  umap_coords           <- Embeddings(seurat_obj, reduction = "umap")[colnames(cds), ]
  reducedDims(cds)$UMAP <- umap_coords
  colData(cds)$cell_type <- Idents(seurat_obj)

  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds)

  # Select root cells: top-10 stem cells ranked by CytoTRACE2 score
  stem_idx <- which(colData(cds)$cell_type == "Stem cell")
  if (length(stem_idx) == 0) {
    stop("No 'Stem cell' cells found in this subset — cannot define trajectory root.")
  }
  stem_scores    <- colData(cds)$CytoTRACE2_Score[stem_idx]
  top_root_cells <- colnames(cds)[stem_idx[order(stem_scores, decreasing = TRUE)][seq_len(min(10, length(stem_idx)))]]

  cds <- order_cells(cds, root_cells = top_root_cells)

  pt <- pseudotime(cds)
  pt[is.infinite(pt)] <- NA
  return(pt)
}

# Run per library and merge into a single vector aligned to the full object
DefaultAssay(hev) <- "RNA"
hev.list          <- SplitObject(hev, split.by = "orig.ident")

pseudotime_all           <- setNames(rep(NA_real_, ncol(hev)), colnames(hev))
pt_hev                   <- run_pseudotime(hev.list$hev,   counts_layer = "counts.hev.1")
pt_hevpp                 <- run_pseudotime(hev.list$hevpp, counts_layer = "counts.hevpp.2")
pseudotime_all[names(pt_hev)]   <- pt_hev
pseudotime_all[names(pt_hevpp)] <- pt_hevpp

hev$pseudotime <- pseudotime_all

# Suppl. Fig. S11C — Pseudotime on UMAP
p3 <- FeaturePlot(hev, features = "pseudotime", reduction = "umap") +
  scale_color_viridis_c(option = "magma", na.value = "lightgrey") +
  ggtitle("Monocle3 pseudotime") &
  NoAxes()

ggsave(file.path(out_dir, "S11C.tiff"), plot = p3, device = "tiff",
       width = 8, height = 6, units = "in", dpi = 300)
