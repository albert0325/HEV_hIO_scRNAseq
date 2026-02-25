# =============================================================================
# 02_annotation.R
# Cell type annotation using a custom SingleR reference built from
# curated intestinal marker genes, followed by UMAP visualisation
# and HEV infection rate quantification per cell type.
#
# The marker gene table (Suppl. Table S1) is provided as an Excel file
# in the data/ directory. SingleR scores are visualised as a heatmap
# to allow inspection of annotation confidence.
#
# Input:  hev              — Seurat object from 01_quality_control.R
#         data/hIOs_marker_v3_singleR.xlsx
# Output: outputs/fig3A.tiff  — Dot plot of marker gene expression
#         outputs/fig3B.tiff  — UMAP coloured by cell type
#         outputs/fig3C.tiff  — UMAP coloured by infection status
#         outputs/fig3D.tiff  — Infection rate per cell type (bar chart)
#         outputs/S9B.tiff    — Violin plots of QC metrics by cell type
#         outputs/S9C.tiff    — SingleR score heatmap
#
# =============================================================================

source("R/utils.R")
library(Seurat)
library(SingleR)
library(SummarizedExperiment)
library(readxl)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(tidyverse)

out_dir <- "outputs/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Step 1: Build a custom pseudo-bulk reference from curated marker genes
# Each cell type is represented as a binary vector (1 = marker, 0 = absent).
# This encodes the domain knowledge in Suppl. Table S1 into a format
# compatible with SingleR's reference-based classification.
# -----------------------------------------------------------------------------

marker_df   <- read_excel("data/hIOs_marker_v3_singleR.xlsx")
marker_list <- split(marker_df$Marker, marker_df$CellType)
all_genes   <- unique(marker_df$Marker)

ref_mat <- sapply(marker_list, function(genes) as.integer(all_genes %in% genes))
rownames(ref_mat) <- all_genes

ref_se <- SummarizedExperiment(
  assays  = list(logcounts = ref_mat),
  colData = data.frame(label = colnames(ref_mat))
)

# -----------------------------------------------------------------------------
# Step 2: Run SingleR annotation
# SCT-normalised expression is used as the query; classification is
# performed against the custom intestinal reference built above.
# -----------------------------------------------------------------------------

query_expr      <- GetAssayData(hev, layer = "data", assay = "SCT")
singleR_results <- SingleR(
  test   = as.matrix(query_expr),
  ref    = ref_se,
  labels = ref_se$label
)

hev$SingleR_label <- singleR_results$labels
all_celltypes     <- unique(hev$SingleR_label)
colors            <- custom_colors[intersect(names(custom_colors), all_celltypes)]

# -----------------------------------------------------------------------------
# Step 3: SingleR score heatmap (Suppl. Fig. S9C)
# Scores are scaled across cell types to highlight relative confidence.
# Cells are split by their pruned label to reveal classification clarity.
# -----------------------------------------------------------------------------

scores_mat_t <- t(scale(singleR_results$scores))
score_range  <- quantile(scores_mat_t, probs = c(0.15, 0.98), na.rm = TRUE)
col_fun      <- colorRamp2(
  seq(score_range[1], score_range[2], length.out = 100),
  viridis(100)
)

cell_labels <- as.character(singleR_results$pruned.labels)
top_anno    <- HeatmapAnnotation(
  Labels                = cell_labels,
  col                   = list(Labels = custom_colors),
  show_annotation_name  = FALSE
)

ht <- Heatmap(
  matrix           = scores_mat_t,
  name             = "score",
  col              = col_fun,
  top_annotation   = top_anno,
  show_column_names= FALSE,
  show_row_names   = TRUE,
  cluster_columns  = FALSE,
  cluster_rows     = TRUE,
  column_split     = cell_labels,
  column_title     = " ",
  row_title        = "",
  heatmap_legend_param = list(title = "Score")
)

tiff(file.path(out_dir, "S9C.tiff"), width = 4000, height = 2500, res = 300)
draw(ht)
dev.off()

# -----------------------------------------------------------------------------
# Step 4: Marker gene dot plot (Fig. 3A)
# Genes are ordered to mirror the cell type groupings in Suppl. Table S1.
# scale = FALSE preserves raw average expression values on the colour axis,
# making cross-cell-type comparisons interpretable.
# -----------------------------------------------------------------------------

genes_of_interest <- unique(c(
  # Enterocyte (EC)
  "CDX2", "FABP1", "CA2", "KRT20", "MTTP", "ACE2", "HNF4A",
  # TA (transit-amplifying)
  "MKI67", "PCNA", "TOP2A", "CDK1", "CCNB1", "CCNB2", "CCNA2", "AURKA",
  # Goblet
  "MUC2", "ATOH1", "SPINK4", "FCGBP", "REG4", "CLCA1",
  # Paneth
  "LYZ", "MMP7", "PLA2G2A", "SGK1", "KLF4", "DUOX2",
  # Enteroendocrine (EE)
  "CHGA", "CHGB", "NEUROD1", "NEUROG3", "TPH1", "SCG2",
  # Intestinal Stem cell (ISC)
  "LGR5", "OLFM4", "ASCL2", "SOX9", "HES1", "CD44"
))
genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(hev)]

Idents(hev)       <- "SingleR_label"
DefaultAssay(hev) <- "integrated"
desired_order     <- unique(marker_df$CellType)

p <- DotPlot(
  hev,
  features = genes_of_interest,
  group.by = "SingleR_label",
  dot.scale = 12,
  scale     = FALSE
) +
  scale_color_viridis_c(option = "viridis", limits = c(0, 12), oob = scales::squish) +
  scale_y_discrete(limits = desired_order) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 17),
    axis.text.y  = element_text(size = 15, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

ggsave(file.path(out_dir, "fig3A.tiff"), plot = p, device = "tiff",
       width = 20, height = 6, units = "in", dpi = 300)

# -----------------------------------------------------------------------------
# Step 5: QC violin plots by cell type (Suppl. Fig. S9B)
# Plotting post-annotation confirms that QC metric distributions are
# comparable across annotated cell types, ruling out annotation bias
# driven by residual low-quality cells.
# -----------------------------------------------------------------------------

p <- VlnPlot(hev, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             ncol = 3, pt.size = 0)

ggsave(file.path(out_dir, "S9B.tiff"), plot = p, device = "tiff",
       width = 18, height = 6, units = "in", dpi = 300)

# -----------------------------------------------------------------------------
# Step 6: UMAP coloured by cell type (Fig. 3B)
# -----------------------------------------------------------------------------

p <- DimPlot(hev, group.by = "SingleR_label", reduction = "umap", cols = colors) +
  NoAxes()

ggsave(file.path(out_dir, "fig3B.tiff"), plot = p, device = "tiff",
       width = 8, height = 6, units = "in", dpi = 300)

# -----------------------------------------------------------------------------
# Step 7: UMAP coloured by HEV infection status (Fig. 3C)
# -----------------------------------------------------------------------------

p <- DimPlot(
  hev,
  group.by  = "infected_status",
  reduction = "umap",
  cols      = c("Positive" = "red", "Negative" = "grey")
) + NoAxes()

ggsave(file.path(out_dir, "fig3C.tiff"), plot = p, device = "tiff",
       width = 6, height = 6, units = "in", dpi = 300)

# -----------------------------------------------------------------------------
# Step 8: HEV infection rate per cell type (Fig. 3D)
# Proportions are computed from all cells (infected + uninfected) to avoid
# ascertainment bias. Labels show both percentage and absolute counts.
# -----------------------------------------------------------------------------

infection_rate <- hev@meta.data %>%
  group_by(SingleR_label) %>%
  summarise(
    total_cells       = n(),
    infected_cells    = sum(infected),
    proportion_infected = infected_cells / total_cells,
    .groups = "drop"
  ) %>%
  arrange(desc(proportion_infected))

p <- ggplot(infection_rate,
            aes(x = reorder(SingleR_label, proportion_infected),
                y = proportion_infected,
                fill = SingleR_label)) +
  geom_col() +
  geom_text(
    aes(label = sprintf("%.1f%%\n(%d/%d)",
                        proportion_infected * 100,
                        infected_cells,
                        total_cells)),
    hjust = -0.1, size = 5
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, max(infection_rate$proportion_infected) * 1.2)
  ) +
  scale_fill_manual(values = custom_colors) +
  labs(
    x     = "Cell type",
    y     = "Proportion infected",
    title = "HEV infection rate by cell type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    axis.text.y  = element_text(size = 20),
    axis.title.y = element_text(size = 25),
    legend.position = "none"
  )

ggsave(file.path(out_dir, "fig3D.tiff"), plot = p, device = "tiff",
       width = 10, height = 6, units = "in", dpi = 300)
