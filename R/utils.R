# =============================================================================
# utils.R
# Shared colour palettes and helper functions
#
# Source this file at the top of each analysis script:
#   source("R/utils.R")
#
# =============================================================================

library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)

# -----------------------------------------------------------------------------
# Colour palettes
# Cell type colours are consistent across all figures.
# -----------------------------------------------------------------------------

custom_colors <- c(
  "TA"                  = "#7C1D6F",
  "Enterocyte"          = "#FCE498",
  "Stem cell"           = "#D78287",
  "Enteroendocrine cell"= "#A8D8EA",
  "Goblet cell"         = "#B7C370",
  "Paneth cell"         = "#008042"
)

# -----------------------------------------------------------------------------
# plot_signature_score()
#
# Generates a boxplot comparing a gene signature score across cell types
# and treatment conditions, with Wilcoxon rank-sum statistics relative to mock.
#
# Arguments:
#   seurat_obj      : Seurat object containing metadata columns below.
#   score_colname   : Column name in seurat_obj@meta.data holding the score.
#   score_label     : Human-readable label for the y-axis and plot title.
#   celltype_colname: Column name for cell type labels (default: "SingleR_label").
#   out_dir         : Output directory path for saving the PDF.
#
# Returns:
#   A named list: list(overlay_plot = p2), where p2 is the ggplot object.
#   Returns NULL for overlay_plot when only one treatment is present.
# -----------------------------------------------------------------------------

plot_signature_score <- function(seurat_obj,
                                 score_colname,
                                 score_label,
                                 celltype_colname = "SingleR_label",
                                 out_dir) {

  plot_data <- seurat_obj@meta.data %>%
    dplyr::select(all_of(celltype_colname), treatments, all_of(score_colname)) %>%
    dplyr::rename(
      signature_score = all_of(score_colname),
      celltype        = all_of(celltype_colname)
    )

  # Place "mock" first; all other conditions follow in their natural order
  all_treatments     <- unique(plot_data$treatments)
  ordered_treatments <- c("mock", setdiff(all_treatments, "mock"))
  plot_data$treatments <- factor(plot_data$treatments, levels = ordered_treatments)

  # Order cell types by ascending mean score in the mock condition
  mock_order <- plot_data %>%
    dplyr::filter(treatments == "mock") %>%
    group_by(celltype) %>%
    dplyr::summarise(mean_score = mean(signature_score, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_score) %>%
    pull(celltype)

  plot_data$celltype <- factor(plot_data$celltype, levels = mock_order)

  # Only produce a multi-treatment comparison plot when >1 condition is present
  p2 <- NULL
  n_treatments <- length(ordered_treatments)

  if (n_treatments > 1) {

    # Assign colours: grey for mock, then colorblind-friendly palette
    if (n_treatments == 2) {
      treatment_colors <- c("grey", "#E69F00")
    } else if (n_treatments == 3) {
      treatment_colors <- c("grey", "#E69F00", "#56B4E9")
    } else {
      extra_colors     <- rainbow(n_treatments - 3, start = 0.3, end = 0.9)
      treatment_colors <- c("grey", "#E69F00", "#56B4E9", extra_colors)
    }
    names(treatment_colors) <- ordered_treatments

    # Wilcoxon rank-sum test vs. mock, Bonferroni-adjusted
    stat_test <- plot_data %>%
      group_by(celltype) %>%
      wilcox_test(signature_score ~ treatments, ref.group = "mock", alternative = "less") %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance() %>%
      add_xy_position(x = "celltype", dodge = 0.8)

    p2 <- ggplot(plot_data, aes(x = celltype, y = signature_score, fill = treatments)) +
      geom_boxplot(
        alpha        = 1,
        outlier.size = 0.5,
        color        = "black",
        linewidth    = 0.5,
        width        = 0.5
      ) +
      stat_pvalue_manual(
        stat_test,
        label         = "p.adj.signif",
        tip.length    = 0.01,
        step.increase = 0.1,
        size          = 6
      ) +
      scale_fill_manual(values = treatment_colors) +
      labs(
        title = paste(score_label, "— Treatment Comparison"),
        x     = "Cell type",
        y     = score_label,
        fill  = "Treatment"
      ) +
      theme_bw() +
      theme(
        strip.placement   = "outside",
        strip.background  = element_blank(),
        axis.text.x       = element_text(angle = 45, hjust = 1, size = 18),
        axis.title.x      = element_text(size = 15),
        axis.text.y       = element_text(size = 15),
        axis.title.y      = element_text(size = 15),
        legend.position   = "top",
        legend.text       = element_text(size = 16),
        legend.title      = element_text(size = 16),
        panel.spacing     = unit(0.6, "lines"),
        aspect.ratio      = 1,
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank()
      )

    ggsave(
      filename = file.path(out_dir, paste0("boxplot_", score_colname, ".pdf")),
      plot     = p2,
      device   = "pdf",
      width    = 11,
      height   = 7,
      units    = "in",
      dpi      = 300
    )
  }

  return(list(overlay_plot = p2))
}
