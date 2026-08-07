# 05_tage_spread_figure.R
#
# Figure: per-sample tAge spread by timepoint, faceted by cell type, with a
# global Kruskal-Wallis test per facet and pairwise Wilcoxon brackets
# (each timepoint vs 'none'). Uses the tage_temporal_by_celltype.csv already
# computed by 04_tage_by_celltype.R (scale_counts=FALSE, corrected version).

source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})

tage_temporal <- read.csv(file.path(RERUN_DIR, "tage_temporal_by_celltype.csv"))
tage_temporal$time_after_treatment <- factor(tage_temporal$time_after_treatment,
                                              levels = c("none", "4_days", "10_days", "20_days"))
tage_temporal$cell_type <- factor(tage_temporal$cell_type,
                                   levels = c("Fibroblast", "Keratinocyte", "Melanocyte"))

make_figure <- function(model_col, model_label, file_name) {
  df <- tage_temporal
  df$tAge <- df[[model_col]]

  # vs-'none' pairwise Wilcoxon brackets, BH-adjusted (geom_pwc, not
  # stat_compare_means(comparisons=...) -- the latter only ever shows raw
  # p-values in this ggpubr version despite accepting p.adjust.method)
  p <- ggplot(df, aes(x = time_after_treatment, y = tAge)) +
    geom_violin(aes(fill = time_after_treatment), alpha = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    geom_jitter(width = 0.08, size = 1.3, alpha = 0.6) +
    facet_wrap(~cell_type, nrow = 1) +
    geom_pwc(
      method = "wilcox.test", ref.group = "none", p.adjust.method = "BH",
      label = "p.signif", tip.length = 0.01, step.increase = 0.08
    ) +
    stat_compare_means(
      method = "kruskal.test", label.y.npc = "top", label.x.npc = "left",
      aes(label = paste0("KW ", after_stat(p.format)))
    ) +
    scale_fill_brewer(palette = "Blues") +
    theme_bw() +
    theme(legend.position = "none", strip.text = element_text(face = "bold")) +
    labs(
      x = NULL, y = paste0("tAge (", model_label, ")"),
      title = paste0("Per-sample transcriptomic age by timepoint post-irradiation (", model_label, " model)"),
      subtitle = paste0(
        "Points = individual samples (n=6/timepoint/cell type). Brackets: pairwise Wilcoxon vs 'none', BH-adjusted (ns/*/**/*** = p≥.05/<.05/<.01/<.001).\n",
        "KW = Kruskal-Wallis across all 4 timepoints, per cell type."
      )
    )

  ggsave(file.path(RERUN_DIR, file_name), p, width = 12, height = 5.5, dpi = 300)
  cat(sprintf("Saved %s\n", file.path(RERUN_DIR, file_name)))
  p
}

p1 <- make_figure("scaled_diff_EN_tAge", "scaled_diff EN model", "tage_temporal_spread_figure_scaleddiff.png")
p2 <- make_figure("yugene_diff_EN_tAge", "yugene_diff EN model", "tage_temporal_spread_figure_yugenediff.png")

cat("\nDone.\n")
