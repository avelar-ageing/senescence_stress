# 05_tage_spread_figure.R
#
# Figure: per-sample tAge spread by timepoint, faceted by cell type, with a
# global Kruskal-Wallis test per facet and pairwise Wilcoxon brackets
# (each timepoint vs 'none'). Uses the tage_temporal_by_celltype.csv already
# computed by 04_tage_by_celltype.R (scale_counts=FALSE, corrected version).
#
# P-VALUE CORRECTION: brackets use the p.adj column from
# tage_temporal_wilcoxon_vs_none.csv -- i.e. BH-adjusted across the FULL
# 18-test family (3 cell types x 3 timepoints x 2 models), the same
# correction reported in that CSV -- NOT geom_pwc()'s/stat_compare_means()'s
# own per-facet (3-test) adjustment, which gives different, less conservative
# numbers for the same comparisons. Uses stat_pvalue_manual() with
# precomputed p.adj so the figure and the CSV table always agree.

source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
})

tage_temporal <- read.csv(file.path(RERUN_DIR, "tage_temporal_by_celltype.csv"))
tage_temporal$time_after_treatment <- factor(tage_temporal$time_after_treatment,
                                              levels = c("none", "4_days", "10_days", "20_days"))
tage_temporal$cell_type <- factor(tage_temporal$cell_type,
                                   levels = c("Fibroblast", "Keratinocyte", "Melanocyte"))

# Precomputed, family-wise (18-test) BH-corrected p-values -- same numbers as
# reported in tage_temporal_wilcoxon_vs_none.csv, not recomputed here.
wilcox_adj <- read.csv(file.path(RERUN_DIR, "tage_temporal_wilcoxon_vs_none.csv"))

kw_pvals <- tage_temporal %>%
  group_by(cell_type) %>%
  summarise(scaled_diff = kruskal.test(scaled_diff_EN_tAge ~ time_after_treatment)$p.value,
            yugene_diff = kruskal.test(yugene_diff_EN_tAge ~ time_after_treatment)$p.value,
            .groups = "drop")

make_figure <- function(model_col, model_name, model_label, file_name) {
  df <- tage_temporal
  df$tAge <- df[[model_col]]

  # Bracket data: one row per (cell_type, timepoint), 'none' vs that timepoint,
  # using the family-wise-adjusted p.adj already computed for this exact model.
  bracket_df <- wilcox_adj %>%
    filter(model == model_name) %>%
    transmute(cell_type = factor(cell_type, levels = levels(tage_temporal$cell_type)),
              group1 = "none", group2 = timepoint, p.adj = p.adj) %>%
    add_significance("p.adj") %>%
    arrange(cell_type, group2)

  # Stack bracket heights per facet, above that facet's max point.
  y_max <- df %>% group_by(cell_type) %>% summarise(ymax = max(tAge), .groups = "drop")
  bracket_df <- bracket_df %>%
    left_join(y_max, by = "cell_type") %>%
    group_by(cell_type) %>%
    mutate(step = row_number(), y.position = ymax + step * (0.12 * diff(range(df$tAge)))) %>%
    ungroup()

  kw_labels <- kw_pvals %>%
    transmute(cell_type, label = paste0("KW ", formatC(.data[[model_name]], format = "g", digits = 2)))
  kw_labels <- kw_labels %>% left_join(y_max, by = "cell_type") %>%
    mutate(y.position = ymax + 0.12 * diff(range(df$tAge)) * 0.4)

  p <- ggplot(df, aes(x = time_after_treatment, y = tAge)) +
    geom_violin(aes(fill = time_after_treatment), alpha = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    geom_jitter(width = 0.08, size = 1.3, alpha = 0.6) +
    facet_wrap(~cell_type, nrow = 1) +
    stat_pvalue_manual(bracket_df, label = "p.adj.signif", tip.length = 0.01) +
    geom_text(data = kw_labels, aes(x = 1.5, y = y.position, label = label),
              inherit.aes = FALSE, hjust = 0, size = 3.5) +
    scale_fill_brewer(palette = "Blues") +
    theme_bw() +
    theme(legend.position = "none", strip.text = element_text(face = "bold")) +
    labs(
      x = NULL, y = paste0("tAge (", model_label, ")"),
      title = paste0("Per-sample transcriptomic age by timepoint post-irradiation (", model_label, " model)"),
      subtitle = paste0(
        "Points = individual samples (n=6/timepoint/cell type). Brackets: pairwise Wilcoxon vs 'none', BH-adjusted across all\n",
        "18 tests (3 cell types x 3 timepoints x 2 models) -- same p.adj as tage_temporal_wilcoxon_vs_none.csv (ns/*/**/*** = p≥.05/<.05/<.01/<.001).\n",
        "KW = Kruskal-Wallis across all 4 timepoints, per cell type (uncorrected, one test per facet)."
      )
    )

  ggsave(file.path(RERUN_DIR, file_name), p, width = 12, height = 6.2, dpi = 300)
  cat(sprintf("Saved %s\n", file.path(RERUN_DIR, file_name)))
  print(bracket_df[, c("cell_type", "group1", "group2", "p.adj", "p.adj.signif")])
  p
}

p1 <- make_figure("scaled_diff_EN_tAge", "scaled_diff", "scaled_diff EN model", "tage_temporal_spread_figure_scaleddiff.png")
p2 <- make_figure("yugene_diff_EN_tAge", "yugene_diff", "yugene_diff EN model", "tage_temporal_spread_figure_yugenediff.png")

cat("\nDone.\n")
