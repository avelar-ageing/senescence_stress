# 06_temporal_universal_tage_figure.R
#
# Figure for the temporal-tAge results section: tAge trajectory over time
# post-irradiation, both EN models, faceted by cell type, styled to match
# meta_analysis/08_universal_tage_figure.R (violin+jitter only, brackets vs
# 'none' with non-significant ones hidden, no title/subtitle, large text).
# Reads only the already-computed outputs of 04_tage_by_celltype.R -- does
# not re-run tAge prediction.

source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
})

TIME_LEVELS <- c("none", "4_days", "10_days", "20_days")
CT_LEVELS <- c("Fibroblast", "Keratinocyte", "Melanocyte")

tage_temporal <- read.csv(file.path(RERUN_DIR, "tage_temporal_by_celltype.csv"))
tage_temporal$time_after_treatment <- factor(tage_temporal$time_after_treatment, levels = TIME_LEVELS)
tage_temporal$cell_type <- factor(tage_temporal$cell_type, levels = CT_LEVELS)
# All-pairs family (07_tage_temporal_pairwise.R): every timepoint vs every
# other, BH-adjusted across 3 cell types x 6 pairs x 2 models = 36 tests.
# This is the canonical family -- it subsumes the vs-baseline comparisons
# (none vs 4/10/20 days are 3 of the 6 pairs) and additionally covers the
# between-timepoint contrasts the trajectory claims depend on.
pairwise <- read.csv(file.path(RERUN_DIR, "tage_temporal_pairwise_all_timepoints.csv"))

plot_df <- tage_temporal %>%
  dplyr::select(cell_type, time_after_treatment, scaled_diff_EN_tAge, yugene_diff_EN_tAge) %>%
  tidyr::pivot_longer(cols = c(scaled_diff_EN_tAge, yugene_diff_EN_tAge),
                       names_to = "model", values_to = "tAge") %>%
  mutate(model = ifelse(model == "scaled_diff_EN_tAge", "scaled_diff", "yugene_diff"),
         model_label = ifelse(model == "scaled_diff", "Scaled Difference", "YuGene"))

sig_symbol <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

# geom_violin(trim=FALSE) draws kernel-density tails that extend past the real
# data range, so max(tAge) understates how tall each violin actually renders
# (this is what left brackets overlapping the Keratinocyte violins). Build the
# violin layer once and read its true rendered y extent per model row.
violin_extent <- function(model_name) {
  d <- plot_df[plot_df$model == model_name, ]
  gb <- ggplot_build(
    ggplot(d, aes(x = time_after_treatment, y = tAge)) +
      geom_violin(trim = FALSE) + facet_wrap(~cell_type)
  )
  max(gb$data[[1]]$ymax, na.rm = TRUE)
}

# Brackets: every significant timepoint pair, per cell type, per model, from
# the 36-test all-pairs BH family. Non-significant comparisons dropped
# entirely. Brackets are stacked shortest-span-first so short contrasts sit
# low and wide ones sit above them, minimising visual crossing.
make_stat_df <- function(model_name) {
  sub <- pairwise[pairwise$model == model_name, ]
  sub$label <- sig_symbol(sub$p.adj)
  sub <- sub[sub$label != "ns", ]
  sub$cell_type <- factor(sub$cell_type, levels = CT_LEVELS)
  sub$span <- abs(match(sub$timepoint_2, TIME_LEVELS) - match(sub$timepoint_1, TIME_LEVELS))
  sub <- sub[order(sub$cell_type, sub$span,
                    match(sub$timepoint_1, TIME_LEVELS)), ]

  # facet_grid(scales="free_y") gives each ROW (model) one shared y scale, so
  # bracket heights must come from that shared row-wide extent -- not from each
  # cell type's own range, which put panels on different scales (inconsistent
  # spacing) and placed brackets below neighbouring violins that extend higher
  # on the shared axis (overlap). Uses the rendered violin extent, not the data
  # max, so the trim=FALSE density tails are cleared too.
  row_vals <- plot_df$tAge[plot_df$model == model_name]
  y_max <- violin_extent(model_name); y_range <- diff(range(row_vals))
  sub <- sub %>% group_by(cell_type) %>%
    mutate(y.position = y_max + y_range * (0.06 + 0.14 * row_number())) %>% ungroup()

  data.frame(model = model_name,
             model_label = ifelse(model_name == "scaled_diff", "Scaled Difference", "YuGene"),
             cell_type = sub$cell_type, group1 = sub$timepoint_1, group2 = sub$timepoint_2,
             label = sub$label, y.position = sub$y.position)
}
stat_df <- rbind(make_stat_df("scaled_diff"), make_stat_df("yugene_diff"))

p <- ggplot(plot_df, aes(x = time_after_treatment, y = tAge, fill = time_after_treatment)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_jitter(width = 0.08, size = 1.2, alpha = 0.5, colour = "black") +
  # Median trajectory across timepoints. Median (not mean) to match the
  # Wilcoxon tests and the median differences reported in
  # tage_temporal_pairwise_all_timepoints.csv. group=1 connects across the
  # discrete x within each panel.
  stat_summary(aes(group = 1), fun = median, geom = "line",
               colour = "grey20", linewidth = 0.9) +
  stat_summary(aes(group = 1), fun = median, geom = "point",
               colour = "grey20", size = 2.4) +
  stat_pvalue_manual(stat_df, label = "label", xmin = "group1", xmax = "group2",
                      y.position = "y.position", tip.length = 0, bracket.size = 0.5, size = 5) +
  # Relabel for display only -- the underlying factor levels stay as
  # none/4_days/... because the bracket table matches against those exact
  # strings from tage_temporal_pairwise_all_timepoints.csv.
  scale_x_discrete(labels = c("none" = "None", "4_days" = "4 Days",
                              "10_days" = "10 Days", "20_days" = "20 Days")) +
  facet_grid(model_label ~ cell_type, scales = "free_y") +
  theme_bw(base_size = 18) +
  theme(legend.position = "none",
        strip.text = element_text(face = "plain", size = 18),
        strip.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_text(size = 18)) +
  labs(x = "Days Post-Irradiation", y = "tAge")

ggsave(file.path(RERUN_DIR, "figure_temporal_universal_tage.png"), p, width = 13, height = 10, dpi = 300)
cat(sprintf("Saved -> %s\n", file.path(RERUN_DIR, "figure_temporal_universal_tage.png")))
