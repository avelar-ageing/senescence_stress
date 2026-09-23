# 06_temporal_universal_tage_figure.R
#
# Figure for the temporal-tAge results section: tAge trajectory over time
# post-irradiation, faceted by cell type, styled to match
# meta_analysis/08_universal_tage_figure.R (violin+jitter only, brackets vs
# 'none' with non-significant ones hidden, no title/subtitle, large text).
# Reads only already-computed per-sample tAge -- does not re-run prediction.
#
# UPDATED 2026-08-25: the MORTALITY clock is now a third row. Section 2.2.3 is a
# whole-transcriptome section and so reports both instruments; the figure was
# built before the mortality clock existed and showed only the two chronological
# models, which no longer matches the text.
#
# TWO BH FAMILIES, DELIBERATELY SEPARATE. The chronological brackets come from the
# 36-test all-pairs family of 07_tage_temporal_pairwise.R (3 cell types x 6 pairs x
# 2 models); the mortality brackets come from the 18-test family of
# 08_mortality_temporal.py (3 x 6 x 1). Each clock's stars are adjusted within its
# own analysis, which is how they are reported in the text; they are NOT pooled
# into one 54-test family here, because the two analyses were run and reported
# separately and re-pooling them would change published p-values.
#
# P-VALUE FLOOR. Every comparison is 6 against 6, so the smallest attainable
# two-sided Wilcoxon p is 2/choose(12,6) = 0.00216 and no bracket can ever show
# "***" after adjustment. Absent stars mean "not resolvable at n = 6", not "no
# effect".

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
# Plot what the tests describe: the keratinocyte batch is removed from the
# chronological clocks here exactly as in 07_tage_temporal_pairwise.R, so the
# points and the significance brackets refer to the same values. The mortality
# panel needs nothing - 08_mortality_temporal.py already writes the centred
# value into mortality_tAge (mortality_tAge_raw keeps the uncorrected one).
source(file.path("temporal_analysis", "R_keratinocyte_batch.R"))
tage_temporal <- batch_centre_keratinocytes(
  tage_temporal, c("scaled_diff_EN_tAge", "yugene_diff_EN_tAge"))
mort_temporal <- read.csv(file.path(RERUN_DIR, "mortality_temporal_by_celltype.csv"))
# All-pairs family (07_tage_temporal_pairwise.R): every timepoint vs every
# other, BH-adjusted across 3 cell types x 6 pairs x 2 models = 36 tests.
# This is the canonical family -- it subsumes the vs-baseline comparisons
# (none vs 4/10/20 days are 3 of the 6 pairs) and additionally covers the
# between-timepoint contrasts the trajectory claims depend on.
# One consolidated test table for all three clocks, pooled into a single BH
# family of 54 by 07_tage_temporal_tests.R (previously two files with two
# different adjusted floors for the same 6v6 design).
pairwise <- read.csv(file.path(RERUN_DIR, "tage_temporal_tests.csv"))
pairwise <- pairwise[pairwise$test == "pairwise_timepoints", ]
pairwise <- data.frame(cell_type = pairwise$cell_type,
                       timepoint_1 = pairwise$group_1,
                       timepoint_2 = pairwise$group_2,
                       model = pairwise$model, p.adj = pairwise$p_adj,
                       diff = pairwise$diff)

MODEL_LEVELS <- c("scaled_diff", "yugene_diff", "mortality")
# CLOCK ORDER (2026-09-01): YuGene, Scaled Difference, Mortality everywhere -
# figures, printed tables and the manuscript text all follow this order.
MODEL_LABELS <- c(yugene_diff = "Chronological, YuGene",
                  scaled_diff = "Chronological, Scaled Difference",
                  mortality = "Mortality")

plot_df <- bind_rows(
  tage_temporal %>%
    dplyr::select(cell_type, time_after_treatment,
                  scaled_diff_EN_tAge, yugene_diff_EN_tAge) %>%
    tidyr::pivot_longer(cols = c(scaled_diff_EN_tAge, yugene_diff_EN_tAge),
                        names_to = "model", values_to = "tAge") %>%
    mutate(model = ifelse(model == "scaled_diff_EN_tAge", "scaled_diff", "yugene_diff")),
  mort_temporal %>%
    transmute(cell_type, time_after_treatment = timepoint,
              model = "mortality", tAge = mortality_tAge)) %>%
  mutate(time_after_treatment = factor(time_after_treatment, levels = TIME_LEVELS),
         cell_type = factor(cell_type, levels = CT_LEVELS),
         model = factor(model, levels = MODEL_LEVELS),
         model_label = factor(MODEL_LABELS[as.character(model)],
                              levels = MODEL_LABELS))
stopifnot(!any(is.na(plot_df$tAge)), !any(is.na(plot_df$time_after_treatment)))

# BASELINE-SUBTRACT, per cell type and clock, so the axis label is true and the
# panels are comparable across cell types. Without it the untreated medians are
# non-zero and differ by cell type - trivially on YuGene and mortality (< 0.8 and
# < 0.04 units) but by 4.20 for keratinocytes and -1.99 for melanocytes on the
# scaled difference model. On that row the raw heights therefore invert the
# ordering: keratinocytes look highest at 10 days when their elevation (+17.9) is
# below the fibroblasts' (+20.6). The median is subtracted, not the mean, to
# match the Wilcoxon tests and the median differences in tage_temporal_tests.csv.
# Every significance bracket is a WITHIN-cell-type comparison, so subtracting a
# per-cell-type constant leaves all of them unchanged.
baselines <- plot_df %>%
  dplyr::filter(time_after_treatment == "none") %>%
  dplyr::group_by(cell_type, model) %>%
  dplyr::summarise(base = median(tAge), .groups = "drop")
plot_df <- plot_df %>%
  dplyr::left_join(baselines, by = c("cell_type", "model")) %>%
  dplyr::mutate(tAge = tAge - base) %>%
  dplyr::select(-base)
cat("\nbaseline subtracted (median of each cell type's untreated samples):\n")
print(as.data.frame(baselines %>% tidyr::pivot_wider(names_from = model, values_from = base)),
      row.names = FALSE, digits = 3)

cat(sprintf("\nsamples per clock: %s\n",
            paste(table(plot_df$model), collapse = " / ")))

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
  # Magnitude ALONGSIDE the asterisk, matching 19_within_study_sample_figure.R.
  # 34 of the 54 contrasts sit at the 6v6 rank-test floor, so the asterisks are
  # saturated: they establish that the groups separate but cannot distinguish a
  # +0.41 mortality shift from a +31.2 YuGene one, both marked "**". The median
  # difference is the only quantity that can. Two decimals below 5 units so the
  # mortality row stays readable, one above.
  sub$stars <- sig_symbol(sub$p.adj)
  sub$label <- ifelse(abs(sub$diff) < 5,
                      sprintf("%+.2f %s", sub$diff, sub$stars),
                      sprintf("%+.1f %s", sub$diff, sub$stars))
  sub <- sub[sub$stars != "ns", ]
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
             model_label = factor(MODEL_LABELS[[model_name]], levels = MODEL_LABELS),
             cell_type = sub$cell_type, group1 = sub$timepoint_1, group2 = sub$timepoint_2,
             label = sub$label, diff = sub$diff, y.position = sub$y.position)
}
stat_df <- do.call(rbind, lapply(MODEL_LEVELS, make_stat_df))

# FOURTH COLUMN: the three cell types overlaid, on each row's shared y scale, so
# the cross-cell-type comparison the text makes ("keratinocytes show the smallest
# response") can be read directly instead of by eye across three panels. Median
# trajectories only - the distributions are in the per-cell-type columns.
CT_ORDER <- levels(plot_df$cell_type)
overlay <- plot_df %>% mutate(ct_orig = cell_type,
                              cell_type = factor("All three", levels = "All three"))
plot_df$ct_orig <- plot_df$cell_type
both_df <- bind_rows(plot_df, overlay) %>%
  mutate(cell_type = factor(as.character(cell_type), levels = c(CT_ORDER, "All three")),
         ct_orig = factor(as.character(ct_orig), levels = CT_ORDER))
CT_COLS <- setNames(c("#D55E00", "#0072B2", "#009E73"), CT_ORDER)

p <- ggplot(both_df, aes(x = time_after_treatment, y = tAge)) +
  geom_violin(data = ~ dplyr::filter(.x, cell_type != "All three"),
              aes(fill = time_after_treatment), alpha = 0.6, trim = FALSE) +
  scale_fill_discrete(name = "Days post-irradiation",
                      labels = c("none" = "None", "4_days" = "4 days",
                                 "10_days" = "10 days", "20_days" = "20 days")) +
  geom_jitter(data = ~ dplyr::filter(.x, cell_type != "All three"),
              width = 0.08, size = 1.2, alpha = 0.5, colour = "black") +
  stat_summary(data = ~ dplyr::filter(.x, cell_type == "All three"),
               aes(group = ct_orig, colour = ct_orig), fun = median,
               geom = "line", linewidth = 1.1) +
  stat_summary(data = ~ dplyr::filter(.x, cell_type == "All three"),
               aes(group = ct_orig, colour = ct_orig), fun = median,
               geom = "point", size = 2.6) +
  scale_colour_manual(values = CT_COLS, name = NULL) +
  # Median trajectory across timepoints. Median (not mean) to match the
  # Wilcoxon tests and the median differences reported in
  # tage_temporal_pairwise_all_timepoints.csv. group=1 connects across the
  # discrete x within each panel.
  stat_summary(data = ~ dplyr::filter(.x, cell_type != "All three"),
               aes(group = 1), fun = median, geom = "line",
               colour = "grey45", linewidth = 0.7, alpha = 0.65) +
  stat_summary(data = ~ dplyr::filter(.x, cell_type != "All three"),
               aes(group = 1), fun = median, geom = "point",
               colour = "grey45", size = 1.8, alpha = 0.75) +
  stat_pvalue_manual(stat_df, label = "label", xmin = "group1", xmax = "group2",
                      y.position = "y.position", tip.length = 0, bracket.size = 0.5, size = 3.6) +
  # Relabel for display only -- the underlying factor levels stay as
  # none/4_days/... because the bracket table matches against those exact
  # strings from tage_temporal_pairwise_all_timepoints.csv.
  scale_x_discrete(labels = c("none" = "None", "4_days" = "4 days",
                              "10_days" = "10 days", "20_days" = "20 days")) +
  facet_grid(model_label ~ cell_type, scales = "free_y") +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "plain", size = 18),
        strip.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_text(size = 18)) +
  labs(x = "Days Post-Irradiation",
       y = "treated tAge - untreated tAge, within each cell type")

ggsave(file.path(RERUN_DIR, "figure_temporal_universal_tage.png"), p,
       width = 17, height = 14.5, dpi = 300)
cat(sprintf("Saved -> %s\n", file.path(RERUN_DIR, "figure_temporal_universal_tage.png")))
