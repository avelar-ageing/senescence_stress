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
wilcox_vs_none <- read.csv(file.path(RERUN_DIR, "tage_temporal_wilcoxon_vs_none.csv"))

plot_df <- tage_temporal %>%
  dplyr::select(cell_type, time_after_treatment, scaled_diff_EN_tAge, yugene_diff_EN_tAge) %>%
  tidyr::pivot_longer(cols = c(scaled_diff_EN_tAge, yugene_diff_EN_tAge),
                       names_to = "model", values_to = "tAge") %>%
  mutate(model = ifelse(model == "scaled_diff_EN_tAge", "scaled_diff", "yugene_diff"),
         model_label = ifelse(model == "scaled_diff", "Scaled difference EN model", "YuGene EN model"))

sig_symbol <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

# Brackets: each timepoint vs 'none', per cell type, per model -- using the
# already-computed 18-test (3 cell types x 3 timepoints x 2 models) BH family
# from tage_temporal_wilcoxon_vs_none.csv. Non-significant comparisons
# dropped entirely.
make_stat_df <- function(model_name) {
  sub <- wilcox_vs_none[wilcox_vs_none$model == model_name, ]
  sub$label <- sig_symbol(sub$p.adj)
  sub <- sub[sub$label != "ns", ]
  sub$cell_type <- factor(sub$cell_type, levels = CT_LEVELS)
  sub$timepoint <- factor(sub$timepoint, levels = TIME_LEVELS)
  sub <- sub[order(sub$cell_type, sub$timepoint), ]

  y_max <- plot_df %>% filter(model == model_name) %>% group_by(cell_type) %>%
    summarise(y_max = max(tAge), y_range = diff(range(tAge)), .groups = "drop")
  # left_join (not merge/base R, which silently re-sorts by the join key and
  # so scrambled the none-vs-4d/10d/20d stacking order) preserves sub's
  # existing cell_type/timepoint order.
  sub <- dplyr::left_join(sub, y_max, by = "cell_type")
  sub <- sub %>% group_by(cell_type) %>%
    mutate(y.position = y_max + y_range * (0.25 + 0.18 * row_number())) %>% ungroup()

  data.frame(model = model_name,
             model_label = ifelse(model_name == "scaled_diff", "Scaled difference EN model", "YuGene EN model"),
             cell_type = sub$cell_type, group1 = "none", group2 = sub$timepoint,
             label = sub$label, y.position = sub$y.position)
}
stat_df <- rbind(make_stat_df("scaled_diff"), make_stat_df("yugene_diff"))

p <- ggplot(plot_df, aes(x = time_after_treatment, y = tAge, fill = time_after_treatment)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_jitter(width = 0.08, size = 1.2, alpha = 0.5, colour = "black") +
  stat_pvalue_manual(stat_df, label = "label", xmin = "group1", xmax = "group2",
                      y.position = "y.position", tip.length = 0, bracket.size = 0.5, size = 5) +
  facet_grid(model_label ~ cell_type, scales = "free_y") +
  theme_bw(base_size = 18) +
  theme(legend.position = "none",
        strip.text = element_text(face = "plain", size = 18),
        strip.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_text(size = 18)) +
  labs(x = NULL, y = "tAge")

ggsave(file.path(RERUN_DIR, "figure_temporal_universal_tage.png"), p, width = 13, height = 10, dpi = 300)
cat(sprintf("Saved -> %s\n", file.path(RERUN_DIR, "figure_temporal_universal_tage.png")))
