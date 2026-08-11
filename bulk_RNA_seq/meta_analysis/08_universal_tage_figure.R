# 08_universal_tage_figure.R
#
# Figure for section 2.1.5.1 (Differences in Universal Transcriptomic Age):
# tAge spread by condition, both EN models, with BH-adjusted Wilcoxon
# significance vs Proliferating annotated (family = 5 conditions x 2 models,
# matching tage_wilcoxon_vs_proliferating.csv / 05_tage_all_conditions.R).
# Reads only the already-computed outputs of 05_tage_all_conditions.R --
# does not re-run tAge prediction.

source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
})

tage_result <- read.csv(file.path(RERUN_DIR, "tage_all_conditions.csv"))
tage_result$condition <- factor(tage_result$condition,
                                 levels = c("Proliferating", "CICQ", "SSCQ", "RS", "SIPS", "OIS"))
wilcox_vs_prolif <- read.csv(file.path(RERUN_DIR, "tage_wilcoxon_vs_proliferating.csv"))

plot_df <- tage_result %>%
  dplyr::select(condition, scaled_diff_EN_tAge, yugene_diff_EN_tAge) %>%
  tidyr::pivot_longer(cols = c(scaled_diff_EN_tAge, yugene_diff_EN_tAge),
                       names_to = "model", values_to = "tAge") %>%
  mutate(model = ifelse(model == "scaled_diff_EN_tAge", "scaled_diff", "yugene_diff"),
         model_label = ifelse(model == "scaled_diff", "Scaled difference EN model", "YuGene EN model"))

sig_symbol <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

# Build per-facet significance annotations vs Proliferating (group1), placed
# well above each facet's max value so brackets clear the violins, and
# dropping non-significant comparisons entirely (no line, no label).
make_stat_df <- function(model_name) {
  sub <- wilcox_vs_prolif[wilcox_vs_prolif$model == paste0(model_name, "_EN_tAge") |
                            wilcox_vs_prolif$model == model_name, ]
  if (nrow(sub) == 0) sub <- wilcox_vs_prolif[grepl(model_name, wilcox_vs_prolif$model), ]
  y_max <- max(plot_df$tAge[plot_df$model == model_name])
  y_range <- diff(range(plot_df$tAge[plot_df$model == model_name]))
  data.frame(
    model = model_name,
    model_label = ifelse(model_name == "scaled_diff", "Scaled difference EN model", "YuGene EN model"),
    group1 = "Proliferating", group2 = sub$condition,
    p.adj = sub$p.adj, label = sig_symbol(sub$p.adj),
    y.position = y_max + y_range * (0.25 + 0.16 * seq_len(nrow(sub)))
  )
}
stat_df <- rbind(make_stat_df("scaled_diff"), make_stat_df("yugene_diff"))
stat_df <- stat_df[stat_df$label != "ns", ]

p <- ggplot(plot_df, aes(x = condition, y = tAge, fill = condition)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_jitter(width = 0.08, size = 0.6, alpha = 0.35, colour = "black") +
  stat_pvalue_manual(stat_df, label = "label", xmin = "group1", xmax = "group2",
                      y.position = "y.position", tip.length = 0, bracket.size = 0.5, size = 6) +
  facet_wrap(~model_label, ncol = 1, scales = "free_y") +
  theme_bw(base_size = 20) +
  theme(legend.position = "none",
        strip.text = element_text(face = "plain", size = 20),
        strip.background = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20)) +
  labs(x = NULL, y = "tAge")

ggsave(file.path(RERUN_DIR, "figure_universal_tage_differences.png"), p, width = 9, height = 11, dpi = 300)
cat(sprintf("Saved -> %s\n", file.path(RERUN_DIR, "figure_universal_tage_differences.png")))
