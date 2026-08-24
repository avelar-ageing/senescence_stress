# 08_universal_tage_figure.R
#
# Figure for section 2.1.5.1 (Differences in Universal Transcriptomic Age):
# tAge spread by condition, both EN models, with BH-adjusted Wilcoxon
# significance annotated for (a) each condition vs Proliferating (family = 5
# conditions x 2 models, tage_wilcoxon_vs_proliferating.csv) and (b) SSCQ vs
# RS/SIPS/OIS (from the full 15-pair x 2-model family,
# tage_pairwise_all_conditions.csv) -- the pairwise result the text discusses
# (SSCQ significantly lower than SIPS/RS/OIS on YuGene, SIPS only on scaled).
# Reads only the already-computed outputs of 05_tage_all_conditions.R --
# does not re-run tAge prediction.

source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
})

COND_LEVELS <- c("Proliferating", "CICQ", "SSCQ", "RS", "SIPS", "OIS")
tage_result <- read.csv(file.path(RERUN_DIR, "tage_all_conditions.csv"))
tage_result$condition <- factor(tage_result$condition, levels = COND_LEVELS)
wilcox_vs_prolif <- read.csv(file.path(RERUN_DIR, "tage_wilcoxon_vs_proliferating.csv"))
pairwise <- read.csv(file.path(RERUN_DIR, "tage_pairwise_all_conditions.csv"))

plot_df <- tage_result %>%
  dplyr::select(condition, scaled_diff_EN_tAge, yugene_diff_EN_tAge) %>%
  tidyr::pivot_longer(cols = c(scaled_diff_EN_tAge, yugene_diff_EN_tAge),
                       names_to = "model", values_to = "tAge") %>%
  mutate(model = ifelse(model == "scaled_diff_EN_tAge", "scaled_diff", "yugene_diff"),
         model_label = ifelse(model == "scaled_diff", "Scaled Difference", "YuGene"))

sig_symbol <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
pos_of <- function(x) match(x, COND_LEVELS)

# Build per-facet significance annotations: each condition vs Proliferating,
# plus SSCQ vs RS/SIPS/OIS. Non-significant comparisons are dropped entirely
# (no line, no label). Brackets are stacked by span (number of x-categories
# crossed) so short ones sit low and long ones sit high, minimizing overlap.
make_stat_df <- function(model_name) {
  sub_vp <- wilcox_vs_prolif[wilcox_vs_prolif$model == paste0(model_name, "_EN_tAge") |
                               wilcox_vs_prolif$model == model_name, ]
  if (nrow(sub_vp) == 0) sub_vp <- wilcox_vs_prolif[grepl(model_name, wilcox_vs_prolif$model), ]
  vp <- data.frame(group1 = "Proliferating", group2 = sub_vp$condition, p.adj = sub_vp$p.adj)

  sub_pw <- pairwise[pairwise$model == model_name &
                       ((pairwise$condition_1 == "SSCQ" & pairwise$condition_2 %in% c("RS", "SIPS", "OIS")) |
                        (pairwise$condition_2 == "SSCQ" & pairwise$condition_1 %in% c("RS", "SIPS", "OIS"))), ]
  pw <- data.frame(group1 = "SSCQ",
                    group2 = ifelse(sub_pw$condition_1 == "SSCQ", sub_pw$condition_2, sub_pw$condition_1),
                    p.adj = sub_pw$p.adj)

  all_cmp <- rbind(vp, pw)
  all_cmp$label <- sig_symbol(all_cmp$p.adj)
  all_cmp <- all_cmp[all_cmp$label != "ns", ]
  all_cmp$span <- abs(pos_of(all_cmp$group2) - pos_of(all_cmp$group1))
  all_cmp <- all_cmp[order(all_cmp$span), ]

  y_max <- max(plot_df$tAge[plot_df$model == model_name])
  y_range <- diff(range(plot_df$tAge[plot_df$model == model_name]))
  all_cmp$model <- model_name
  all_cmp$model_label <- ifelse(model_name == "scaled_diff", "Scaled Difference", "YuGene")
  all_cmp$y.position <- y_max + y_range * (0.25 + 0.16 * seq_len(nrow(all_cmp)))
  all_cmp
}
stat_df <- rbind(make_stat_df("scaled_diff"), make_stat_df("yugene_diff"))

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
