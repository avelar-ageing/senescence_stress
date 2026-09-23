# 24_figure_cancellation.R
#
# Limitation figure: the net tAge difference is a small residual of much larger
# opposing per-gene contributions.
#
# WHY THIS FIGURE EXISTS. The set-level sections rest on decomposing a condition's
# tAge difference into per-gene, and then per-set, contributions. That is exact
# arithmetic, but it invites a reading the data does not support - that a condition
# produces a coherent shift which the sets then partition. It does not. Summed
# positive and negative gene contributions are an order of magnitude larger than
# their difference, so the reported net is what survives near-complete
# cancellation, and roughly half of measured clock genes move against it.
#
# This bounds how much weight any set-level contribution can carry, and it is the
# reason the yardstick asks whether a set exceeds its expected SHARE rather than
# whether its contribution is non-zero.
#
# LEVEL. This figure is TRANSCRIPTOME-WIDE: the bars sum over all measured clock
# genes and the diamond is the condition's whole-transcriptome net. The same
# cancellation occurs WITHIN individual gene sets, at a median of six-fold - see
# withinset_cancellation.csv from exploratory/23 - which is the level the text
# quotes when discussing what a set's own contribution means. The two must not be
# conflated: the within-set figures are about fortyfold smaller in absolute units.
#
# Reads exploratory/23's output; makes no calculation of its own beyond layout.
# Output: figure_contribution_cancellation.png

source("R/config.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })

d <- read.csv(file.path(RERUN_DIR, "decomposition_diagnostics.csv"))
d$arm_label <- ifelse(d$arm == "cross_sectional", "Arrest conditions",
                      "Irradiation time course")
ORD <- c("CICQ", "SSCQ", "RS", "SIPS", "OIS",
         paste0(rep(c("Fibroblast", "Keratinocyte", "Melanocyte"), each = 3),
                "_", c("4_days", "10_days", "20_days")))
d$group <- factor(d$group, levels = rev(ORD[ORD %in% d$group]))
d$pretty <- gsub("_", " ", as.character(d$group))
d$pretty <- factor(d$pretty, levels = gsub("_", " ", levels(d$group)))

bars <- d %>%
  select(pretty, arm_label, sum_positive, sum_negative) %>%
  pivot_longer(c(sum_positive, sum_negative), names_to = "side", values_to = "value") %>%
  mutate(side = ifelse(side == "sum_positive", "Genes raising the score",
                       "Genes lowering the score"))

# SHARE, NOT FOLD-RATIO. "34x cancellation" is total movement over the net, so it
# is 1/share and explodes as the net approaches zero: keratinocytes at 20 days read
# 87x because their net is 0.19, not because they cancel harder. Total movement is
# in fact near-constant across groups (8-17 units), so the ratio mostly restates
# the net printed beside it. The share surviving is bounded, is what the text
# quotes ("about 5% of all the movement"), and matches 28_figure_withinset_*.R,
# whose header rejected the ratio for the same reasons.
lab <- d %>% mutate(surviving_pct = 100 * abs(net) / (abs(sum_positive) + abs(sum_negative)),
                    txt = sprintf("net %+.2f  (%.0f%% of the movement)",
                                  net, surviving_pct))
xr <- max(abs(c(d$sum_positive, d$sum_negative)))

p <- ggplot(bars, aes(x = value, y = pretty, fill = side)) +
  geom_col(width = 0.62) +
  geom_point(data = d, aes(x = net, y = pretty), inherit.aes = FALSE,
             shape = 18, size = 3.6, colour = "black") +
  geom_text(data = lab, aes(x = xr * 1.07, y = pretty, label = txt),
            inherit.aes = FALSE, hjust = 0, size = 4.5) +
  geom_vline(xintercept = 0, colour = "grey30") +
  facet_grid(arm_label ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Genes raising the score" = "#B2182B",
                               "Genes lowering the score" = "#2166AC"), name = NULL) +
  scale_x_continuous(limits = c(-xr * 1.05, xr * 2.15)) +
  theme_bw(base_size = 17) +
  theme(legend.position = "top",
        strip.text = element_text(face = "plain", size = 16),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  labs(x = "Summed contribution of all measured clock genes to the mortality-score difference (diamond = net)",
       y = NULL)

out <- file.path(RERUN_DIR, "figure_contribution_cancellation.png")
ggsave(out, p, width = 12.5, height = 9.5, dpi = 300)
cat(sprintf("Saved -> %s\n", out))

cat("\n== what the figure shows ==\n")
print(d %>% transmute(arm, group, sum_positive = round(sum_positive, 2),
                      sum_negative = round(sum_negative, 2), net = round(net, 3),
                      cancellation = round(cancellation_ratio),
                      pct_genes_with_net = round(100 * frac_genes_in_net_direction, 1),
                      genes_for_whole_net = n_genes_for_the_whole_net),
      row.names = FALSE)
