# 40_figure_signal_concentration.R
#
# What the gene-set decomposition can and cannot see, on statistics that need no
# rank cut. Both panels come from 42_coverage_threshold_free.py.
#
# The earlier version of this figure plotted "genes needed to reach half the net"
# and the gene-set coverage of those genes. Both were biased. The count's length is
# set by how small the residue is, not by how concentrated the movement is (SSCQ
# cancels 71-fold and gets 3 genes, RS 12-fold and gets 20), which manufactured a
# quiescence/senescence coverage difference that disappears on a fixed-size list.
# And ranking by |contribution| selects on |coefficient|, which is itself related to
# annotation: annotated clock genes carry 1.05x the median |coefficient| and are
# 35.2% of the top decile against 29.3% overall.
#
# (a) SPREAD. n_eff = (sum|c|)^2 / sum(c^2) over every measured clock gene: the
#     Genes needed for half the TOTAL MOVEMENT, sum|c|, taken largest first: 651 to 971
#     of about 9,100 measured. Half the NET took 3 to 20 genes, but that count shrinks
#     with cancellation and is not a concentration statistic. This is, and it replaces
#     any claim that a few genes carry the shift. n_eff is in the CSV alongside it.
#
# (b) COVERAGE. Share of each group's total |contribution| sitting in genes that are
#     in none of the 50 sets, against a null that permutes set membership within
#     |coefficient| deciles, so the coefficient relationship above is held fixed. Shown
#     as the ANNOTATED share, the part a set-level decomposition can see: 36-39%
#     observed against 32% expected, every group at the 20,000-draw floor. The sets
#     carry more of the movement than their weight predicts, and still see under 40%.
source("R/config.R")
suppressMessages({library(dplyr); library(ggplot2); library(patchwork)})

ROWS <- c("CICQ", "SSCQ", "RS", "SIPS", "OIS",
          "Fibroblast 4 days", "Fibroblast 10 days", "Fibroblast 20 days",
          "Keratinocyte 4 days", "Keratinocyte 10 days", "Keratinocyte 20 days",
          "Melanocyte 4 days", "Melanocyte 10 days", "Melanocyte 20 days")
d <- read.csv(file.path(RERUN_DIR, "coverage_threshold_free.csv")) %>%
  mutate(dataset = ifelse(dataset == "arrest", "Arrest conditions", "Irradiation time course"),
         label = gsub("_", " ", group),
         label = factor(label, levels = rev(ROWS)))
stopifnot(!any(is.na(d$label)), nrow(d) == length(ROWS))

# ---- (a) how much of the clock do the gene sets account for? -------------
# Part-to-whole across three measures of the same model, so the reader sees the gain
# from counting genes to weighting them to measuring movement, and that the
# unaccounted majority survives all three. Palette #3A6EA5/#C97B2B passes the six
# colour checks (lightness band, chroma floor, CVD and normal-vision separation,
# contrast) and avoids panel (b)'s red and green.
acc <- read.csv(file.path(RERUN_DIR, "clock_unaccounted_shares.csv"))
LV <- c("by gene count", "by clock weight", "by tAge movement")
acc <- acc %>%
  mutate(rowlab = c("by gene count\n3,071 of 10,487 features",
                    "by clock weight\nsum of |coefficient|",
                    "by tAge movement\nsum of |contribution|")[match(measure, LV)],
         rowlab = factor(rowlab, levels = rev(rowlab[match(LV, measure)])))
pa <- acc %>%
  tidyr::pivot_longer(c(pct_in_hallmark, pct_not_in_hallmark),
                      names_to = "part", values_to = "pct") %>%
  mutate(part = factor(ifelse(part == "pct_in_hallmark", "in at least one of the 50 sets",
                              "in none of the 50 sets"),
                       levels = c("in none of the 50 sets", "in at least one of the 50 sets"))) %>%
  ggplot(aes(x = rowlab, y = pct, fill = part)) +
  geom_col(width = 0.6, colour = "white", linewidth = 1.1) +   # 2px surface gap
  geom_text(aes(label = sprintf("%.0f%%", pct)), position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 4.1) +
  coord_flip() +
  scale_fill_manual(values = c("in none of the 50 sets" = "#C97B2B",
                               "in at least one of the 50 sets" = "#3A6EA5"), name = NULL,
                    breaks = c("in none of the 50 sets", "in at least one of the 50 sets")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid = element_blank(),
        panel.border = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10, hjust = 0)) +
  labs(x = NULL, y = NULL)

# ---- (b) how much of the movement is outside the gene sets? --------------
# the null is tight (sd 0.7-1.0), so it is drawn as a +/- 3 sd ribbon per group
pb <- ggplot(d %>% mutate(ratio = ann_pct / ann_null_mean,
                          band = 3 * null_sd / ann_null_mean), aes(x = label)) +
  geom_hline(yintercept = 1, colour = "grey30") +
  geom_linerange(aes(ymin = 1 - band, ymax = 1 + band),
                 linewidth = 2.6, colour = "grey78") +
  geom_point(aes(y = ratio), shape = 21, size = 3.4, fill = "#333333", colour = "white",
             stroke = 0.5) +
  geom_text(aes(y = ratio, label = sprintf("%.2f", ratio)), hjust = -0.4, size = 3) +
  coord_flip(ylim = c(0.90, 1.30)) +
  scale_y_continuous(breaks = c(0.9, 1.0, 1.1, 1.2, 1.3)) +
  facet_grid(dataset ~ ., scales = "free_y", space = "free_y") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 9), axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = "grey70"),
        strip.text.y = element_text(size = 9.5, angle = -90)) +
  labs(x = NULL,
       y = paste0("absolute movement in the 50 sets, relative to what\n",
                  "their weight predicts   (1 = as predicted)"))

p <- pa + pb + plot_layout(widths = c(1, 1.2)) +
  plot_annotation(tag_levels = "a")
ggsave(file.path(RERUN_DIR, "figure_signal_concentration.png"), p,
       width = 14.5, height = 5.0, dpi = 300)
cat("Saved -> figure_signal_concentration.png\n")
print(acc, row.names = FALSE)
print(d %>% transmute(dataset, group, n_measured,                       n_half = n_half_movement, ann_pct = round(ann_pct, 1),
                      ann_null = round(ann_null_mean, 1), z_ann = round(z_ann, 1), p_emp), row.names = FALSE)
