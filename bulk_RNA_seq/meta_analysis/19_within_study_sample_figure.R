# 19_within_study_sample_figure.R
#
# Per-SAMPLE version of the within-study effect figure, replacing the per-study
# version in 17_universal_tage_figures.R.
#
# WHY. The per-study figure plotted one point per study with the point SIZE
# encoding precision weight, and gave no legend for it, so the most prominent
# visual channel carried information the reader could not decode. It also hid the
# sample sizes and the within-group spread that the tests actually run on.
#
# WHAT IS PLOTTED. One point per sample. The x value is that sample's tAge minus
# the MEAN OF ITS OWN STUDY'S proliferating controls. Raw per-sample values are not
# comparable between studies - that is the whole reason for the within-study design
# (control medians span 87.8 units) - but deviations from a sample's own control
# mean are, and they are exactly what the within-study estimate aggregates. The
# Proliferating row shows the controls themselves, scattering around zero by
# construction, which gives the reader the within-study noise floor to judge the
# arrest effects against.
#
# Two colourings are produced, since they answer different questions:
#   _bycondition  colour = arrest condition. The main figure.
#   _bystudy      colour = study. Shows whether an effect rests on agreement
#                 between studies or on one dominant study; no legend, since 34
#                 studies cannot be read off a key, so it is a structure plot.
#
# All three clocks appear as rows with free x scales, because mortality units are
# roughly fortyfold smaller than chronological ones.
#
# Output: figure_within_study_samples_bycondition.png
#         figure_within_study_samples_bystudy.png

source("R/config.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

COND <- c("CICQ", "SSCQ", "RS", "SIPS", "OIS")
ROWS <- c("Proliferating", COND)
# Okabe-Ito, colourblind-safe
PAL <- c(CICQ = "#0072B2", SSCQ = "#56B4E9", RS = "#D55E00",
         SIPS = "#E69F00", OIS = "#CC79A7", Proliferating = "grey55")
sig_symbol <- function(p) ifelse(p < 0.001, "***",
                          ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

chron <- read.csv(file.path(RERUN_DIR, "tage_all_conditions.csv"))
mort  <- read.csv(file.path(RERUN_DIR, "mortality_tage.csv"))
long <- bind_rows(
  chron %>% transmute(external_id, study, condition,
                      clock = "Chronological, Scaled Difference",
                      value = scaled_diff_EN_tAge),
  chron %>% transmute(external_id, study, condition,
                      clock = "Chronological, YuGene",
                      value = yugene_diff_EN_tAge),
  mort  %>% transmute(external_id, study, condition,
                      clock = "Mortality", value = mortality_tAge))

# centre every sample on the mean of its OWN study's proliferating controls
ctrl <- long %>% filter(condition == "Proliferating") %>%
  group_by(clock, study) %>% summarise(ref = mean(value), .groups = "drop")
d <- long %>% inner_join(ctrl, by = c("clock", "study")) %>%
  mutate(dev = value - ref,
         condition = factor(condition, levels = rev(ROWS)),
         clock = factor(clock, levels = c("Chronological, Scaled Difference",
                                          "Chronological, YuGene", "Mortality")))
cat(sprintf("samples plotted per clock: %d of %d (studies lacking internal controls dropped)\n",
            nrow(d) / 3, nrow(long) / 3))

# combined within-study estimates, for the reference bars and labels
cw <- read.csv(file.path(RERUN_DIR, "condition_within_study.csv")) %>%
  filter(test == "condition_within_study_stratified") %>%
  transmute(condition, clock = ifelse(model == "scaled_diff",
              "Chronological, Scaled Difference", "Chronological, YuGene"),
            est = diff_within_study, p = p_perm_adj)
mw <- read.csv(file.path(RERUN_DIR, "mortality_within_study.csv")) %>%
  filter(test == "condition_within_study") %>%
  transmute(condition, clock = "Mortality", est = diff_within_study, p = p_perm_adj)
est <- bind_rows(cw, mw) %>%
  mutate(condition = factor(condition, levels = rev(ROWS)),
         clock = factor(clock, levels = levels(d$clock)),
         lab = ifelse(abs(est) < 5, sprintf("%+.2f %s", est, sig_symbol(p)),
                                    sprintf("%+.1f %s", est, sig_symbol(p))))
pos <- d %>% group_by(clock) %>%
  summarise(lo = min(dev), hi = max(dev), .groups = "drop") %>%
  mutate(span = hi - lo, xlab = hi + span * 0.06, xend = hi + span * 0.30)
est <- left_join(est, pos, by = "clock")

base <- function(colour_by) {
  ggplot(d, aes(x = dev, y = condition)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_jitter(aes(colour = .data[[colour_by]]), height = 0.22, size = 2.1,
                alpha = 0.75) +
    geom_segment(data = est, aes(x = est, xend = est,
                                 y = as.numeric(condition) - 0.34,
                                 yend = as.numeric(condition) + 0.34),
                 linewidth = 1.1, colour = "black", inherit.aes = FALSE) +
    geom_text(data = est, aes(x = xlab, y = condition, label = lab),
              size = 4.9, hjust = 0, inherit.aes = FALSE) +
    geom_blank(data = est, aes(x = xend, y = condition), inherit.aes = FALSE) +
    facet_wrap(~clock, ncol = 1, scales = "free_x") +
    theme_bw(base_size = 18) +
    theme(strip.text = element_text(face = "plain", size = 18),
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 17)) +
    labs(x = "Difference from the mean of the same study's proliferating controls (each clock's own units)",
         y = NULL)
}

p1 <- base("condition") +
  scale_colour_manual(values = PAL, guide = "none")
ggsave(file.path(RERUN_DIR, "figure_within_study_samples_bycondition.png"), p1,
       width = 12, height = 12.5, dpi = 300)

p2 <- base("study") +
  scale_colour_viridis_d(option = "turbo", guide = "none")
ggsave(file.path(RERUN_DIR, "figure_within_study_samples_bystudy.png"), p2,
       width = 12, height = 12.5, dpi = 300)

cat("Saved:\n  figure_within_study_samples_bycondition.png\n",
    " figure_within_study_samples_bystudy.png\n")
