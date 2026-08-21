# 17_universal_tage_figures.R
#
# Replaces 08_universal_tage_figure.R for section 2.1.5.1.
#
# WHY THE OLD FIGURE WAS WITHDRAWN. It plotted every sample's raw tAge as a violin
# per condition, with significance from the pooled vs-Proliferating family. Both
# halves are now known to be wrong. The significance came from pooled contrasts
# that 13_condition_within_study.R showed to be confounded by control baseline
# (CICQ +47.6 pooled, +5.9 within study). More importantly the VIOLINS themselves
# carry the confound: a condition's spread is dominated by which studies and cell
# strains contribute to it, so re-labelling the stars would have left a figure
# whose visual impression still said CICQ was old. The distribution being drawn is
# not the quantity being tested.
#
# WHAT REPLACES IT. Two panels, matching what is now actually estimated.
#
#   Figure A - per-study effects (figure_universal_tage_within_study.png).
#   One point per study per condition, each the difference between that study's
#   arrested and its own proliferating samples, sized by precision weight; the
#   combined within-study estimate as a vertical line; BH-adjusted permutation p
#   annotated. This shows the data the test actually uses, and makes the
#   between-study spread visible instead of hiding it inside a violin.
#
#   Figure B - why the pooled contrast fails
#   (figure_proliferating_baseline_heterogeneity.png). tAge of the 91 untreated
#   proliferating controls alone, by cell strain, ordered by median. This is the
#   evidence for the design change and belongs in the paper alongside it.
#
# Reads only existing outputs: condition_within_study.csv, tage_all_conditions.csv,
# immortalisation_annotation_corrected.csv.

source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2)
})

COND <- c("CICQ", "SSCQ", "RS", "SIPS", "OIS")
MODEL_LAB <- c(scaled_diff = "Scaled Difference", yugene_diff = "YuGene")
sig_symbol <- function(p) ifelse(p < 0.001, "***",
                          ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

cw <- read.csv(file.path(RERUN_DIR, "condition_within_study.csv"))
per <- cw %>% filter(test == "per_study_contrast") %>%
  mutate(condition = factor(condition, levels = rev(COND)),
         model_label = MODEL_LAB[model],
         weight = n_test * n_control / (n_test + n_control))
comb <- cw %>% filter(test == "condition_within_study_stratified") %>%
  mutate(condition = factor(condition, levels = rev(COND)),
         model_label = MODEL_LAB[model],
         lab = sprintf("%+.1f %s", diff_within_study, sig_symbol(p_perm_adj)))

# ---- Figure A -------------------------------------------------------------
xr <- range(per$diff, na.rm = TRUE)
pad <- diff(xr) * 0.20
A <- ggplot(per, aes(x = diff, y = condition)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(aes(size = weight), colour = "grey30", alpha = 0.55) +
  geom_segment(data = comb,
               aes(x = diff_within_study, xend = diff_within_study,
                   y = as.numeric(condition) - 0.30,
                   yend = as.numeric(condition) + 0.30),
               linewidth = 1.1, colour = "black", inherit.aes = FALSE) +
  geom_text(data = comb, aes(x = xr[2] + pad * 0.62, y = condition, label = lab),
            size = 5.6, hjust = 0) +
  facet_wrap(~model_label, ncol = 1) +
  scale_size_continuous(range = c(1.6, 7), guide = "none") +
  scale_x_continuous(limits = c(xr[1] - pad * 0.15, xr[2] + pad * 1.9)) +
  theme_bw(base_size = 20) +
  theme(strip.text = element_text(face = "plain", size = 20),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20)) +
  labs(x = "tAge difference from same-study proliferating controls", y = NULL)
ggsave(file.path(RERUN_DIR, "figure_universal_tage_within_study.png"), A,
       width = 11, height = 10, dpi = 300)

# ---- Figure B -------------------------------------------------------------
d <- read.csv(file.path(RERUN_DIR, "tage_all_conditions.csv"))
ann <- read.csv(file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv"))
d$line <- ann$cell_line_resolved[match(d$external_id, ann$external_id)]
p <- d %>% filter(condition == "Proliferating") %>%
  tidyr::pivot_longer(c(scaled_diff_EN_tAge, yugene_diff_EN_tAge),
                      names_to = "model", values_to = "tAge") %>%
  mutate(model_label = ifelse(grepl("^scaled", model), "Scaled Difference", "YuGene"))
ord <- p %>% filter(model_label == "Scaled Difference") %>%
  group_by(line) %>% summarise(m = median(tAge), .groups = "drop") %>% arrange(m)
p$line <- factor(p$line, levels = ord$line)
B <- ggplot(p, aes(x = tAge, y = line)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(size = 2, alpha = 0.5, colour = "grey30") +
  stat_summary(fun = median, geom = "point", shape = 124, size = 8, colour = "black") +
  facet_wrap(~model_label, ncol = 2, scales = "free_x") +
  theme_bw(base_size = 18) +
  theme(strip.text = element_text(face = "plain", size = 18),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18)) +
  labs(x = "tAge of untreated proliferating controls", y = NULL)
ggsave(file.path(RERUN_DIR, "figure_proliferating_baseline_heterogeneity.png"), B,
       width = 13, height = 7, dpi = 300)

cat("Saved:\n  figure_universal_tage_within_study.png\n",
    " figure_proliferating_baseline_heterogeneity.png\n")

# retire the superseded figure so it cannot be picked up by mistake
old <- file.path(RERUN_DIR, "figure_universal_tage_differences.png")
if (file.exists(old)) {
  file.rename(old, file.path(RERUN_DIR, "SUPERSEDED_figure_universal_tage_differences.png"))
  cat("Renamed old pooled figure -> SUPERSEDED_...png\n")
}
