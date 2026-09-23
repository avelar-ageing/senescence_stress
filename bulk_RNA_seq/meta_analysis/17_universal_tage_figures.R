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
# WHAT REPLACES IT. Two panels were produced; only Figure B is still in use.
#
#   Figure A - per-study effects. WITHDRAWN 2026-08-25, see the note at the
#   ggsave site below; superseded by meta_analysis/19. Originally described as:
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
# CLOCK ORDER (2026-09-01): YuGene, Scaled Difference, Mortality everywhere -
# figures, printed tables and the manuscript text all follow this order.
MODEL_LAB <- c(yugene_diff = "Chronological, YuGene",
               scaled_diff = "Chronological, Scaled Difference",
               mortality   = "Mortality")
sig_symbol <- function(p) ifelse(p < 0.001, "***",
                          ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

cw <- read.csv(file.path(RERUN_DIR, "condition_within_study.csv"))

# The mortality clock has no per-study rows saved, so derive them here from the
# per-sample values, using the same per-study difference of means and the same
# precision weight as the chronological arm.
mt <- read.csv(file.path(RERUN_DIR, "mortality_tage.csv"))
mort_per <- do.call(rbind, lapply(COND, function(cond) {
  s <- mt[mt$condition %in% c(cond, "Proliferating"), ]
  do.call(rbind, lapply(split(s, s$study), function(g) {
    x <- g$mortality_tAge[g$condition == cond]
    y <- g$mortality_tAge[g$condition == "Proliferating"]
    if (!length(x) || !length(y)) return(NULL)
    data.frame(test = "per_study_contrast", condition = cond, model = "mortality",
               study = g$study[1], n_test = length(x), n_control = length(y),
               diff = mean(x) - mean(y))
  }))
}))
mort_comb <- do.call(rbind, lapply(COND, function(cond) {
  d <- mort_per[mort_per$condition == cond, ]
  w <- d$n_test * d$n_control / (d$n_test + d$n_control)
  m <- read.csv(file.path(RERUN_DIR, "mortality_within_study.csv"))
  m <- m[m$test == "condition_within_study" & m$condition == cond, ]
  data.frame(test = "condition_within_study_stratified", condition = cond,
             model = "mortality", diff_within_study = sum(w * d$diff) / sum(w),
             p_perm_adj = m$p_perm_adj[1])
}))
cw <- dplyr::bind_rows(cw, mort_per, mort_comb)

per <- cw %>% filter(test == "per_study_contrast") %>%
  mutate(condition = factor(condition, levels = rev(COND)),
         model_label = factor(MODEL_LAB[model], levels = unname(MODEL_LAB)),
         weight = n_test * n_control / (n_test + n_control))
comb <- cw %>% filter(test == "condition_within_study_stratified") %>%
  mutate(condition = factor(condition, levels = rev(COND)),
         model_label = factor(MODEL_LAB[model], levels = unname(MODEL_LAB)),
         # 2 decimals where the clock's units are small (mortality), 1 where they are large
         lab = ifelse(abs(diff_within_study) < 5,
                      sprintf("%+.2f %s", diff_within_study, sig_symbol(p_perm_adj)),
                      sprintf("%+.1f %s", diff_within_study, sig_symbol(p_perm_adj))))

# ---- Figure A -------------------------------------------------------------
# labels placed per facet, because the mortality clock's units are ~40x smaller
lab_pos <- per %>% group_by(model_label) %>%
  summarise(xmax = max(diff, na.rm = TRUE), xmin = min(diff, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(span = xmax - xmin, xlab = xmax + span * 0.14, xhi = xmax + span * 0.40)
comb <- dplyr::left_join(comb, lab_pos, by = "model_label")
A <- ggplot(per, aes(x = diff, y = condition)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(aes(size = weight), colour = "grey30", alpha = 0.55) +
  geom_segment(data = comb,
               aes(x = diff_within_study, xend = diff_within_study,
                   y = as.numeric(condition) - 0.30,
                   yend = as.numeric(condition) + 0.30),
               linewidth = 1.1, colour = "black", inherit.aes = FALSE) +
  geom_text(data = comb, aes(x = xlab, y = condition, label = lab),
            size = 5.2, hjust = 0) +
  geom_blank(data = comb, aes(x = xhi, y = condition)) +
  facet_wrap(~model_label, ncol = 1, scales = "free_x") +
  scale_size_continuous(range = c(1.6, 7), guide = "none") +
  theme_bw(base_size = 20) +
  theme(strip.text = element_text(face = "plain", size = 20),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20)) +
  labs(x = "Difference from same-study proliferating controls (each clock's own units)",
       y = NULL)
# FIGURE A WITHDRAWN 2026-08-25 and no longer written. It encoded each study's
# precision weight as point SIZE with no legend for it, so the most prominent
# visual channel could not be decoded, and it hid the sample sizes and
# within-group spread the tests are run on. Replaced by the per-sample version in
# meta_analysis/19_within_study_sample_figure.R, which plots one point per sample
# centred on its own study's control mean. The code above is left in place because
# object A feeds nothing else and the per-study aggregation it performs documents
# the estimator; the ggsave is deliberately removed rather than commented, so a
# rerun cannot recreate the withdrawn file.
invisible(A)

# ---- Figure B -------------------------------------------------------------
d <- read.csv(file.path(RERUN_DIR, "tage_all_conditions.csv"))
ann <- read.csv(file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv"))
d$line <- ann$cell_line_resolved[match(d$external_id, ann$external_id)]
# ALL THREE CLOCKS (2026-09-01). This figure carried the two chronological models
# only, because it predates the mortality clock entering the pipeline -- the same
# omission that was fixed in meta_analysis/11 and 16. The text quotes mortality
# baselines for these strains (HCA2 +1.0 against IMR90 -0.5) and the +1.19 mortality
# foreskin-versus-lung contrast, so the figure has to show that clock too.
d$mortality_tAge <- mt$mortality_tAge[match(d$external_id, mt$external_id)]
stopifnot(!any(is.na(d$mortality_tAge)))
CLOCK_LAB <- c(yugene_diff_EN_tAge = "YuGene",
               scaled_diff_EN_tAge = "Scaled Difference",
               mortality_tAge      = "Mortality")
p <- d %>% filter(condition == "Proliferating") %>%
  tidyr::pivot_longer(c(scaled_diff_EN_tAge, yugene_diff_EN_tAge, mortality_tAge),
                      names_to = "model", values_to = "tAge") %>%
  mutate(model_label = factor(CLOCK_LAB[model], levels = unname(CLOCK_LAB)))
ord <- p %>% filter(model_label == "Scaled Difference") %>%
  group_by(line) %>% summarise(m = median(tAge), .groups = "drop") %>% arrange(m)

# n>=3 MARKED ON THE FIGURE (2026-08-31). The text quotes only strains with at least
# three proliferating controls -- "it ran from IMR90 at -15.3 to HCA2 at +47.1" in 2.1.5,
# and the +33.3 to +47.1 foreskin range in the Discussion -- while this figure plots every
# strain that has any controls at all. Without the distinction a reader sees foreskin 2DD
# above HCA2 and HDF 10-2 below the stated foreskin floor, both on one or two samples, and
# the figure looks like it contradicts the text. Each label now carries its n and the
# strains below the threshold are drawn in light grey.
# TISSUE AS COLOUR (2026-09-01). The strains are already on the y axis; adding
# tissue as colour shows both levels at once, which is the point: tissue separates
# the controls, and strains still differ within a tissue. Palette #0072B2/#D55E00 is
# the Okabe-Ito blue/vermillion pair, validated for CVD separation (protan dE 21.9,
# tritan 30.9, normal 31.2) and for contrast against the panel surface.
p$tissue <- ann$tissue_verified[match(p$external_id, ann$external_id)]
p$study  <- ann$study[match(p$external_id, ann$external_id)]
# ENCODING (2026-09-01): colour = study, shape = tissue. Colour carries the study
# because the question a reader has at a strain like IMR90 (32 samples, 10 studies)
# is whether its spread is one lab or many; tissue only has two levels, so it rides
# on shape, which stays legible at this point size. With ~29 studies colour cannot
# identify a particular study reliably - it shows how many contribute and whether
# they cluster, which is what the figure is for.
STUDIES <- sort(unique(p$study))
STUDY_PAL <- setNames(grDevices::hcl.colors(length(STUDIES), "Dark 3"), STUDIES)
n_line <- p %>% filter(model_label == "Scaled Difference") %>% count(line, name = "n_ctrl")
p <- dplyr::left_join(p, n_line, by = "line")
p$enough <- p$n_ctrl >= 3
lab <- setNames(sprintf("%s (n=%d)", n_line$line, n_line$n_ctrl), n_line$line)
p$line <- factor(p$line, levels = ord$line, labels = lab[as.character(ord$line)])

B <- ggplot(p, aes(x = tAge, y = line, colour = study, shape = tissue, alpha = enough)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(size = 2.3) +
  stat_summary(fun = median, geom = "point", shape = 124, size = 8,
               colour = "grey15", show.legend = FALSE) +
  scale_colour_manual(values = STUDY_PAL, name = NULL) +
  scale_shape_manual(values = c(Foreskin = 17, Lung = 16), name = NULL) +
  scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0.3),
                     breaks = c("TRUE", "FALSE"),
                     labels = c("at least 3 controls", "fewer than 3 controls"),
                     name = NULL) +
  facet_wrap(~model_label, ncol = 3, scales = "free_x") +
  theme_bw(base_size = 18) +
  theme(strip.text = element_text(face = "plain", size = 18),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.position = "bottom",
        legend.box = "vertical", legend.box.just = "left",
        legend.text = element_text(size = 13)) +
  guides(shape  = guide_legend(order = 1, nrow = 1,
                               override.aes = list(size = 4, alpha = 1, colour = "grey25")),
         alpha  = guide_legend(order = 2, nrow = 1, override.aes = list(size = 4, colour = "grey25")),
         colour = guide_legend(order = 3, nrow = 4, byrow = TRUE,
                               override.aes = list(size = 3, alpha = 1))) +
  labs(x = "tAge of untreated proliferating controls", y = NULL)

# ---- the same figure restricted to the strains the text actually quotes --------
# The version above plots every strain that has any proliferating controls, with the
# sub-threshold ones faded. This one drops them, so the figure contains exactly the
# n>=3 strains the text quotes (-15.3 IMR90 to +47.1 HCA2, and the +55.1 foreskin-
# versus-lung contrast). With the single-sample strains gone the alpha legend is
# redundant and tissue separation on the scaled-difference clock is total: all five
# lung strains below all three foreskin strains, no interleaving.
p3 <- p[p$n_ctrl >= 3, ]
p3$line <- droplevels(p3$line)
B3 <- ggplot(p3, aes(x = tAge, y = line, colour = study, shape = tissue)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(size = 2.7, alpha = 0.85) +
  stat_summary(fun = median, geom = "point", shape = 124, size = 8,
               colour = "grey15", show.legend = FALSE) +
  scale_colour_manual(values = STUDY_PAL, name = NULL) +
  scale_shape_manual(values = c(Foreskin = 17, Lung = 16), name = NULL) +
  facet_wrap(~model_label, ncol = 3, scales = "free_x") +
  theme_bw(base_size = 18) +
  theme(strip.text = element_text(face = "plain", size = 18),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.position = "bottom",
        legend.box = "vertical", legend.box.just = "left",
        legend.text = element_text(size = 13)) +
  guides(shape  = guide_legend(order = 1, nrow = 1, override.aes = list(size = 4, colour = "grey25")),
         colour = guide_legend(order = 2, nrow = 4, byrow = TRUE,
                               override.aes = list(size = 3, alpha = 1))) +
  labs(x = "tAge of untreated proliferating controls", y = NULL)
# One figure, two panels: every strain (a) and the n >= 3 subset (b). Previously
# two files showing the same quantity at two inclusion thresholds.
suppressMessages(library(patchwork))
BB <- (B + labs(tag = "a")) / (B3 + labs(tag = "b")) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 20))
ggsave(file.path(RERUN_DIR, "figure_proliferating_baseline_heterogeneity.png"), BB,
       width = 16, height = 22, dpi = 300, limitsize = FALSE)
cat("Saved -> figure_proliferating_baseline_heterogeneity.png (a: all strains, b: n >= 3)\n")
cat(sprintf("  n>=3 version: %d strains, %d samples\n",
            length(unique(p3$line)), nrow(p3) / length(unique(p3$model_label))))

cat("Saved:\n  figure_proliferating_baseline_heterogeneity.png\n",
    " (Figure A withdrawn 2026-08-25 - see the note above; use meta_analysis/19)\n")

# The retirement step that used to live here renamed the old pooled figure to
# SUPERSEDED_figure_universal_tage_differences.png. Both files were deleted on
# 2026-08-25, and meta_analysis/08 - which produced the original - now carries a
# SUPERSEDED banner, so there is nothing left to retire and re-creating a
# SUPERSEDED_ file would only reintroduce a stale figure to the output directory.
if (file.exists(file.path(RERUN_DIR, "figure_universal_tage_differences.png"))) {
  warning("figure_universal_tage_differences.png is back in rerun_outputs; ",
          "meta_analysis/08 is superseded and should not have been run")
}
