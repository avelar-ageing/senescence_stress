# 25_gse175533_figures.R
#
# Figures for the GSE175533 immortalisation test (meta_analysis/20-22).
#
# TWO FIGURES, because the result has two halves that a single panel conflates.
#
#   figure_gse175533_htert.png
#     tAge against population doublings, one row per clock, one colour per arm,
#     with a least-squares fit and its slope per arm. This is the SLOPE half: the
#     parental arm climbs, the hTERT arm does not, which is the authors' own
#     finding reproduced on our pipeline. It also shows the LEVEL half in the same
#     picture - the hTERT points sit below the parental ones - and makes visible
#     why the two arms cannot be matched on doublings while both divide, since the
#     hTERT course begins at PD 46, the parental replicative limit. Zero on the y
#     axis is the mean of the parental dividing samples, so the parental cloud
#     starting near zero is construction, not a result.
#
#   figure_immortalisation_reversal.png
#     The same contrast estimated two ways: across our 34 studies, where
#     immortalisation is perfectly nested within study, and within this one
#     laboratory and strain. This is the CONFOUND half. The sign flips, which is
#     the point, so the panel is drawn on a signed axis with zero marked rather
#     than as two bars whose direction has to be read off the numbers.
#
# WHY THE ESTIMATES SHOWN ARE TIMEPOINT-LEVEL. Three replicate libraries from one
# culture at one doubling are not three observations. The reversal figure therefore
# quotes the timepoint-level contrast from script 22 (6 hTERT timepoints against 5
# parental), whose p = 0.0043 is the floor attainable at that size. The trajectory
# figure plots all 48 libraries, because there the point is the shape of the
# trajectory rather than a p-value.
#
# Output: figure_gse175533_htert.png, figure_immortalisation_reversal.png

source("R/config.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

DIR <- file.path(RERUN_DIR, "gse175533")
# CLOCK ORDER (2026-09-01): YuGene, Scaled Difference, Mortality everywhere -
# figures, printed tables and the manuscript text all follow this order.
CLOCKS <- c(chrono_yugene_diff = "Chronological, YuGene",
            chrono_scaled_diff = "Chronological, Scaled Difference",
            mortality = "Mortality")
ARM_PAL <- c(parental = "#0072B2", hTERT = "#D55E00")   # Okabe-Ito
ARM_LAB <- c(parental = "WI-38 parental", hTERT = "WI-38 hTERT")

fmt_p <- function(p) ifelse(p < 1e-4, "p < 1e-4", sprintf("p = %.3g", p))

d <- read.csv(file.path(DIR, "gse175533_tage.csv"))
ct <- read.csv(file.path(DIR, "gse175533_contrasts.csv"))

long <- d %>%
  tidyr::pivot_longer(all_of(names(CLOCKS)), names_to = "clock", values_to = "tAge") %>%
  mutate(clock = factor(CLOCKS[clock], levels = CLOCKS),
         arm = factor(arm, levels = names(ARM_PAL)))

# ---- figure 1: trajectory ------------------------------------------------
# slopes over population doublings, from script 22's D1 (recomputed here only to
# place the annotation; the values are asserted against the CSV below)
sl <- ct %>% filter(contrast == "D1_within_arm_slope", x == "population_doublings") %>%
  transmute(clock = factor(CLOCKS[clock], levels = CLOCKS), arm, slope, p)
chk <- long %>% group_by(clock, arm) %>%
  summarise(s = coef(lm(tAge ~ population_doublings))[2], .groups = "drop")
stopifnot(all(abs(chk$s - sl$slope[match(paste(chk$clock, chk$arm),
                                         paste(sl$clock, sl$arm))]) < 1e-8))
cat("slope annotations verified against gse175533_contrasts.csv\n")

lab <- sl %>%
  mutate(txt = sprintf("%s: %+.3f per doubling (%s)", ARM_LAB[as.character(arm)],
                       slope, fmt_p(p))) %>%
  group_by(clock) %>%
  summarise(txt = paste(txt, collapse = "\n"), .groups = "drop") %>%
  left_join(long %>% group_by(clock) %>%
              summarise(y = max(tAge) + 0.16 * diff(range(tAge)),
                        ytop = max(tAge) + 0.42 * diff(range(tAge)), .groups = "drop"),
            by = "clock")

p1 <- ggplot(long, aes(population_doublings, tAge, colour = arm, fill = arm)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_vline(xintercept = 46, linetype = "dotted", colour = "grey35") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.14,
              linewidth = 1.1) +
  geom_point(size = 2.6, alpha = 0.85) +
  geom_text(data = lab, aes(x = 20, y = y, label = txt), hjust = 0, vjust = 0,
            size = 4.6, colour = "black", inherit.aes = FALSE, lineheight = 1.15) +
  geom_blank(data = lab, aes(x = 20, y = ytop), inherit.aes = FALSE) +
  facet_wrap(~clock, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = ARM_PAL, labels = ARM_LAB, name = NULL) +
  scale_fill_manual(values = ARM_PAL, labels = ARM_LAB, name = NULL) +
  theme_bw(base_size = 18) +
  theme(legend.position = "top",
        strip.text = element_text(size = 17), strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 15), axis.title = element_text(size = 17)) +
  labs(x = paste("Population doublings\ndotted line: PD 46, the only doubling",
                 "present in both arms"),
       y = "tAge, difference from the parental dividing mean (each clock's own units)")


# ---- figure 2: the sign reversal ----------------------------------------
ours_chr <- read.csv(file.path(RERUN_DIR, "immortalisation_contrasts.csv")) %>%
  filter(test == "A_naive_within_condition_CONFOUNDED", stratum == "Proliferating") %>%
  transmute(clock = ifelse(model == "scaled_diff", "chrono_scaled_diff",
                                                   "chrono_yugene_diff"),
            est = diff, p = p)
ours_mort <- read.csv(file.path(RERUN_DIR, "mortality_within_study.csv")) %>%
  filter(test == "immortalised_vs_primary", stratum == "Proliferating") %>%
  transmute(clock = "mortality", est = diff, p = p)
here <- ct %>% filter(contrast == "A_dividing_vs_dividing_bytimepoint") %>%
  transmute(clock, est = diff, p)

cmp <- bind_rows(
  bind_rows(ours_chr, ours_mort) %>%
    mutate(design = "Across 34 studies\n(immortalisation nested within study)"),
  here %>% mutate(design = "Within one laboratory and strain\n(GSE175533, timepoint level)")) %>%
  mutate(clock = factor(CLOCKS[clock], levels = rev(CLOCKS)),
         design = factor(design, levels = c(
           "Across 34 studies\n(immortalisation nested within study)",
           "Within one laboratory and strain\n(GSE175533, timepoint level)")),
         lab = sprintf("%+.2f, %s", est, fmt_p(p)))

p2 <- ggplot(cmp, aes(est, clock, colour = design)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_line(aes(group = clock), colour = "grey70", linewidth = 0.9) +
  geom_point(size = 5) +
  geom_text(aes(label = lab), vjust = -1.4, size = 4.4, show.legend = FALSE) +
  facet_wrap(~clock, ncol = 1, scales = "free", strip.position = "left") +
  scale_colour_manual(values = c("grey45", "#D55E00"), name = NULL) +
  scale_y_discrete(labels = NULL, breaks = NULL) +
  # the value labels sit above their points, so the x range needs room at both
  # ends or they are clipped by the panel edge
  scale_x_continuous(expand = expansion(mult = 0.32)) +
  theme_bw(base_size = 17) +
  theme(legend.position = "top", legend.text = element_text(size = 13),
        strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 15),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = 16)) +
  labs(x = "Immortalised minus primary, proliferating cells (each clock's own units)",
       y = NULL)

# One figure, two panels: the GSE175533 hTERT contrast (a) and the cross-study
# reversal it sits inside (b). They are read together and were two files.
suppressMessages(library(patchwork))
PP <- (p1 / p2) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 18))
ggsave(file.path(RERUN_DIR, "figure_immortalisation_htert.png"), PP,
       width = 13, height = 16, dpi = 300, limitsize = FALSE)
cat("Saved -> figure_immortalisation_htert.png (a: GSE175533, b: cross-study reversal)\n")

cat("Saved:\n  figure_gse175533_htert.png\n  figure_immortalisation_reversal.png\n")
print(cmp %>% select(clock, design, est, p) %>% arrange(clock))
