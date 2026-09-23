# 30_per_study_effects_figure.R
#
# One point per STUDY per condition: the per-study difference d_s that the
# within-study estimate aggregates.
#
# WHY THIS SCRIPT EXISTS. 2.1.5.1 claims the conditions are elevated "in the same
# direction in every contributing study". That is a claim about the sign of every
# d_s, and no existing figure shows it. figure_within_study_samples_bycondition.png
# plots one point per SAMPLE, so within-study scatter dominates and study-level
# agreement cannot be read off it; the _bystudy variant colours by study but has no
# legend by design, so it shows structure rather than which study went which way.
# This figure plots d_s itself, so the reader can count the signs.
#
# Per-study effects are recomputed here rather than read from
# condition_within_study.csv because that file stores them for the two
# chronological clocks only - mortality_within_study.csv has no per-study rows - and
# the figure should carry all three.
#
# WHAT IS PLOTTED. One point per SAMPLE: that sample's tAge minus the mean of its
# OWN study's proliferating controls. Raw tAge is not comparable between studies -
# control medians span 87.8 units, which is the whole reason for the within-study
# design - but deviations from a sample's own control mean are, and they are exactly
# what the within-study estimate aggregates. Colour is study, so the reader can see
# whether a condition's elevation rests on agreement across studies or on one
# dominant study. The Proliferating row is the controls themselves, scattering around
# zero by construction, which gives the within-study noise floor to judge the arrest
# rows against. The red diamond is the precision-weighted mean T quoted in the text,
# with w_s = n_t n_c/(n_t+n_c). One facet per clock, free x, because mortality is on
# a different scale.
#
# Output: rerun_outputs/figure_per_study_effects.png
#         rerun_outputs/per_study_effects.csv

source("R/config.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

d   <- read.csv(file.path(RERUN_DIR, "tage_all_conditions.csv"))
ann <- read.csv(file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv"))
mt  <- read.csv(file.path(RERUN_DIR, "mortality_tage.csv"))
d$study <- ann$study[match(d$external_id, ann$external_id)]
d$mortality_tAge <- mt$mortality_tAge[match(d$external_id, mt$external_id)]

CLOCKS <- c("YuGene" = "yugene_diff_EN_tAge",
            "Scaled difference" = "scaled_diff_EN_tAge",
            "Mortality" = "mortality_tAge")
CONDS <- c("RS", "SIPS", "OIS", "CICQ", "SSCQ")

rows <- list(); wrows <- list()
for (cl in names(CLOCKS)) {
  v <- CLOCKS[[cl]]
  for (st in unique(d$study)) {
    q  <- d[d$study == st, ]
    co <- q[[v]][q$condition == "Proliferating"]
    if (length(co) == 0) next
    base <- mean(co)
    for (cond in c("Proliferating", CONDS)) {
      te <- q[[v]][q$condition == cond]
      if (length(te) == 0) next
      rows[[length(rows) + 1]] <- data.frame(
        clock = cl, condition = cond, study = st,
        centred = te - base)
      if (cond != "Proliferating")
        wrows[[length(wrows) + 1]] <- data.frame(
          clock = cl, condition = cond, study = st,
          d_s = mean(te) - base,
          w = length(te) * length(co) / (length(te) + length(co)))
    }
  }
}
ps <- bind_rows(rows); pw <- bind_rows(wrows)
LEV <- rev(c(CONDS, "Proliferating"))
ps$condition <- factor(ps$condition, levels = LEV)
pw$condition <- factor(pw$condition, levels = LEV)
ps$clock <- factor(ps$clock, levels = names(CLOCKS))
pw$clock <- factor(pw$clock, levels = names(CLOCKS))

wm <- pw %>% group_by(clock, condition) %>%
  summarise(T = sum(w * d_s) / sum(w), n_studies = n(),
            n_pos = sum(d_s > 0), .groups = "drop")

cat("== per-study sign counts (what the text claims) ==\n")
print(as.data.frame(wm %>% mutate(agree = sprintf("%d/%d", n_pos, n_studies)) %>%
        select(clock, condition, T, agree)), row.names = FALSE, digits = 3)
cat(sprintf("\nsamples plotted: %d   studies: %d\n", nrow(ps), length(unique(ps$study))))

studies <- sort(unique(ps$study))
pal <- setNames(grDevices::hcl.colors(length(studies), "Dark 3"), studies)

p <- ggplot(ps, aes(x = centred, y = condition)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(aes(colour = study), position = position_jitter(height = 0.22, width = 0),
             size = 1.9, alpha = 0.85) +
  geom_point(data = wm, aes(x = T, y = condition), shape = 18, size = 5.5,
             colour = "firebrick", inherit.aes = FALSE) +
  scale_colour_manual(values = pal, name = NULL) +
  facet_wrap(~clock, ncol = 1, scales = "free_x") +
  guides(colour = guide_legend(ncol = 6, override.aes = list(size = 3, alpha = 1))) +
  theme_bw(base_size = 15) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), legend.position = "bottom",
        legend.text = element_text(size = 9), legend.key.height = unit(9, "pt")) +
  labs(x = "sample tAge - proliferating mean tAge", y = NULL)
ggsave(file.path(RERUN_DIR, "figure_per_study_effects.png"), p,
       width = 11, height = 13, dpi = 300)
write.csv(ps, file.path(RERUN_DIR, "per_study_effects.csv"), row.names = FALSE)
cat(sprintf("\nSaved -> %s\n         %s\n",
    file.path(RERUN_DIR, "figure_per_study_effects.png"),
    file.path(RERUN_DIR, "per_study_effects.csv")))
