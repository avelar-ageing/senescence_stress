# 26_mortality_pathway_heatmaps.R
#
# Gene-set contribution heatmaps on the MORTALITY clock, replacing the
# chronological-clock versions in exploratory/08, 09, 10 and 19.
#
# WHY THIS REPLACES THEM. Sections 2.1.5.2 and 2.2.4 are set-level and report the
# mortality clock only. That is not a presentational choice: the chronological
# clock is a sparse elastic net (1,839 of 10,487 features non-zero), so a Hallmark
# set retains ~20 weighted genes and the allocation among co-expressed genes is
# close to arbitrary, whereas the mortality clock is pure ridge with every feature
# non-zero, giving ~125 weighted genes per set. Decomposing a prediction by gene
# set is defensible on the second and marginal on the first. Keeping figures built
# on the chronological decomposition alongside a mortality-only text invited the
# reader to compare panels that were never the evidence.
#
# NO REPRESENTATION GATE. Scripts 08-10 admitted a set only if its five
# largest-contributing genes carried <= 65% of its total contribution under BOTH
# chronological models, which kept 17 of 50. That threshold was calibrated on the
# sparse clock, where the median top-5 share is 64-73%. On the ridge clock the
# median is 25% and 49 of 50 sets pass, so the gate discriminates nothing here and
# only silently drops sets. All 50 are shown.
#
# EFFECT SIZE. Fill is contrib_diff_within_study = the set's summed
# coefficient-weighted contribution in the arrested group minus the control group,
# in the clock's own units - not Cohen's d, which divides by a within-group SD that
# differs several-fold between the meta-analysis and the tight 6-sample temporal
# design and so is not comparable across the two panels.
#
# STARS AND THEIR FLOOR. BH-adjusted within each analysis as run in
# 22_mortality_pathway_decomposition.py: 250 tests for the meta-analysis (50 sets x
# 5 conditions) and 450 for the temporal arm (50 x 3 cell types x 3 timepoints).
# Every temporal comparison is 6 against 6, so the smallest attainable two-sided
# Wilcoxon p is 2/choose(12,6) = 0.00216 and nothing there can reach padj < 0.001;
# "***" is impossible in the temporal panel by construction, and its absence is a
# statement about n, not about effect size.
#
# Output: figure_pathway_heatmap_mortality_meta.png
#         figure_temporal_pathway_heatmap_mortality.png

source("R/config.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

COND_LEVELS <- c("CICQ", "SSCQ", "RS", "SIPS", "OIS")
CT_LEVELS <- c("Fibroblast", "Keratinocyte", "Melanocyte")
TIME_LEVELS <- c("4_days", "10_days", "20_days")

df <- read.csv(file.path(RERUN_DIR, "mortality_partial_tage_ALL.csv"),
               check.names = FALSE)
sig_symbol <- function(p) ifelse(p < 0.001, "***",
                          ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))

# rows ordered by the set's largest absolute contribution anywhere in the panel,
# so the sets that carry the signal are adjacent rather than scattered
# alphabetically; the order is stated in the caption so it is not mistaken for
# clustering
order_rows <- function(d) {
  o <- d %>% group_by(pathway_short) %>%
    summarise(m = max(abs(contrib_diff_within_study)), .groups = "drop") %>%
    arrange(m)
  factor(d$pathway_short, levels = o$pathway_short)
}

heat <- function(d, xvar, facet, file, width, height, xlab_angle = 0) {
  cap <- as.numeric(quantile(abs(d$contrib_diff_within_study), 0.97, na.rm = TRUE))
  p <- ggplot(d, aes(x = .data[[xvar]], y = pathway_short,
                     fill = contrib_diff_within_study)) +
    geom_tile(colour = "grey80") +
    geom_text(aes(label = stars), size = 3.2) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-cap, cap),
                         oob = scales::squish, name = "Mortality\ntAge units") +
    theme_minimal(base_size = 16) +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 12, angle = xlab_angle,
                                     hjust = if (xlab_angle) 1 else 0.5),
          axis.title = element_blank(),
          strip.text = element_text(size = 13),
          panel.grid = element_blank(),
          panel.spacing.x = unit(0.25, "lines"))
  if (!is.null(facet)) p <- p + facet_grid(cols = vars(.data[[facet]]))
  # panels are combined below into one figure; no per-panel file is written
  cat(sprintf("  panel %-52s (%d sets, fill capped at +/-%.3f)\n", file,
              length(unique(d$pathway_short)), cap))
  invisible(p)
}

# ---- meta-analysis: 50 sets x 5 conditions ------------------------------
m <- df %>% filter(analysis == "meta_analysis") %>%
  mutate(pathway_short = gsub("^HALLMARK ", "", pathway),
         condition = factor(label, levels = COND_LEVELS),
         stars = sig_symbol(p_adj))
m$pathway_short <- order_rows(m)
stopifnot(nrow(m) == 250)
pm <- heat(m, "condition", NULL, "figure_pathway_heatmap_mortality_meta.png", 8.5, 12)

# ---- temporal: 50 sets x 3 cell types x 3 timepoints --------------------
t <- df %>% filter(analysis == "temporal_bytimepoint") %>%
  mutate(pathway_short = gsub("^HALLMARK ", "", pathway),
         cell_type = factor(cell_type, levels = CT_LEVELS),
         timepoint = factor(timepoint, levels = TIME_LEVELS,
                            labels = c("4d", "10d", "20d")),
         stars = sig_symbol(p_adj))
t$pathway_short <- order_rows(t)
stopifnot(nrow(t) == 450)
pt <- heat(t, "timepoint", "cell_type", "figure_temporal_pathway_heatmap_mortality.png",
     11, 12)

# One figure, two panels: the arrest conditions (a) and the time course (b). The
# panels share a fill scale and a row order, so they are read against each other;
# they were previously two files.
suppressMessages(library(patchwork))
HH <- (pm + labs(subtitle = "Arrest conditions")) +
      (pt + labs(subtitle = "Irradiation time course")) +
      plot_layout(widths = c(5, 9), guides = "collect") +
      plot_annotation(tag_levels = "a") &
      theme(plot.tag = element_text(face = "bold", size = 18))
ggsave(file.path(RERUN_DIR, "figure_pathway_heatmap_mortality.png"), HH,
       width = 19, height = 13, dpi = 300, limitsize = FALSE)
cat("Saved -> figure_pathway_heatmap_mortality.png (a: arrest conditions, b: time course)\n")

cat("\nsets reaching padj < 0.05 at least once:\n")
cat(sprintf("  meta-analysis : %d of 50\n",
            n_distinct(m$pathway[m$p_adj < 0.05])))
cat(sprintf("  temporal      : %d of 50\n",
            n_distinct(t$pathway[t$p_adj < 0.05])))
