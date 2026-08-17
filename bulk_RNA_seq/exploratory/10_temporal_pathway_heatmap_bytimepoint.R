# 10_temporal_pathway_heatmap_bytimepoint.R
#
# Per-timepoint version of 09_temporal_pathway_heatmap_both_models.R: instead
# of pooling all irradiated samples per cell type, each timepoint (4/10/20
# days vs that cell type's own untreated baseline, 6 vs 6) gets its own
# column pair. 3 cell types x 3 timepoints x 2 models = 18 columns, all 50
# pathways. Shows whether a pathway's contribution to tAge is stable across
# the time course or specific to one timepoint -- which the pooled figure
# averages away.
#
# Same geom_tile design as 08/09 (alphabetical rows, no clustering).
#
# P-VALUE FLOOR (important when reading the stars): each comparison is 6 vs 6,
# so the smallest two-sided Wilcoxon p attainable is 2/choose(12,6) = 0.00216
# even under complete rank separation. After BH across this 900-test family
# nothing can fall below padj ~0.0046, so "***" (padj<0.001) is impossible
# here by construction and only */** ever appear -- that is a property of the
# sample size, not evidence that effects are weaker than in the pooled figure
# (09), where n is doubled. There are also only 19 distinct raw p-values
# across all 900 tests. Read the colour (Cohen's d) as the primary signal at
# this resolution; the stars are coarse.

# EFFECT SIZE: fill is contrib_diff = mean(test) - mean(control) in
# species-adjusted tAge units, NOT Cohen's d. d divides by the pooled
# within-group SD, which makes the tight 6-sample temporal design produce
# values ~10x the meta-analysis for near-identical real differences, and is
# not comparable across analyses. See 05_consolidate_partial_scores.R header.

source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

TIME_LEVELS <- c("4_days", "10_days", "20_days")
CT_LEVELS <- c("Fibroblast", "Keratinocyte", "Melanocyte")

df <- read.csv(file.path(RERUN_DIR, "partial_tage_ALL.csv"), check.names = FALSE)
bt <- df[df$analysis == "temporal_bytimepoint", ]
bt$pathway_short <- gsub("^HALLMARK ", "", bt$pathway)
# ---------------------------------------------------------------------------
# REPRESENTATION FILTER (added 2026-08-17)
# A gene set can only support a pathway-level claim if the clock's prediction
# for it is spread across many genes rather than carried by one or two. Sets are
# gated in exploratory/14_pathway_representation.py: a set is plotted only if its
# five largest-contributing genes carry <=65% of its total contribution under BOTH
# models (17 of 50 sets). Row labels carry the evidence: clock genes in the set,
# the % with a non-zero coefficient under each model (zero-coefficient genes
# contribute exactly nothing, ever), and that worst-case top-5 share.
# ---------------------------------------------------------------------------
rep_csv <- file.path(RERUN_DIR, "pathway_representation.csv")
if (!file.exists(rep_csv)) {
  stop("pathway_representation.csv not found -- run exploratory/14_pathway_representation.py first.")
}
rep <- read.csv(rep_csv, check.names = FALSE)
rep <- rep[rep$tier == "INTERPRETABLE", ]
rep$label_full <- sprintf("%s  (%d; %.0f/%.0f%% nz; top5 %.0f%%)",
                          gsub("^HALLMARK ", "", rep$pathway),
                          rep$n_clock_scaled,
                          100 * rep$n_nonzero_scaled / rep$n_clock_scaled,
                          100 * rep$n_nonzero_yugene / rep$n_clock_yugene,
                          100 * rep$top5_max)
cat(sprintf("Gate (top5 <= 65%% both models): keeping %d of 50 gene sets\n", nrow(rep)))


sig_symbol <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))

plot_df <- bt %>%
  filter(pathway %in% rep$pathway) %>%
  mutate(
    pathway_short = rep$label_full[match(pathway, rep$pathway)],
    pathway_short = factor(pathway_short, levels = sort(unique(pathway_short), decreasing = TRUE)),
    cell_type = factor(cell_type, levels = CT_LEVELS),
    timepoint = factor(timepoint, levels = TIME_LEVELS,
                       labels = c("4d", "10d", "20d")),
    model_label = factor(ifelse(model == "scaled", "Scaled", "YuGene"),
                          levels = c("Scaled", "YuGene")),
    stars = sig_symbol(p_adj)
  )

cap <- ceiling(quantile(abs(plot_df$contrib_diff), 0.97, na.rm = TRUE))

p <- ggplot(plot_df, aes(x = model_label, y = pathway_short, fill = contrib_diff)) +
  geom_tile(colour = "grey80") +
  geom_text(aes(label = stars), size = 3) +
  facet_grid(cols = vars(cell_type, timepoint)) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-cap, cap), oob = scales::squish, name = "tAge units") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.ticks.y = element_line(colour = "grey40"),
    axis.title = element_blank(),
    strip.text = element_text(size = 13, face = "plain"),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    panel.grid = element_blank(),
    panel.spacing.x = unit(0.25, "lines")
  )

ggsave(file.path(RERUN_DIR, "figure_temporal_pathway_heatmap_bytimepoint.png"), p,
       width = 16, height = 7.5, dpi = 300)
cat(sprintf("Saved -> %s\n", file.path(RERUN_DIR, "figure_temporal_pathway_heatmap_bytimepoint.png")))

n_sig <- plot_df %>% group_by(pathway_short) %>% summarise(n_sig = sum(p_adj < 0.05))
cat("\nDistribution of n significant tests (of 18: 3 cell types x 3 timepoints x 2 models):\n")
print(table(n_sig$n_sig))
