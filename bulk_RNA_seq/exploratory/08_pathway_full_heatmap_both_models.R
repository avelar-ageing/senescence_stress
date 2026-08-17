# 08_pathway_full_heatmap_both_models.R
#
# Figure for section 2.1.5.2 (Pathway-level Transcriptomic Age Differences):
# comprehensive partial-tAge heatmap, meta-analysis only, ALL 50 pathways
# (every pathway reaches significance in at least one of the 5 conditions x
# 2 models = 10 tests, so no pathway is dropped), both EN models shown
# side-by-side per condition so model-(dis)agreement is visible directly in
# the figure rather than only in partial_tage_ALL.csv. Complements
# 06_pathway_effect_heatmap.R (yugene-only, recurrence-filtered cross-
# analysis view) and 07_pathway_divergence_meta_conditions.R (yugene-only,
# meta-analysis divergence ranking) rather than replacing them.
#
# Plain geom_tile grid, not pheatmap: with clustering off (rows are
# alphabetical, a lookup reference -- see script history/discussion) and no
# dendrogram, pheatmap wasn't buying anything but its annotation strips,
# which facet_wrap reproduces directly and more controllably.

source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

df <- read.csv(file.path(RERUN_DIR, "partial_tage_ALL.csv"), check.names = FALSE)
meta <- df[df$analysis == "meta_analysis", ]
meta$pathway_short <- gsub("^HALLMARK ", "", meta$pathway)
# ---------------------------------------------------------------------------
# REPRESENTATION FILTER (added 2026-08-17)
# A gene set can only support a pathway-level claim if the clock's prediction
# for it is spread across many genes rather than carried by one or two. Sets are
# tiered by effective number of contributing genes (eff_n = 1/sum(share^2)) in
# exploratory/14_pathway_representation.py; only WELL (eff_n>=14) and MODERATE
# (>=9) sets are plotted. Row labels carry the evidence: clock genes in the set,
# then the % with a non-zero coefficient under each model (scaled/yugene) --
# genes with a zero coefficient contribute exactly nothing, ever.
# ---------------------------------------------------------------------------
rep_csv <- file.path(RERUN_DIR, "pathway_representation.csv")
if (!file.exists(rep_csv)) {
  stop("pathway_representation.csv not found -- run exploratory/14_pathway_representation.py first.")
}
rep <- read.csv(rep_csv, check.names = FALSE)
rep <- rep[rep$tier %in% c("WELL", "MODERATE"), ]
rep$label_full <- sprintf("%s  (%d; %.0f%%/%.0f%%)",
                          gsub("^HALLMARK ", "", rep$pathway),
                          rep$n_clock_scaled,
                          100 * rep$n_nonzero_scaled / rep$n_clock_scaled,
                          100 * rep$n_nonzero_yugene / rep$n_clock_yugene)
cat(sprintf("Representation filter: keeping %d of 50 gene sets (WELL/MODERATE)\n", nrow(rep)))


sig_symbol <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))

plot_df <- meta %>%
  filter(pathway %in% rep$pathway) %>%
  mutate(
    pathway_short = rep$label_full[match(pathway, rep$pathway)],
    pathway_short = factor(pathway_short, levels = sort(unique(pathway_short), decreasing = TRUE)),
    label = factor(label, levels = c("CICQ", "SSCQ", "RS", "SIPS", "OIS")),
    model_label = ifelse(model == "scaled", "Scaled", "YuGene"),
    model_label = factor(model_label, levels = c("Scaled", "YuGene")),
    stars = sig_symbol(p_adj)
  )

cap <- min(6, ceiling(quantile(abs(plot_df$cohens_d), 0.97, na.rm = TRUE)))

p <- ggplot(plot_df, aes(x = model_label, y = pathway_short, fill = cohens_d)) +
  geom_tile(colour = "grey80") +
  geom_text(aes(label = stars), size = 4) +
  facet_grid(~label) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-cap, cap), oob = scales::squish, name = "Cohen's d") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.ticks.y = element_line(colour = "grey40"),
    axis.title = element_blank(),
    strip.text = element_text(size = 16, face = "plain"),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    panel.grid = element_blank(),
    panel.spacing = unit(0.4, "lines")
  )

ggsave(file.path(RERUN_DIR, "figure_pathway_heatmap_all_both_models.png"), p, width = 12, height = 7.5, dpi = 300)
cat(sprintf("Saved -> %s\n", file.path(RERUN_DIR, "figure_pathway_heatmap_all_both_models.png")))

n_sig <- meta %>% group_by(pathway_short) %>% summarise(n_sig = sum(p_adj < 0.05))
cat("\nDistribution of n significant tests (of 10: 5 conditions x 2 models), all 50 pathways:\n")
print(table(n_sig$n_sig))
