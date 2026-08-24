# 19_pathway_heatmap_within_study.R
#
# Replaces 08_pathway_full_heatmap_both_models.R as the figure for 2.1.5.2.
#
# WHY. Script 08 fills each tile with contrib_diff from partial_tage_ALL.csv,
# which is a POOLED test-minus-control difference. meta_analysis/13 showed pooled
# meta-analysis contrasts are confounded by control baseline, and
# exploratory/18 showed the pathway contributions inherit it: 15-16 of 50 sets
# change sign in CICQ and 8-11 in RS once each sample is compared to its own
# study's controls. The old tiles are therefore the wrong numbers, and in several
# cells the wrong colour.
#
# WHAT CHANGES. Fill is contrib_diff_within_study from
# partial_tage_within_study.csv. Stars are the within-study permutation FDR.
# Cells where the set is reportable - significant on BOTH models AND the two
# models agree in sign - are outlined in black; when the normalisations disagree
# in direction the effect is not identifiable and is deliberately left unmarked
# even if both are individually significant (3 such cells: CICQ XENOBIOTIC
# METABOLISM, SSCQ EMT, SIPS ESTROGEN RESPONSE LATE).
#
# UNCHANGED. Visual grammar is kept so the figure remains comparable with the
# temporal heatmaps (09/10): geom_tile, alphabetical rows, one facet per
# condition, both models side by side, fill in species-adjusted tAge units rather
# than Cohen's d, and row labels carrying the representation evidence.
# 06_pathway_effect_heatmap.R and 07_pathway_divergence_meta_conditions.R are
# exploratory-only and their meta rows are superseded by this script; the temporal
# heatmaps are NOT affected, the time course being a single study.

source("R/config.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

w <- read.csv(file.path(RERUN_DIR, "partial_tage_within_study.csv"), check.names = FALSE)
rep <- read.csv(file.path(RERUN_DIR, "pathway_representation.csv"), check.names = FALSE)
rep <- rep[rep$tier == "INTERPRETABLE", ]
rep$label_full <- sprintf("%s  (%d; %.0f/%.0f%% nz; top5 %.0f%%)",
                          gsub("^HALLMARK ", "", rep$pathway), rep$n_clock_scaled,
                          100 * rep$n_nonzero_scaled / rep$n_clock_scaled,
                          100 * rep$n_nonzero_yugene / rep$n_clock_yugene,
                          100 * rep$top5_max)
cat(sprintf("Gate (top5 <= 65%% both models): %d of 50 gene sets\n", nrow(rep)))

sig_symbol <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                          ifelse(p < 0.05, "*", "")))

# reportable = significant on both models AND same sign on both
flag <- w %>% filter(pathway %in% rep$pathway) %>%
  group_by(label, pathway) %>%
  summarise(both_sig = all(p_perm_adj < 0.05),
            same_sign = prod(sign(contrib_diff_within_study)) > 0,
            .groups = "drop") %>%
  mutate(reportable = both_sig & same_sign)
cat(sprintf("reportable cells (both models sig + sign agree): %d of %d\n",
            sum(flag$reportable), nrow(flag)))
cat(sprintf("  dropped for between-model sign disagreement: %d\n",
            sum(flag$both_sig & !flag$same_sign)))

plot_df <- w %>% filter(pathway %in% rep$pathway) %>%
  left_join(flag[, c("label", "pathway", "reportable")], by = c("label", "pathway")) %>%
  mutate(pathway_short = rep$label_full[match(pathway, rep$pathway)],
         pathway_short = factor(pathway_short,
                                levels = sort(unique(pathway_short), decreasing = TRUE)),
         label = factor(label, levels = c("CICQ", "SSCQ", "RS", "SIPS", "OIS")),
         model_label = factor(ifelse(model == "scaled", "Scaled", "YuGene"),
                              levels = c("Scaled", "YuGene")),
         stars = sig_symbol(p_perm_adj))

cap <- ceiling(quantile(abs(plot_df$contrib_diff_within_study), 0.97, na.rm = TRUE))

p <- ggplot(plot_df, aes(x = model_label, y = pathway_short,
                         fill = contrib_diff_within_study)) +
  geom_tile(colour = "grey80") +
  geom_tile(data = subset(plot_df, reportable), colour = "black", linewidth = 0.7,
            fill = NA) +
  geom_text(aes(label = stars), size = 4) +
  facet_grid(~label) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-cap, cap), oob = scales::squish,
                       name = "tAge units") +
  theme_minimal(base_size = 16) +
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.ticks.y = element_line(colour = "grey40"),
        axis.title = element_blank(),
        strip.text = element_text(size = 16, face = "plain"),
        legend.text = element_text(size = 13), legend.title = element_text(size = 15),
        panel.grid = element_blank(), panel.spacing = unit(0.4, "lines"))

out <- file.path(RERUN_DIR, "figure_pathway_heatmap_within_study.png")
ggsave(out, p, width = 12, height = 7.5, dpi = 300)
cat(sprintf("Saved -> %s\n", out))

old <- file.path(RERUN_DIR, "figure_pathway_heatmap_all_both_models.png")
if (file.exists(old)) {
  file.rename(old, file.path(RERUN_DIR,
              "SUPERSEDED_figure_pathway_heatmap_all_both_models.png"))
  cat("Renamed old pooled heatmap -> SUPERSEDED_...png\n")
}
